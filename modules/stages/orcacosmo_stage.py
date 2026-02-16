import os
import shutil
import subprocess
import json
import time

import pandas as pd

from modules.stages.base_stage import BaseStage
from modules.parsers.cosmo.orca_log_parser import OrcaLogParser
from modules.parsers.cosmo.orca_cpcm_parser import OrcaCpcmParser
from modules.parsers.cosmo.orca_cpcm_corr_parser import OrcaCpcmCorrParser

from modules.orcacosmo_reconstructors.orcacosmo_orchestrator import OrcaCosmoOrchestrator
from modules.utils.atomic_write import AtomicWriter


class OrcacosmoStage(BaseStage):
    """
    ORCA COSMO stage.

    Input:
        stage_input = energies.json from Optimisation (or Pruning)
        Each entry contains:
          - lookup_id
          - inchi_key
          - conf_num (optional)
          - xyz_path
          - energy
          - smiles
          - provenance
          - optimisation_history

    Output (canonical stage_output):
        orcacosmo_summary.json

    Auxiliary outputs:
        orcacosmo_summary.csv
        item_to_lookup_mapping.json
        raw_outputs/<lookup_id>.{log,cpcm,cpcm_corr}
        parsed_outputs/<lookup_id>.{log,cpcm,cpcm_corr}.json
        orcacosmo_outputs/<lookup_id>.orcacosmo
    """

    # ------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------
    def execute(self):
        self.log_header("Starting ORCA COSMO Stage")

        # Canonical stage output
        self.set_stage_output("orcacosmo_summary.json")

        self._load_stage_config()
        self._load_cpcm_data()
        self._prepare_directories()
        self._load_xyzs_and_entries()

        lookup_ids = self._discover_lookup_ids()
        self.set_items(lookup_ids)

        self.successful_outputs = []
        self.item_index = 0
        self.item_to_lookup = {}
        self.stats = {
            "successful": 0,
            "failed": 0,
            "fallback_used": 0,
            "missing_xyz": 0,
        }

        stage_start = time.perf_counter()

        for lookup_id in list(self.job.pending_items):
            xyz_path = os.path.join(self.inputs_dir, f"{lookup_id}.xyz")

            if not os.path.exists(xyz_path):
                self._skip_missing_xyz(lookup_id)
                continue

            item_num = self._assign_item_number(lookup_id)

            try:
                self._process_single(lookup_id, xyz_path, item_num)
                self.update_progress(lookup_id)
                self.stats["successful"] += 1
            except Exception as e:
                self._handle_failure(lookup_id, e)
                self.stats["failed"] += 1

        stage_end = time.perf_counter()
        self._write_summary(stage_start, stage_end)
        self._log_warning_summary()

        self.job.mark_complete()
        self.log_header("ORCA COSMO Stage Complete")

    # ------------------------------------------------------------
    # Stage config
    # ------------------------------------------------------------
    def _load_stage_config(self):
        cfg = self.config
        if not cfg:
            self.fail("OrcacosmoStage: missing config in job.config")

        orca_cfg = cfg.get("orca_6_1_1") or cfg.get("orca", {}) or cfg.get("orca_6_0_0")
        self.orca_command = orca_cfg.get("executable")
        self.orca_home = orca_cfg.get("home")

        if not self.orca_command:
            self.fail("Missing ORCA executable in config")

        const_cfg = cfg.get("constant_files", {})
        self.chemistry_dir = const_cfg.get("chemistry_dir")
        self.metadata_dir = const_cfg.get("metadata_dir")

        if not self.chemistry_dir:
            self.fail("Missing chemistry_dir in config['constant_files']")
        if not self.metadata_dir:
            self.fail("Missing metadata_dir in config['constant_files']")

        self.enable_fallback = cfg.get("enable_tzvp_fallback", False)
        self.strict_mode = bool(cfg.get("cosmo", {}).get("strict", False))

        # Determine max parallel processes for ORCA (%pal nprocs)
        max_procs = None

        resources = getattr(self, "resources", None)
        if resources is not None and getattr(resources, "cpus", None):
            try:
                max_procs = int(resources.cpus)
            except Exception:
                max_procs = None

        if max_procs is None:
            params = getattr(self, "parameters", {}) or {}
            try:
                max_procs = int(params.get("resources", {}).get("cpus", 0) or 0)
            except Exception:
                max_procs = None

        if not max_procs or max_procs <= 0:
            max_procs = os.cpu_count() or 1

        physical_cores = os.cpu_count() or 1
        self.max_procs = max(1, min(max_procs, physical_cores))

        self.log(
            f"[INFO] ORCA parallel nprocs set to {self.max_procs} "
            f"(requested={max_procs}, physical={physical_cores})"
        )

        # ORCA version (best-effort)
        self.orca_version = "unknown"
        try:
            out = subprocess.check_output([self.orca_command, "-h"], stderr=subprocess.STDOUT)
            self.orca_version = out.decode(errors="ignore").splitlines()[0].strip()
        except Exception:
            pass

    # ------------------------------------------------------------
    # CPCM radii
    # ------------------------------------------------------------
    def _load_cpcm_data(self):
        cpcm_file = os.path.join(self.chemistry_dir, "cpcm_radii.json")
        if not os.path.exists(cpcm_file):
            self.fail(f"Missing CPCM radii JSON: {cpcm_file}")

        with open(cpcm_file) as f:
            cpcm = json.load(f)

        self.cpcm_radii = cpcm["radii"]
        self.cpcm_cut_area = cpcm["cut_area"]
        self.cpcm_file = cpcm_file

        self.log(f"[INFO] Loaded CPCM radii from {cpcm_file}")

    # ------------------------------------------------------------
    # Directory setup
    # ------------------------------------------------------------
    def _prepare_directories(self):
        self.orca_inputs_dir = os.path.join(self.outputs_dir, "orca_inputs")
        self.parsed_outputs_dir = os.path.join(self.outputs_dir, "parsed_outputs")
        self.raw_outputs_dir = os.path.join(self.outputs_dir, "raw_outputs")
        self.tmp_exec = os.path.join(os.path.dirname(self.outputs_dir), "tmp_exec")
        self.orcacosmo_outputs_dir = os.path.join(self.outputs_dir, "orcacosmo_outputs")

        for d in [
            self.orca_inputs_dir,
            self.parsed_outputs_dir,
            self.raw_outputs_dir,
            self.tmp_exec,
            self.orcacosmo_outputs_dir,
        ]:
            os.makedirs(d, exist_ok=True)

        self.log("[INFO] Output directories prepared")

    # ------------------------------------------------------------
    # Load XYZs + full entries from stage_input energies.json
    # ------------------------------------------------------------
    def _load_xyzs_and_entries(self):
        energies_file = self.require_file(
            self.get_stage_input(),
            "stage_input energies.json"
        )
        self.log(f"Stage input: {energies_file}")

        # Copy energies.json into inputs for provenance
        inputs_energies = os.path.join(self.inputs_dir, "energies.json")
        os.makedirs(self.inputs_dir, exist_ok=True)
        shutil.copy(energies_file, inputs_energies)

        with open(energies_file) as f:
            entries = json.load(f)

        # Preserve full entries exactly as given
        self.entries = entries
        self.entry_map = {e["lookup_id"]: e for e in entries}

        # Copy XYZs into inputs/
        for entry in entries:
            lookup_id = entry.get("lookup_id")
            xyz_src = entry.get("xyz_path")

            if lookup_id and xyz_src and os.path.exists(xyz_src):
                dst = os.path.join(self.inputs_dir, f"{lookup_id}.xyz")
                shutil.copy(xyz_src, dst)

        self.log("[INFO] XYZs + optimisation entries loaded from stage_input")

    # ------------------------------------------------------------
    # Lookup discovery
    # ------------------------------------------------------------
    def _discover_lookup_ids(self):
        ids = [
            os.path.splitext(f)[0]
            for f in os.listdir(self.inputs_dir)
            if f.endswith(".xyz")
        ]
        ids.sort()
        self.log(f"Discovered {len(ids)} XYZ structures")
        return ids

    # ------------------------------------------------------------
    # Helpers for execute()
    # ------------------------------------------------------------
    def _skip_missing_xyz(self, lookup_id):
        self.log(f"[WARNING] Missing XYZ for {lookup_id}, skipping.")
        self.update_progress(lookup_id, success=False)
        self.stats["missing_xyz"] += 1

    def _assign_item_number(self, lookup_id):
        item_num = self.item_index
        self.item_to_lookup[item_num] = lookup_id
        self.item_index += 1
        return item_num

    def _handle_failure(self, lookup_id, e):
        self.log(f"[ERROR] {lookup_id}: {e}")
        self.update_progress(lookup_id, success=False)
        if self.strict_mode:
            self.fail(f"Strict mode enabled; aborting due to failure on {lookup_id}")

    # ------------------------------------------------------------
    # Main processing
    # ------------------------------------------------------------
    def _process_single(self, lookup_id, xyz_path, item_num):
        total = len(self.job.items)

        self.log_item(item_num, total, f"Starting ORCA run (lookup_id={lookup_id})")

        workdir = self._prepare_workdir(item_num)
        local_xyz = self._copy_xyz_to_workdir(lookup_id, xyz_path, workdir)

        start = time.perf_counter()

        # Run ORCA (with fallback)
        method_used, fallback_triggered = self._run_orca_with_fallback(
            lookup_id, item_num, workdir
        )

        elapsed = time.perf_counter() - start
        self.log_item_step(item_num, f"ORCA finished in {elapsed:.1f} s")

        # Copy raw outputs
        self._copy_raw_outputs(lookup_id, item_num, workdir)
        self.log_item_step(item_num, "Copied raw outputs")

        # Parse
        log_json = self._parse_log_json(lookup_id)
        cpcm_json = self._parse_cpcm_json(lookup_id)
        cpcm_corr_json = self._parse_cpcm_corr_json(lookup_id)
        self.log_item_step(item_num, "Parsed log/cpcm/cpcm_corr")

        # Write parser JSONs
        self._write_parser_json(lookup_id, "log", log_json)
        self._write_parser_json(lookup_id, "cpcm", cpcm_json)
        self._write_parser_json(lookup_id, "cpcm_corr", cpcm_corr_json)

        # Build bundle
        bundle = self._build_bundle(
            lookup_id,
            method_used,
            fallback_triggered,
            log_json,
            cpcm_json,
            cpcm_corr_json,
        )

        bundle_path = os.path.join(self.parsed_outputs_dir, f"{lookup_id}_bundle.json")
        with AtomicWriter(bundle_path) as f:
            json.dump(bundle, f, indent=2)

        # Orchestrate
        orchestrator = OrcaCosmoOrchestrator(bundle)
        orcacosmo_text = orchestrator.reconstruct()

        orcacosmo_path = os.path.join(self.orcacosmo_outputs_dir, f"{lookup_id}.orcacosmo")
        with open(orcacosmo_path, "w") as f:
            f.write(orcacosmo_text)

        self.log_item_step(item_num, f"Wrote .orcacosmo → {orcacosmo_path}")

        # Record success (no parsed JSON stored)
        self._record_success(
            lookup_id=lookup_id,
            item_num=item_num,
            method_used=method_used,
            fallback_triggered=fallback_triggered,
            elapsed_seconds=elapsed,
            log_json=None,
            cpcm_json=None,
            cpcm_corr_json=None,
            orcacosmo_path=orcacosmo_path,
            bundle_path=bundle_path,
        )

        # Append lightweight summary entry
        entry = self.entry_map[lookup_id]
        self._append_summary_entry(
            lookup_id=lookup_id,
            inchi_key=entry["inchi_key"],
            orcacosmo_path=orcacosmo_path,
            method_used=method_used,
            fallback_triggered=fallback_triggered,
        )

        self.log_item(item_num, total, f"Completed successfully in {elapsed:.1f} s")

    # ------------------------------------------------------------
    # Workdir helpers
    # ------------------------------------------------------------
    def _prepare_workdir(self, item_num):
        workdir = os.path.join(self.tmp_exec, f"item{item_num:03d}")
        os.makedirs(workdir, exist_ok=True)
        return workdir

    def _copy_xyz_to_workdir(self, lookup_id, xyz_path, workdir):
        local_xyz = os.path.join(workdir, f"{lookup_id}.xyz")
        shutil.copy(xyz_path, local_xyz)
        return local_xyz

    # ------------------------------------------------------------
    # ORCA input writers
    # ------------------------------------------------------------
    def _write_orca_input_tzvpd(self, path, xyz_name, item_num):
        with open(path, "w") as f:
            f.write("%MaxCore 2000\n\n")

            f.write(f'%base "item{item_num:03d}"\n\n')

            f.write("%cpcm\n")
            for Z, r in self.cpcm_radii.items():
                f.write(f"radius[{Z}]  {r}\n")
            f.write(f"cut_area {self.cpcm_cut_area}\n")
            f.write("end\n\n")

            f.write("! CPCM BP86 def2-TZVPD SP\n\n")

            f.write("%elprop\n")
            f.write("  Polar 1\n")
            f.write("  Polaratom 1\n")
            f.write("end\n\n")

            f.write(f"* xyzfile 0 1 {xyz_name}\n")

    def _write_orca_input_tzvp(self, path, xyz_name, item_num):
        with open(path, "w") as f:
            f.write("%MaxCore 2000\n\n")

            f.write(f'%base "item{item_num:03d}"\n\n')

            f.write("%cpcm\n")
            for Z, r in self.cpcm_radii.items():
                f.write(f"radius[{Z}]  {r}\n")
            f.write(f"cut_area {self.cpcm_cut_area}\n")
            f.write("end\n\n")

            f.write("! CPCM BP86 def2-TZVP SP\n\n")

            f.write("%elprop\n")
            f.write("  Polar 1\n")
            f.write("  Polaratom 1\n")
            f.write("end\n\n")

            f.write(f"* xyzfile 0 1 {xyz_name}\n")

    # ------------------------------------------------------------
    # Run ORCA safely
    # ------------------------------------------------------------
    def _run_orca(self, workdir, inp_file, log_path):
        env = os.environ.copy()
        orca_dir = os.path.dirname(self.orca_command)
        env["LD_LIBRARY_PATH"] = f"{orca_dir}:{env.get('LD_LIBRARY_PATH','')}"

        # Set OpenMP threads for parallelization
        env["OMP_NUM_THREADS"] = str(self.max_procs)

        # Disable MPI to avoid slot allocation issues
        for key in list(env.keys()):
            if 'OMPI' in key or 'MPI' in key or 'PRTE' in key:
                del env[key]

        self.log(f"[INFO] Running ORCA with OMP_NUM_THREADS={self.max_procs}")

        with open(log_path, "w") as f:
            subprocess.run(
                [self.orca_command, inp_file],
                cwd=workdir,
                stdout=f,
                stderr=subprocess.STDOUT,
                check=False,
                env=env,
            )

    def _run_orca_with_fallback(self, lookup_id, item_num, workdir):
        base = f"item{item_num:03d}"
        inp_tzvpd = f"{base}_tzvpd.inp"
        inp_tzvp = f"{base}_tzvp.inp"
        log = os.path.join(workdir, f"{base}.log")

        # TZVPD
        inp_path = os.path.join(workdir, inp_tzvpd)
        self._write_orca_input_tzvpd(inp_path, f"{lookup_id}.xyz", item_num)
        shutil.copy(inp_path, os.path.join(self.orca_inputs_dir, f"{lookup_id}_tzvpd.inp"))
        self.log(f"[INFO] Wrote ORCA TZVPD input: {inp_path}")

        self._run_orca(workdir, inp_path, log)

        cpcm = os.path.join(workdir, f"{base}.cpcm")
        cpcm_corr = os.path.join(workdir, f"{base}.cpcm_corr")

        if self._validate_cpcm(cpcm, cpcm_corr, log):
            return "TZVPD", False

        # Fallback
        if not self.enable_fallback:
            self.fail(f"TZVPD failed for {lookup_id} and fallback is disabled.")

        self._cleanup_failed_run(workdir, base)

        inp_path = os.path.join(workdir, inp_tzvp)
        self._write_orca_input_tzvp(inp_path, f"{lookup_id}.xyz", item_num)
        shutil.copy(inp_path, os.path.join(self.orca_inputs_dir, f"{lookup_id}_tzvp.inp"))
        self.log(f"[INFO] Wrote ORCA TZVP fallback input: {inp_path}")

        self._run_orca(workdir, inp_path, log)

        if not self._validate_cpcm(cpcm, cpcm_corr, log):
            self.fail(f"Both TZVPD and TZVP failed for {lookup_id}")

        return "TZVP", True

    def _cleanup_failed_run(self, workdir, base):
        for ext in ['.log', '.cpcm', '.cpcm_corr', '.gbw', '.property.txt']:
            p = os.path.join(workdir, f"{base}{ext}")
            if os.path.exists(p):
                os.remove(p)

    # ------------------------------------------------------------
    # Raw output copying
    # ------------------------------------------------------------
    def _copy_raw_outputs(self, lookup_id, item_num, workdir):
        base_name = f"item{item_num:03d}"

        for ext in ["log", "cpcm", "cpcm_corr"]:
            src = os.path.join(workdir, f"{base_name}.{ext}")
            dst = os.path.join(self.raw_outputs_dir, f"{lookup_id}.{ext}")

            if not os.path.exists(src):
                self.log(f"[WARNING] Missing raw output: {src}")
                continue

            shutil.copy(src, dst)
            self.log(f"[COPY] {src} -> {dst}")

    # ------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------
    def _parse_log_json(self, lookup_id):
        path = os.path.join(self.raw_outputs_dir, f"{lookup_id}.log")
        if not os.path.exists(path):
            self.fail(f"Cannot parse log - file not found: {path}")

        self.log(f"[PARSE] Log: {path}")
        return OrcaLogParser(path).parse()

    def _parse_cpcm_json(self, lookup_id):
        path = os.path.join(self.raw_outputs_dir, f"{lookup_id}.cpcm")
        if not os.path.exists(path):
            self.fail(f"Cannot parse cpcm - file not found: {path}")

        self.log(f"[PARSE] CPCM: {path}")
        return OrcaCpcmParser(path).parse()

    def _parse_cpcm_corr_json(self, lookup_id):
        path = os.path.join(self.raw_outputs_dir, f"{lookup_id}.cpcm_corr")
        if not os.path.exists(path):
            self.fail(f"Cannot parse cpcm_corr - file not found: {path}")

        self.log(f"[PARSE] CPCM_CORR: {path}")
        return OrcaCpcmCorrParser(path).parse()

    # ------------------------------------------------------------
    # Parser JSON writers (for orchestrator)
    # ------------------------------------------------------------
    def _write_parser_json(self, lookup_id, kind, data):
        out = os.path.join(self.parsed_outputs_dir, f"{lookup_id}.{kind}.json")
        with AtomicWriter(out) as f:
            json.dump(data, f, indent=2)
        self.log(f"[WRITE] Parser JSON: {out}")

    # ------------------------------------------------------------
    # Build bundle for orchestrator (with provenance)
    # ------------------------------------------------------------
    def _build_bundle(self, lookup_id, method_used, fallback_triggered,
                      log_json, cpcm_json, cpcm_corr_json):
        entry = self.entry_map.get(lookup_id)
        if not entry:
            self.fail(f"No optimisation entry found for lookup_id: {lookup_id}")

        inchi_key = entry.get("inchi_key")
        smiles = entry.get("smiles")

        return {
            "meta": {
                "lookup_id": lookup_id,
                "inchi_key": inchi_key,
                "smiles": smiles,
                "method_used": method_used,
                "fallback_triggered": fallback_triggered,
                "orca_input_path": os.path.join(
                    self.orca_inputs_dir,
                    f"{lookup_id}_{method_used.lower()}.inp"
                ),
                "cpcm_radii_source": self.cpcm_file,
                "orca_version": self.orca_version,
            },
            "paths": {
                "log": os.path.join(self.parsed_outputs_dir, f"{lookup_id}.log.json"),
                "cpcm": os.path.join(self.parsed_outputs_dir, f"{lookup_id}.cpcm.json"),
                "cpcm_corr": os.path.join(self.parsed_outputs_dir, f"{lookup_id}.cpcm_corr.json"),
                "xyz": os.path.join(self.inputs_dir, f"{lookup_id}.xyz"),
            },
            "optimisation_entry": entry,
        }

    # ------------------------------------------------------------
    # Record success (preserve original entry + add cosmo block)
    # ------------------------------------------------------------
    def _record_success(
        self,
        lookup_id,
        item_num,
        method_used,
        fallback_triggered,
        elapsed_seconds,
        log_json,
        cpcm_json,
        cpcm_corr_json,
        orcacosmo_path,
        bundle_path,
    ):
        original = self.entry_map.get(lookup_id, {}).copy()

        cosmo_block = {
            "item_number": item_num,
            "method_used": method_used,
            "fallback_triggered": fallback_triggered,
            "elapsed_seconds": elapsed_seconds,
            "orcacosmo_path": orcacosmo_path,
            "bundle_path": bundle_path,
            "raw_outputs": {
                "log": os.path.join(self.raw_outputs_dir, f"{lookup_id}.log"),
                "cpcm": os.path.join(self.raw_outputs_dir, f"{lookup_id}.cpcm"),
                "cpcm_corr": os.path.join(self.raw_outputs_dir, f"{lookup_id}.cpcm_corr"),
            },
            "provenance": {
                "orca_version": self.orca_version,
                "cpcm_radii_source": self.cpcm_file,
                "engine_command_base": "orca",
                "item_number": item_num,
                "last_modified": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            },
        }

        original["cosmo"] = cosmo_block
        self.successful_outputs.append(original)

    # ------------------------------------------------------------
    # Summary writer (big JSON + CSV + mapping)
    # ------------------------------------------------------------
    def _write_summary(self, stage_start, stage_end):
        results_path = self.get_stage_output()
        mapping_path = os.path.join(self.outputs_dir, "item_to_lookup_mapping.json")
        summary_csv = os.path.join(self.outputs_dir, "orcacosmo_summary.csv")

        # Mapping
        with AtomicWriter(mapping_path) as f:
            json.dump(self.item_to_lookup, f, indent=2)

        # CSV
        rows = []
        for entry in self.successful_outputs:
            cosmo = entry["cosmo"]
            rows.append({
                "lookup_id": entry["lookup_id"],
                "inchi_key": entry["inchi_key"],
                "smiles": entry["smiles"],
                "energy": entry["energy"],
                "method_used": cosmo["method_used"],
                "fallback_triggered": cosmo["fallback_triggered"],
                "elapsed_seconds": cosmo["elapsed_seconds"],
                "orcacosmo_path": cosmo["orcacosmo_path"],
                "item_number": cosmo["item_number"],
            })

        if rows:
            pd.DataFrame(rows).to_csv(summary_csv, index=False)

        self.log(
            f"Wrote ORCA COSMO results JSON with {len(self.successful_outputs)} entries: {results_path}"
        )
        self.log(
            f"Wrote item-to-lookup mapping with {len(self.item_to_lookup)} entries: {mapping_path}"
        )
        if rows:
            self.log(f"Wrote ORCA COSMO summary CSV: {summary_csv}")

        self.log(
            f"Stage elapsed time: {stage_end - stage_start:.2f} s "
            f"(successful={self.stats['successful']}, failed={self.stats['failed']})"
        )

    # ------------------------------------------------------------
    # Append lightweight summary entry
    # ------------------------------------------------------------
    def _append_summary_entry(self, lookup_id, inchi_key, orcacosmo_path, method_used, fallback_triggered):
        """Append a lightweight summary entry as each item finishes."""
        summary_path = os.path.join(self.outputs_dir, "orcacosmo_summary.json")

        # Load existing summary if present
        if os.path.exists(summary_path):
            with open(summary_path) as f:
                data = json.load(f)
        else:
            data = []

        data.append({
            "lookup_id": lookup_id,
            "inchi_key": inchi_key,
            "orcacosmo_path": orcacosmo_path,
            "method_used": method_used,
            "fallback_triggered": fallback_triggered,
        })

        with AtomicWriter(summary_path) as f:
            json.dump(data, f, indent=2)

    # ------------------------------------------------------------
    # CPCM validation
    # ------------------------------------------------------------
    def _validate_cpcm(self, cpcm_path, cpcm_corr_path, log_path):
        if not os.path.exists(cpcm_path):
            self.log(f"[VALIDATION] Missing cpcm file: {cpcm_path}")
            return False
        if not os.path.exists(cpcm_corr_path):
            self.log(f"[VALIDATION] Missing cpcm_corr file: {cpcm_corr_path}")
            return False
        if not os.path.exists(log_path):
            self.log(f"[VALIDATION] Missing log file: {log_path}")
            return False

        cpcm_size = os.path.getsize(cpcm_path)
        cpcm_corr_size = os.path.getsize(cpcm_corr_path)
        log_size = os.path.getsize(log_path)

        if cpcm_size < 100:
            self.log(f"[VALIDATION] cpcm file too small: {cpcm_size} bytes")
            return False
        if cpcm_corr_size < 100:
            self.log(f"[VALIDATION] cpcm_corr file too small: {cpcm_corr_size} bytes")
            return False
        if log_size < 1000:
            self.log(f"[VALIDATION] log file too small: {log_size} bytes")
            return False

        self.log(f"[VALIDATION] Basic checks passed for {os.path.basename(cpcm_path)}")
        return True

    # ------------------------------------------------------------
    # Warning summary
    # ------------------------------------------------------------
    def _log_warning_summary(self):
        self.log_header("ORCA COSMO Summary")
        self.log(f"Successful: {self.stats['successful']}")
        self.log(f"Failed: {self.stats['failed']}")
        self.log(f"Fallback used: {self.stats['fallback_used']}")
        self.log(f"Missing XYZ: {self.stats['missing_xyz']}")
        self.log_header("End ORCA COSMO Summary")

    # ------------------------------------------------------------
    # Logging helpers
    # ------------------------------------------------------------
    def log_item(self, item_num, total, message):
        """Top‑level item log, includes item counter."""
        self.log(f"[item{item_num:03d} {item_num+1}/{total}] {message}")

    def log_item_step(self, item_num, message):
        """Indented step log for readability."""
        self.log(f"    → {message}")
