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

    Responsibilities:
      - Load conformer XYZs from an energies.json-like summary_file
      - For each conformer:
          - Run ORCA CPCM (TZVPD, with optional TZVP fallback)
          - Parse LOG, CPCM, CPCM_CORR
          - Build a bundle for the OrcaCosmoOrchestrator
          - Write .orcacosmo file
      - Produce:
          - raw_outputs/<lookup_id>.{log,cpcm,cpcm_corr}
          - orca_inputs/<lookup_id>_{tzvpd|tzvp}.inp
          - orcacosmo_outputs/<lookup_id>.orcacosmo
          - orcacosmo_results.json (all results in one JSON)
          - orcacosmo_summary.csv (human-readable summary)
          - item_to_lookup_mapping.json
    """

    # ------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------
    def execute(self):
        self.log_header("Starting ORCA COSMO Stage")

        self._load_stage_config()
        self._load_cpcm_data()
        self._prepare_directories()
        self._load_xyzs_from_summary()

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
    # XYZ loading from summary_file (energies-like JSON)
    # ------------------------------------------------------------
    def _load_xyzs_from_summary(self):
        summary_file = self.parameters.get("summary_file")
        if not summary_file:
            self.fail("ORCA COSMO stage requires summary_file")

        if not os.path.exists(summary_file):
            self.fail(f"summary_file does not exist: {summary_file}")

        with open(summary_file) as f:
            entries = json.load(f)

        self.energies_entries = entries

        for entry in entries:
            lookup_id = entry.get("lookup_id")
            xyz_src = entry.get("xyz_path")

            if lookup_id and xyz_src and os.path.exists(xyz_src):
                dst = os.path.join(self.inputs_dir, f"{lookup_id}.xyz")
                shutil.copy(xyz_src, dst)

        self.log("[INFO] XYZs loaded from summary_file")

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
        self.log(f"Processing {lookup_id} as item{item_num:03d}")
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
        workdir = self._prepare_workdir(item_num)
        local_xyz = self._copy_xyz_to_workdir(lookup_id, xyz_path, workdir)

        start = time.perf_counter()
        method_used, fallback_triggered = self._run_orca_with_fallback(
            lookup_id, item_num, workdir
        )
        elapsed = time.perf_counter() - start

        self._copy_raw_outputs(lookup_id, item_num, workdir)

        log_json = self._parse_log_json(lookup_id)
        cpcm_json = self._parse_cpcm_json(lookup_id)
        cpcm_corr_json = self._parse_cpcm_corr_json(lookup_id)

        # Keep parser JSONs for orchestrator (internal, not "results")
        self._write_parser_json(lookup_id, "log", log_json)
        self._write_parser_json(lookup_id, "cpcm", cpcm_json)
        self._write_parser_json(lookup_id, "cpcm_corr", cpcm_corr_json)

        bundle = self._build_bundle(lookup_id, method_used, fallback_triggered)

        bundle_path = os.path.join(self.parsed_outputs_dir, f"{lookup_id}_bundle.json")
        with AtomicWriter(bundle_path) as f:
            json.dump(bundle, f, indent=2)

        # Orchestrate
        orchestrator = OrcaCosmoOrchestrator(bundle)
        orcacosmo_text = orchestrator.reconstruct()

        orcacosmo_path = os.path.join(self.orcacosmo_outputs_dir, f"{lookup_id}.orcacosmo")
        with open(orcacosmo_path, "w") as f:
            f.write(orcacosmo_text)
        self.log(f"[WRITE] .orcacosmo written to {orcacosmo_path}")

        if fallback_triggered:
            self.stats["fallback_used"] += 1

        # Record success into in-memory results
        self._record_success(
            lookup_id=lookup_id,
            item_num=item_num,
            method_used=method_used,
            fallback_triggered=fallback_triggered,
            elapsed_seconds=elapsed,
            log_json=log_json,
            cpcm_json=cpcm_json,
            cpcm_corr_json=cpcm_corr_json,
            orcacosmo_path=orcacosmo_path,
            bundle_path=bundle_path,
        )

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
    # Build JSON for orchestrator
    # ------------------------------------------------------------
    def _build_bundle(self, lookup_id, method_used, fallback_triggered):
        inchi_key = self._get_inchi_key(lookup_id)
        smiles = self._smiles_from_inchikey(inchi_key)

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
            },
            "paths": {
                "log": os.path.join(self.parsed_outputs_dir, f"{lookup_id}.log.json"),
                "cpcm": os.path.join(self.parsed_outputs_dir, f"{lookup_id}.cpcm.json"),
                "cpcm_corr": os.path.join(self.parsed_outputs_dir, f"{lookup_id}.cpcm_corr.json"),
                "xyz": os.path.join(self.inputs_dir, f"{lookup_id}.xyz"),
            }
        }

    # ------------------------------------------------------------
    # Record success (in-memory, for big JSON + summary)
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
        inchi_key = self._get_inchi_key(lookup_id)
        smiles = self._smiles_from_inchikey(inchi_key)

        entry = {
            "lookup_id": lookup_id,
            "inchi_key": inchi_key,
            "smiles": smiles,
            "item_number": item_num,
            "orcacosmo_path": orcacosmo_path,
            "bundle_path": bundle_path,
            "raw_outputs": {
                "log": os.path.join(self.raw_outputs_dir, f"{lookup_id}.log"),
                "cpcm": os.path.join(self.raw_outputs_dir, f"{lookup_id}.cpcm"),
                "cpcm_corr": os.path.join(self.raw_outputs_dir, f"{lookup_id}.cpcm_corr"),
            },
            "method_used": method_used,
            "fallback_triggered": fallback_triggered,
            "elapsed_seconds": elapsed_seconds,
            "provenance": {
                "orca_version": self.orca_version,
                "cpcm_radii_source": self.cpcm_file,
                "engine_command_base": "orca",
                "item_number": item_num,
                "last_modified": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            },
            "parsed": {
                "log": log_json,
                "cpcm": cpcm_json,
                "cpcm_corr": cpcm_corr_json,
            },
        }

        self.successful_outputs.append(entry)
        self.log(f"[SUCCESS] Recorded output for {lookup_id}")

    # ------------------------------------------------------------
    # Summary writer (big JSON + CSV + mapping)
    # ------------------------------------------------------------
    def _write_summary(self, stage_start, stage_end):
        results_path = os.path.join(self.outputs_dir, "orcacosmo_results.json")
        mapping_path = os.path.join(self.outputs_dir, "item_to_lookup_mapping.json")
        summary_csv = os.path.join(self.outputs_dir, "orcacosmo_summary.csv")

        # Big JSON with all results
        with AtomicWriter(results_path) as f:
            json.dump(self.successful_outputs, f, indent=2)

        # Item mapping
        with AtomicWriter(mapping_path) as f:
            json.dump(self.item_to_lookup, f, indent=2)

        # Human-readable summary CSV
        rows = []
        for entry in self.successful_outputs:
            rows.append(
                {
                    "lookup_id": entry["lookup_id"],
                    "inchi_key": entry["inchi_key"],
                    "method_used": entry["method_used"],
                    "fallback_triggered": entry["fallback_triggered"],
                    "elapsed_seconds": entry["elapsed_seconds"],
                    "orcacosmo_path": entry["orcacosmo_path"],
                    "item_number": entry["item_number"],
                }
            )

        if rows:
            df = pd.DataFrame(rows)
            df.to_csv(summary_csv, index=False)

        self.log(f"Wrote ORCA COSMO results JSON with {len(self.successful_outputs)} entries: {results_path}")
        self.log(f"Wrote item-to-lookup mapping with {len(self.item_to_lookup)} entries: {mapping_path}")
        if rows:
            self.log(f"Wrote ORCA COSMO summary CSV: {summary_csv}")

        self.log(
            f"Stage elapsed time: {stage_end - stage_start:.2f} s "
            f"(successful={self.stats['successful']}, failed={self.stats['failed']})"
        )

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
    # Metadata helpers
    # ------------------------------------------------------------
    def _smiles_from_inchikey(self, inchikey):
        meta_path = os.path.join(
            self.metadata_dir,
            f"{inchikey}.json",
        )

        if not os.path.exists(meta_path):
            self.fail(f"Molecule metadata JSON not found for InChIKey: {inchikey} ({meta_path})")

        with open(meta_path) as f:
            meta = json.load(f)

        smiles = meta.get("smiles")
        if not smiles:
            self.fail(f"No SMILES field in metadata JSON for InChIKey: {inchikey}")

        return smiles

    def _get_inchi_key(self, lookup_id):
        for entry in self.energies_entries:
            if entry["lookup_id"] == lookup_id:
                return entry["inchi_key"]
        self.fail(f"InChIKey not found in energies summary for {lookup_id}")

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
