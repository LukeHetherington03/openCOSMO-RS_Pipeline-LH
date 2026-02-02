import os
import json
import shutil
import time
import traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter
from modules.solubility_engine.cosmors_wrapper import COSMORSWrapper


class SolubilityStage(BaseStage):
    """
    Modernised SolubilityStage with:
      - Parallel execution (config + args)
      - Renumbering preserved for COSMO‑RS
      - Canonical outputs:
            solubility_results.json
            solubility_summary.csv
            solute/<inchi_key>/mapping.json
      - Provenance blocks
      - Metadata validation
      - Wrapper output validation
      - Strict mode
    """

    # ------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------
    def execute(self):
        self.log_header("Starting Solubility Stage")

        self._prepare_directories()
        self._load_stage_config()
        self._load_orcacosmo_summary()
        self._load_wrapper()

        # List of solutes (one per inchi_key)
        solute_ids = [entry["inchi_key"] for entry in self.orcacosmo_entries]
        self.set_items(solute_ids)

        # Results accumulated here
        self.results = []
        self.stats = {
            "successful": 0,
            "failed": 0,
            "missing_metadata": 0,
            "missing_conformers": 0,
        }

        stage_start = time.perf_counter()

        # ------------------------------------------------------------
        # Parallel or sequential execution
        # ------------------------------------------------------------
        if self.parallel:
            self._run_parallel(solute_ids)
        else:
            self._run_sequential(solute_ids)

        stage_end = time.perf_counter()

        self._write_global_outputs(stage_start, stage_end)
        self._log_warning_summary()

        self.log_header("Solubility Stage Complete")
        self.job.mark_complete()

    # ------------------------------------------------------------
    # Directory setup
    # ------------------------------------------------------------
    def _prepare_directories(self):
        self.cosmo_root = os.path.join(self.outputs_dir, "cosmo")
        self.solute_cosmo_root = os.path.join(self.cosmo_root, "solute")
        self.solvent_cosmo_root = os.path.join(self.cosmo_root, "solvent")

        self.workdir = os.path.join(self.outputs_dir, "work")
        self.results_dir = os.path.join(self.outputs_dir, "results")

        logs_dir = os.path.join(self.job.job_dir, "logs")
        os.makedirs(logs_dir, exist_ok=True)
        self.debug_log_path = os.path.join(logs_dir, "solubility_debug.log")

        for d in [
            self.cosmo_root,
            self.solute_cosmo_root,
            self.solvent_cosmo_root,
            self.workdir,
            self.results_dir,
        ]:
            os.makedirs(d, exist_ok=True)

    # ------------------------------------------------------------
    # Load config
    # ------------------------------------------------------------
    def _load_stage_config(self):
        self.summary_file = self.parameters.get("summary_file")
        if not self.summary_file:
            self.fail("SolubilityStage requires 'summary_file'")
        if not os.path.exists(self.summary_file):
            self.fail(f"summary_file does not exist: {self.summary_file}")

        cfg = self.config
        if not cfg:
            self.fail("SolubilityStage: missing config in job.config")

        const_cfg = cfg.get("constant_files", {})
        self.chemistry_dir = const_cfg.get("chemistry_dir")
        self.metadata_dir = const_cfg.get("metadata_dir")
        self.solvent_dir = const_cfg.get("solvent_dir")

        sol_cfg = cfg.get("solubility", {})
        defaults_file = sol_cfg.get("defaults")

        with open(defaults_file) as f:
            defaults = json.load(f)

        # Defaults from JSON
        self.default_solvent = defaults.get("default_solvent", "water")
        self.default_temperature = defaults.get("default_temperature", 298.15)
        self.default_SORcf = defaults.get("default_SORcf", 1.0)
        self.default_parallel = defaults.get("parallel", True)
        self.default_workers = defaults.get("n_workers", 4)

        # Overrideable via args
        self.solvent_name = self.parameters.get("solvent_name", self.default_solvent)
        self.temperature = self.parameters.get("temperature", self.default_temperature)
        self.SORcf = self.parameters.get("SORcf", self.default_SORcf)
        self.calculations = self.parameters.get("calculations", "mixed_only")

        self.parallel = self.parameters.get("parallel", self.default_parallel)
        self.n_workers = self.parameters.get("n_workers", self.default_workers)

        self.strict_mode = bool(sol_cfg.get("strict", False))

    # ------------------------------------------------------------
    # Load COSMO summary
    # ------------------------------------------------------------
    def _load_orcacosmo_summary(self):
        with open(self.summary_file) as f:
            raw_entries = json.load(f)

        grouped = {}
        for entry in raw_entries:
            inchi_key = entry.get("inchi_key")
            path = entry.get("orcacosmo_path")
            lookup_id = entry.get("lookup_id")
            if not inchi_key or not path:
                continue
            grouped.setdefault(inchi_key, []).append((lookup_id, path))

        self.orcacosmo_entries = [
            {
                "inchi_key": inchi_key,
                "conformers": paths,  # list of (lookup_id, path)
            }
            for inchi_key, paths in grouped.items()
        ]

    # ------------------------------------------------------------
    # Load wrapper
    # ------------------------------------------------------------
    def _load_wrapper(self):
        opencosmo_cfg = self.config.get("opencosmo", {})
        self.wrapper = COSMORSWrapper(
            config_path=self.config["solubility"]["defaults"],
            opencosmo_paths=opencosmo_cfg,
        )

    # ------------------------------------------------------------
    # Metadata loader
    # ------------------------------------------------------------
    def _load_metadata(self, inchi_key):
        global_path = os.path.join(self.metadata_dir, f"{inchi_key}.json")
        job_path = os.path.join(self.inputs_dir, "molecule_metadata", f"{inchi_key}.json")

        if os.path.exists(global_path):
            with open(global_path) as f:
                return json.load(f)

        if os.path.exists(job_path):
            with open(job_path) as f:
                return json.load(f)

        self.stats["missing_metadata"] += 1
        return {
            "lookup_id": inchi_key,
            "melting_temp": "N/A",
            "Hfus": "N/A",
            "Gfus_model": "MyrdalYalkowsky",
            "smiles": "",
        }

    # ------------------------------------------------------------
    # Parallel execution
    # ------------------------------------------------------------
    def _run_parallel(self, solute_ids):
        self.log(f"[INFO] Running in parallel with {self.n_workers} workers")

        with ProcessPoolExecutor(max_workers=self.n_workers) as exe:
            futures = {
                exe.submit(self._process_solute_safe, inchi_key): inchi_key
                for inchi_key in solute_ids
            }

            for fut in as_completed(futures):
                inchi_key = futures[fut]
                try:
                    result = fut.result()
                    if result:
                        self.results.append(result)
                        self.stats["successful"] += 1
                except Exception as e:
                    self.stats["failed"] += 1
                    self.log(f"[ERROR] {inchi_key}: {e}")

    # ------------------------------------------------------------
    # Sequential execution
    # ------------------------------------------------------------
    def _run_sequential(self, solute_ids):
        for inchi_key in solute_ids:
            try:
                result = self._process_solute_safe(inchi_key)
                if result:
                    self.results.append(result)
                    self.stats["successful"] += 1
            except Exception as e:
                self.stats["failed"] += 1
                self.log(f"[ERROR] {inchi_key}: {e}")

    # ------------------------------------------------------------
    # Safe wrapper for parallel execution
    # ------------------------------------------------------------
    def _process_solute_safe(self, inchi_key):
        try:
            return self._process_solute(inchi_key)
        except Exception as e:
            if self.strict_mode:
                raise
            tb = traceback.format_exc()
            self.log(f"[ERROR] {inchi_key}: {e}")
            self.log(f"[DEBUG] Traceback:\n{tb}")
            return None

    # ------------------------------------------------------------
    # Process a single solute
    # ------------------------------------------------------------
    def _process_solute(self, inchi_key):
        entry = next((e for e in self.orcacosmo_entries if e["inchi_key"] == inchi_key), None)
        if not entry:
            raise ValueError(f"No COSMO entry found for {inchi_key}")

        conformers = entry["conformers"]  # list of (lookup_id, path)
        metadata = self._load_metadata(inchi_key)

        if not conformers:
            self.stats["missing_conformers"] += 1
            raise RuntimeError(f"No conformers found for {inchi_key}")

        # ------------------------------------------------------------
        # Renumber solute conformers (COSMO‑RS requirement)
        # ------------------------------------------------------------
        solute_dir, mapping = self._copy_and_renumber_conformers(inchi_key, conformers)

        # ------------------------------------------------------------
        # Renumber solvent conformers
        # ------------------------------------------------------------
        solvent_dir = self._copy_and_renumber_solvent()

        # ------------------------------------------------------------
        # Run COSMO‑RS wrapper
        # ------------------------------------------------------------
        start = time.perf_counter()
        wrapper_out = self.wrapper.run(
            job_dir=Path(self.workdir),
            solute_name=inchi_key,
            solute_smiles=metadata.get("smiles", ""),
            solute_Tm=metadata.get("melting_temp", "N/A"),
            solute_dir=Path(solute_dir),
            solvent_name=self.solvent_name,
            solvent_dir=Path(solvent_dir),
        )
        elapsed = time.perf_counter() - start

        # Validate wrapper output
        if "result" not in wrapper_out or "x_solubility" not in wrapper_out["result"]:
            raise RuntimeError(f"COSMO‑RS returned incomplete result for {inchi_key}")

        sol = wrapper_out["result"]["x_solubility"]

        # ------------------------------------------------------------
        # Write raw COSMO‑RS output
        # ------------------------------------------------------------
        raw_path = os.path.join(self.results_dir, f"{inchi_key}_raw.txt")
        with open(raw_path, "w") as f:
            f.write(wrapper_out["raw_output"])

        # ------------------------------------------------------------
        # Write per‑solute JSON
        # ------------------------------------------------------------
        out_json = os.path.join(self.results_dir, f"{inchi_key}.json")
        with AtomicWriter(out_json) as f:
            json.dump(
                {
                    "inchi_key": inchi_key,
                    "smiles": metadata.get("smiles", ""),
                    "solvent": self.solvent_name,
                    "temperature_K": self.temperature,
                    "solubility_x": sol,
                    "n_solute_confs": wrapper_out["n_solute"],
                    "n_solvent_confs": wrapper_out["n_solvent"],
                    "renumbered_mapping": mapping,
                    "raw_output_path": raw_path,
                    "provenance": {
                        "cosmors_version": wrapper_out.get("version", "unknown"),
                        "parallel": self.parallel,
                        "n_workers": self.n_workers,
                        "elapsed_seconds": elapsed,
                        "last_modified": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    },
                },
                f,
                indent=2,
            )

        # ------------------------------------------------------------
        # Write human summary
        # ------------------------------------------------------------
        human_path = os.path.join(self.results_dir, f"{inchi_key}_summary.txt")
        with open(human_path, "w") as f:
            f.write(
                f"Solubility Summary\n"
                f"------------------\n"
                f"Molecule: {inchi_key}\n"
                f"SMILES: {metadata.get('smiles','')}\n"
                f"Solvent: {self.solvent_name}\n"
                f"Temperature: {self.temperature} K\n"
                f"Mole fraction solubility: {sol}\n"
            )

        return {
            "inchi_key": inchi_key,
            "solubility_x": sol,
            "json": out_json,
            "summary": human_path,
            "raw_output": raw_path,
            "n_confs": wrapper_out["n_solute"],
        }

    # ------------------------------------------------------------
    # Renumber solute conformers (COSMO‑RS requirement)
    # ------------------------------------------------------------
    def _copy_and_renumber_conformers(self, inchi_key, conformers):
        out_dir = os.path.join(self.solute_cosmo_root, inchi_key)
        os.makedirs(out_dir, exist_ok=True)

        mapping = {}
        for i, (lookup_id, path) in enumerate(conformers):
            if not os.path.exists(path):
                continue
            new_name = f"{inchi_key}_c{i:03d}.orcacosmo"
            shutil.copy(path, os.path.join(out_dir, new_name))
            mapping[f"c{i:03d}"] = lookup_id

        if not mapping:
            raise RuntimeError(f"No valid solute conformers found for {inchi_key}")

        # Write mapping.json
        with AtomicWriter(os.path.join(out_dir, "mapping.json")) as f:
            json.dump(mapping, f, indent=2)

        return out_dir, mapping

    # ------------------------------------------------------------
    # Renumber solvent conformers
    # ------------------------------------------------------------
    def _copy_and_renumber_solvent(self):
        solvent_path = os.path.join(self.solvent_dir, self.solvent_name)
        out_dir = os.path.join(self.solvent_cosmo_root, self.solvent_name)
        os.makedirs(out_dir, exist_ok=True)

        confs = sorted(Path(solvent_path).glob("*.orcacosmo"))
        if not confs:
            raise RuntimeError(f"No solvent conformers found in {solvent_path}")

        for i, path in enumerate(confs):
            new_name = f"{self.solvent_name}_c{i:03d}.orcacosmo"
            shutil.copy(path, os.path.join(out_dir, new_name))

        return out_dir

    # ------------------------------------------------------------
    # Write global outputs
    # ------------------------------------------------------------
    def _write_global_outputs(self, stage_start, stage_end):
        results_path = os.path.join(self.outputs_dir, "solubility_results.json")
        summary_csv = os.path.join(self.outputs_dir, "solubility_summary.csv")

        # Big JSON
        with AtomicWriter(results_path) as f:
            json.dump(self.results, f, indent=2)

        # CSV summary
        import pandas as pd
        if self.results:
            df = pd.DataFrame(self.results)
            df.to_csv(summary_csv, index=False)

        self.log(f"Wrote solubility_results.json with {len(self.results)} entries")
        self.log(f"Wrote solubility_summary.csv")

        self.log(
            f"Stage elapsed time: {stage_end - stage_start:.2f} s "
            f"(successful={self.stats['successful']}, failed={self.stats['failed']})"
        )

    # ------------------------------------------------------------
    # Warning summary
    # ------------------------------------------------------------
    def _log_warning_summary(self):
        self.log_header("Solubility Summary")
        self.log(f"Successful: {self.stats['successful']}")
        self.log(f"Failed: {self.stats['failed']}")
        self.log(f"Missing metadata: {self.stats['missing_metadata']}")
        self.log(f"Missing conformers: {self.stats['missing_conformers']}")
        self.log_header("End Solubility Summary")
