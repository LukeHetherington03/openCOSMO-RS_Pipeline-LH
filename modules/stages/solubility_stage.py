import os
import sys
import json
import shutil
from pathlib import Path
import traceback

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter
from modules.solubility_engine.cosmors_wrapper import COSMORSWrapper


class SolubilityStage(BaseStage):

    def execute(self):
        self.log_header("Starting Solubility Stage")

        self._prepare_directories()
        self._load_stage_config()
        self._load_orcacosmo_summary()
        self._load_wrapper()

        solute_ids = [entry["lookup_id"] for entry in self.orcacosmo_entries]
        self.set_items(solute_ids)

        self.successful_calcs = []

        for lookup_id in list(self.job.pending_items):
            try:
                self._process_solute(lookup_id)
                self.update_progress(lookup_id)
            except Exception as e:
                tb = traceback.format_exc()
                self.log(f"[ERROR] {lookup_id}: {e}")
                self._debug(f"Full traceback:\n{tb}")
                self.update_progress(lookup_id, success=False)

        self._write_summary()
        self.log_header("Solubility Stage Complete")
        self.job.mark_complete()

    # ------------------------------------------------------------
    # Debug helper
    # ------------------------------------------------------------
    def _debug(self, msg: str):
        self.log(f"[DEBUG] {msg}")
        try:
            with open(self.debug_log_path, "a") as f:
                f.write(msg + "\n")
        except Exception:
            pass

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

        os.makedirs(self.cosmo_root, exist_ok=True)
        os.makedirs(self.solute_cosmo_root, exist_ok=True)
        os.makedirs(self.solvent_cosmo_root, exist_ok=True)
        os.makedirs(self.workdir, exist_ok=True)
        os.makedirs(self.results_dir, exist_ok=True)

        self._debug(
            f"Prepared directories: cosmo_root={self.cosmo_root}, "
            f"workdir={self.workdir}, results_dir={self.results_dir}"
        )

    # ------------------------------------------------------------
    # Load config
    # ------------------------------------------------------------
    def _load_stage_config(self):
        self.summary_file = self.parameters.get("summary_file")
        if not self.summary_file:
            self.fail("SolubilityStage requires 'summary_file' in parameters")
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

        self.default_solvent = defaults.get("default_solvent", "water")
        self.default_temperature = defaults.get("default_temperature", 298.15)
        self.default_SORcf = defaults.get("default_SORcf", 1.0)

        self.solvent_name = self.parameters.get("solvent_name", self.default_solvent)
        self.temperature = self.parameters.get("temperature", self.default_temperature)
        self.SORcf = self.parameters.get("SORcf", self.default_SORcf)
        self.calculations = self.parameters.get("calculations", "mixed_only")

        self.GLOBAL_METADATA_DIR = self.metadata_dir
        self.GLOBAL_SOLVENT_DIR = self.solvent_dir

        self._debug(
            f"Config loaded: solvent={self.solvent_name}, T={self.temperature}, "
            f"SORcf={self.SORcf}, calculations={self.calculations}"
        )

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
            if not inchi_key or not path:
                continue
            grouped.setdefault(inchi_key, []).append(path)

        self.orcacosmo_entries = [
            {
                "lookup_id": inchi_key,
                "inchi_key": inchi_key,
                "conformer_paths": paths,
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
    # Load metadata
    # ------------------------------------------------------------
    def _load_metadata(self, inchi_key):
        global_path = os.path.join(self.GLOBAL_METADATA_DIR, f"{inchi_key}.json")
        job_path = os.path.join(self.inputs_dir, "molecule_metadata", f"{inchi_key}.json")

        if os.path.exists(global_path):
            with open(global_path) as f:
                return json.load(f)

        if os.path.exists(job_path):
            with open(job_path) as f:
                return json.load(f)

        return {
            "lookup_id": inchi_key,
            "melting_temp": "N/A",
            "Hfus": "N/A",
            "Gfus_model": "MyrdalYalkowsky",
            "smiles": "",
        }

    # ------------------------------------------------------------
    # Process a single solute
    # ------------------------------------------------------------
    def _process_solute(self, lookup_id):
        entry = next((e for e in self.orcacosmo_entries if e["lookup_id"] == lookup_id), None)
        if not entry:
            raise ValueError(f"No COSMO entry found for {lookup_id}")

        inchi_key = entry["inchi_key"]
        conformer_paths = entry["conformer_paths"]
        metadata = self._load_metadata(inchi_key)

        solute_dir = self._copy_and_renumber_conformers(inchi_key, conformer_paths)
        solvent_dir = self._copy_and_renumber_solvent()

        # Run wrapper (now returns raw_output + result)
        wrapper_out = self.wrapper.run(
            job_dir=Path(self.workdir),
            solute_name=inchi_key,
            solute_smiles=metadata.get("smiles", ""),
            solute_Tm=metadata.get("melting_temp", "N/A"),
            solute_dir=Path(solute_dir),
            solvent_name=self.solvent_name,
            solvent_dir=Path(solvent_dir),
        )

        # ------------------------------------------------------------
        # Write raw COSMO‑RS output
        # ------------------------------------------------------------
        raw_path = os.path.join(self.workdir, f"{inchi_key}_raw_output.txt")
        with open(raw_path, "w") as f:
            f.write(wrapper_out["raw_output"])

        # ------------------------------------------------------------
        # Write machine‑readable JSON
        # ------------------------------------------------------------
        out_json = os.path.join(self.results_dir, f"{inchi_key}_solubility.json")
        with AtomicWriter(out_json) as f:
            json.dump(
                {
                    "request_id": self.request.request_id,
                    "job_id": self.job.job_id,
                    "lookup_id": inchi_key,
                    "inchi_key": inchi_key,
                    "smiles": metadata.get("smiles", ""),
                    "solvent": self.solvent_name,
                    "temperature_K": self.temperature,
                    "n_solute_confs": wrapper_out["n_solute"],
                    "n_solvent_confs": wrapper_out["n_solvent"],
                    "solubility_result": wrapper_out["result"],
                    "raw_output_path": raw_path,
                },
                f,
                indent=2,
            )

        # ------------------------------------------------------------
        # Write human‑readable summary
        # ------------------------------------------------------------
        human_path = os.path.join(self.results_dir, f"{inchi_key}_summary.txt")
        sol = wrapper_out["result"].get("x_solubility", None)

        summary_txt = f"""
Solubility Summary
------------------
Request: {self.request.request_id}
Job: {self.job.job_id}

Molecule: {inchi_key}
SMILES: {metadata.get("smiles","")}
Solvent: {self.solvent_name}
Temperature: {self.temperature} K

Mole fraction solubility: {sol}
"""

        with open(human_path, "w") as f:
            f.write(summary_txt.strip() + "\n")

        # Track success
        self.successful_calcs.append({
            "inchi_key": inchi_key,
            "solubility": sol,
            "json": out_json,
            "summary": human_path,
            "raw_output": raw_path,
        })

    # ------------------------------------------------------------
    # Copy solute conformers
    # ------------------------------------------------------------
    def _copy_and_renumber_conformers(self, inchi_key, conformer_paths):
        out_dir = os.path.join(self.solute_cosmo_root, inchi_key)
        os.makedirs(out_dir, exist_ok=True)

        count = 0
        for path in conformer_paths:
            if not os.path.exists(path):
                continue
            new_name = f"{inchi_key}_c{count:03d}.orcacosmo"
            shutil.copy(path, os.path.join(out_dir, new_name))
            count += 1

        if count == 0:
            raise RuntimeError(f"No valid solute conformers found for {inchi_key}")

        return out_dir

    # ------------------------------------------------------------
    # Copy solvent conformers
    # ------------------------------------------------------------
    def _copy_and_renumber_solvent(self):
        solvent_path = os.path.join(self.GLOBAL_SOLVENT_DIR, self.solvent_name)
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
    # Write global summary
    # ------------------------------------------------------------
    def _write_summary(self):
        summary_path = os.path.join(self.outputs_dir, "solubility_summary.json")
        with AtomicWriter(summary_path) as f:
            json.dump(
                {
                    "request_id": self.job.get_request_id(),
                    "job_id": self.get_job_id(),
                    "n_molecules": len(self.successful_calcs),
                    "results": self.successful_calcs,
                },
                f,
                indent=2,
            )
        self.log(f"Wrote solubility summary: {summary_path}")
