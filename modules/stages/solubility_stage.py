import os
import json
import shutil
import traceback
from pathlib import Path

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter
from modules.utils.mixture_inputs_builder import build_mixture_inputs
from modules.solubility_engine.legacy_cpp_wrapper import run_legacy_cosmors


class SolubilityStage(BaseStage):
    """
    Modern Solubility Stage using the new COSMO‑RS Python wrapper.

    Contract:
      - Stage input: orcacosmo_summary.json from OrcacosmoStage
      - For each InChI key:
          results/<inchi_key>/
              inputs/
                  solute_raw/   ← original .orcacosmo (unmodified)
                  solute/       ← renamed <InChIKey>_cNNN.orcacosmo (engine input)
                  solvent/      ← solvent .orcacosmo (as-is)
              mixture_inputs.txt
              raw_output.txt
              solubility.json
              summary.txt
    """

    # ------------------------------------------------------------
    # Stage entrypoint
    # ------------------------------------------------------------
    def execute(self):
        self.log_header("Starting Solubility Stage")

        stage_input = self.get_stage_input()
        if not stage_input:
            self.fail("SolubilityStage requires stage_input (orcacosmo_summary.json)")

        self._prepare_directories()
        self._load_stage_config()
        self._load_orcacosmo_summary(stage_input)

        self.set_stage_output("solubility_results.json")

        solute_ids = [entry["inchi_key"] for entry in self.orcacosmo_entries]
        self.set_items(solute_ids)

        self.successful = []
        self.failed = []

        for inchi_key in list(self.job.pending_items):
            try:
                self._process_solute(inchi_key)
                self.update_progress(inchi_key)
            except Exception as e:
                tb = traceback.format_exc()
                self.log(f"[ERROR] {inchi_key}: {e}")
                self._debug(tb)
                self.failed.append(inchi_key)
                self.update_progress(inchi_key, success=False)

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
        self.results_dir = Path(self.outputs_dir) / "results"

        logs_dir = Path(self.job.job_dir) / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        self.debug_log_path = logs_dir / "solubility_debug.log"

        self.results_dir.mkdir(parents=True, exist_ok=True)

        self._debug(f"Prepared directories: results={self.results_dir}")

    # ------------------------------------------------------------
    # Load config
    # ------------------------------------------------------------
    def _load_stage_config(self):
        cfg = self.config["solubility"]
        defaults_path = cfg["defaults"]

        with open(defaults_path) as f:
            self.defaults = json.load(f)

        const_cfg = self.config.get("constant_files", {})
        self.metadata_dir = const_cfg.get("metadata_dir")
        self.solvent_dir = const_cfg.get("solvent_dir")

        opencosmo_cfg = self.config.get("opencosmo", {})
        self.opencosmo_python_src = opencosmo_cfg["python_src"]
        self.opencosmo_cpp_bindings = opencosmo_cfg["cpp_bindings"]
        self.opencosmo_driver = opencosmo_cfg["python_driver"]

        self.temperature = self.parameters.get(
            "temperature", self.defaults.get("default_temperature", 298.15)
        )
        self.solvent_name = self.parameters.get(
            "solvent_name", self.defaults.get("default_solvent", "water")
        )
        self.verbose = self.parameters.get("verbose", False)

        self._debug(
            f"Config loaded: solvent={self.solvent_name}, T={self.temperature}, "
            f"python_src={self.opencosmo_python_src}, cpp_bindings={self.opencosmo_cpp_bindings}"
        )

    # ------------------------------------------------------------
    # Load ORCA‑COSMO summary
    # ------------------------------------------------------------
    def _load_orcacosmo_summary(self, stage_input):
        stage_input = self.require_file(stage_input, "orcacosmo_summary.json")

        with open(stage_input) as f:
            raw_entries = json.load(f)

        grouped = {}
        for entry in raw_entries:
            inchi_key = entry["inchi_key"]
            path = entry["orcacosmo_path"]
            lookup_id = entry["lookup_id"]

            grouped.setdefault(inchi_key, {"lookup_id": lookup_id, "paths": []})
            grouped[inchi_key]["paths"].append(path)

        self.orcacosmo_entries = [
            {
                "lookup_id": data["lookup_id"],
                "inchi_key": inchi_key,
                "paths": data["paths"],
            }
            for inchi_key, data in grouped.items()
        ]

        self.log(f"Loaded {len(self.orcacosmo_entries)} molecules")

    # ------------------------------------------------------------
    # Process a single solute (per molecule)
    # ------------------------------------------------------------
    def _process_solute(self, inchi_key):
        entry = next(e for e in self.orcacosmo_entries if e["inchi_key"] == inchi_key)
        paths = entry["paths"]
        lookup_id = entry["lookup_id"]

        self.log(f"[INFO] === Solubility calculation started for {inchi_key} ===")

        # Load metadata
        meta = self._load_metadata(inchi_key)
        if self.verbose:
            self.log(f"[VERBOSE] Metadata loaded: {meta}")

        # Per-molecule directory structure
        mol_dir = self.results_dir / inchi_key
        inputs_dir = mol_dir / "inputs"
        solute_raw_dir = inputs_dir / "solute_raw"
        solute_dir = inputs_dir / "solute"
        solvent_dir = inputs_dir / "solvent"

        mol_dir.mkdir(parents=True, exist_ok=True)
        inputs_dir.mkdir(parents=True, exist_ok=True)
        solute_raw_dir.mkdir(parents=True, exist_ok=True)
        solute_dir.mkdir(parents=True, exist_ok=True)
        solvent_dir.mkdir(parents=True, exist_ok=True)

        # Copy solute conformers AS-IS into solute_raw (provenance)
        for p in paths:
            shutil.copy(p, solute_raw_dir / Path(p).name)

        # Copy + rename solute conformers into solute/ for COSMO‑RS engine
        for i, p in enumerate(sorted(paths)):
            new_name = f"{inchi_key}_c{i:03d}.orcacosmo"
            shutil.copy(p, solute_dir / new_name)

        # Copy solvent conformers AS-IS
        solvent_src = Path(self.solvent_dir) / self.solvent_name
        for p in sorted(solvent_src.glob("*.orcacosmo")):
            shutil.copy(p, solvent_dir / p.name)

        # Count conformers
        n_solute = len(list(solute_dir.glob("*.orcacosmo")))
        n_solvent = len(list(solvent_dir.glob("*.orcacosmo")))

        # Build mixture_inputs
        mixture_text = build_mixture_inputs(
            solute_meta={
                "inchi_key": inchi_key,
                "smiles": meta["smiles"],
                "melting_temp": meta["Tm"],
                "Gfus_mode": meta["Gfus_mode"],
                "Hfus": meta["Hfus"],
            },
            solute_dir=solute_dir,
            solvent_dir=solvent_dir,
            n_solute_confs=n_solute,
            n_solvent_confs=n_solvent,
            defaults=self.defaults,
            temperature=self.temperature,
        )

        mix_path = mol_dir / "mixture_inputs.txt"
        mix_path.write_text(mixture_text)

        if self.verbose:
            self.log(f"[VERBOSE] mixture_inputs.txt written to {mix_path}")

        # Run COSMO‑RS
        result = run_legacy_cosmors(
            mixture_text,
            python_src=self.opencosmo_python_src,
            cpp_bindings=self.opencosmo_cpp_bindings,
            driver_script=self.opencosmo_driver,
        )

        # Write raw output
        raw_path = mol_dir / "raw_output.txt"
        raw_path.write_text(result["raw_stdout"] or "")

        # Optionally preserve import validation file per molecule
        validation_src = Path(result.get("import_validation_file", ""))
        if validation_src.exists():
            shutil.copy(validation_src, mol_dir / "import_validation.txt")

        # Write JSON result
        out_json = mol_dir / "solubility.json"
        with AtomicWriter(out_json) as f:
            json.dump(
                {
                    "lookup_id": lookup_id,
                    "inchi_key": inchi_key,
                    "smiles": meta["smiles"],
                    "melting_temp": meta["Tm"],
                    "experimental_solubility": meta["experimental_solubility_mol_frac"],
                    "predicted_solubility": result["solubility"],
                    "temperature_K": self.temperature,
                    "n_solute_confs": n_solute,
                    "n_solvent_confs": n_solvent,
                    "raw_output_path": str(raw_path),
                },
                f,
                indent=2,
            )

        # Write human summary
        summary_path = mol_dir / "summary.txt"
        summary = f"""
Solubility Summary
------------------
Molecule: {inchi_key}
SMILES: {meta['smiles']}
Melting Temp: {meta['Tm']}
Experimental solubility: {meta['experimental_solubility_mol_frac']}
Predicted solubility: {result['solubility']}
Temperature: {self.temperature} K
"""
        summary_path.write_text(summary.strip() + "\n")

        self.successful.append(
            {
                "lookup_id": lookup_id,
                "inchi_key": inchi_key,
                "predicted": result["solubility"],
                "experimental": meta["experimental_solubility_mol_frac"],
                "json": str(out_json),
                "summary": str(summary_path),
            }
        )

        self.log(f"[INFO] === Solubility calculation finished for {inchi_key} ===")

    # ------------------------------------------------------------
    # Metadata loader
    # ------------------------------------------------------------
    def _load_metadata(self, inchi_key: str):
        meta_path = Path(self.metadata_dir) / f"{inchi_key}.json"

        with open(meta_path, "r") as f:
            meta = json.load(f)

        Tm = meta.get("melting_temp", "N/A")
        try:
            Tm = float(Tm)
        except Exception:
            pass

        return {
            "inchi_key": meta["inchi_key"],
            "smiles": meta.get("smiles", ""),
            "Tm": Tm,
            "Gfus_mode": meta.get("Gfus_mode", "MyrdalYalkowsky"),
            "Hfus": meta.get("Hfus", "N/A"),
            "experimental_solubility_mol_frac": meta.get(
                "experimental_solubility_mol_frac", None
            ),
            "mol_name": meta.get("mol_name", meta["inchi_key"]),
        }

    # ------------------------------------------------------------
    # Write global summary
    # ------------------------------------------------------------
    def _write_summary(self):
        results_path = self.get_stage_output()
        with AtomicWriter(results_path) as f:
            json.dump(self.successful, f, indent=2)

        self.log(f"Wrote {results_path} with {len(self.successful)} entries")

        human_path = Path(self.outputs_dir) / "solubility_human_summary.txt"

        lines = []
        lines.append("Solubility Stage Summary")
        lines.append("========================")
        lines.append("")
        lines.append(f"Total molecules: {len(self.successful) + len(self.failed)}")
        lines.append(f"Successful:      {len(self.successful)}")
        lines.append(f"Failed:          {len(self.failed)}")
        lines.append("")

        for entry in self.successful:
            inchi = entry["inchi_key"]
            lookup = entry["lookup_id"]
            pred = entry["predicted"]
            exp = entry["experimental"]

            abs_err = abs(pred - exp) if pred is not None and exp is not None else None

            lines.append("------------------------------------------------------------")
            lines.append(f"{inchi} (lookup: {lookup})")
            lines.append(f"Predicted:    {pred}")
            lines.append(f"Experimental: {exp}")
            lines.append(f"Abs Error:    {abs_err}")
            lines.append(f"JSON:         {entry['json']}")
            lines.append(f"Summary:      {entry['summary']}")
            lines.append("")

        if self.failed:
            lines.append("------------------------------------------------------------")
            lines.append("Failed Molecules")
            lines.append("------------------------------------------------------------")
            for f in self.failed:
                lines.append(f"- {f}")
            lines.append("")

        human_path.write_text("\n".join(lines))
        self.log(f"Wrote human-readable summary: {human_path}")
