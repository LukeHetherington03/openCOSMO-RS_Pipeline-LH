import os
import json
import pandas as pd
import numpy as np

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter


class PruningStage(BaseStage):
    """
    Prunes conformers for each molecule group.

    Features:
      - Molecule-level job items (clean logs)
      - Conformer-level pruning logic
      - Configurable pruning strategies:
            * keep_lowest_n
            * energy_window
            * all (no pruning)
            * rmsd (placeholder)
      - Removes conformers with missing energies (logged)
      - Produces:
            * pruned energies.json
            * pruning_summary.csv (human readable)
    """

    # ------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------
    def execute(self):
        self.log_header("Starting Pruning Stage")

        params = self.parameters

        # ------------------------------------------------------------
        # 🔧 Patch: Map pipeline-style args → pruning-stage args
        # ------------------------------------------------------------
        if "strategy" in params:
            params["pruning_strategy"] = params["strategy"]

        if "strategy_params" in params and isinstance(params["strategy_params"], dict):
            for k, v in params["strategy_params"].items():
                params[k] = v

        # Alias: top_n → keep_lowest_n
        if params.get("pruning_strategy") == "top_n":
            params["pruning_strategy"] = "keep_lowest_n"

        self.log(f"[DEBUG] Using pruning strategy: {params.get('pruning_strategy')}")
        self.log(f"[DEBUG] Strategy parameters: { {k: params[k] for k in params if k not in ('summary_file')} }")

        summary_file = params.get("summary_file")
        if summary_file is None:
            self.fail("PruningStage requires summary_file from previous stage.")

        if not os.path.isfile(summary_file):
            self.fail(f"summary_file does not exist: {summary_file}")

        # Load entries
        entries = self._load_entries(summary_file)

        # Group by molecule
        groups = self._group_by_molecule(entries)

        # Job items = molecule IDs
        molecule_ids = list(groups.keys())
        self.set_items(molecule_ids)

        # Prune each molecule group
        pruned_entries = []
        summary_rows = []

        for mol_id, conformers in groups.items():
            try:
                survivors, summary_row = self._prune_group(mol_id, conformers, params)
                pruned_entries.extend(survivors)
                summary_rows.append(summary_row)

                self.update_progress(mol_id, success=True)
                self.log(f"[INFO] Pruning complete for {mol_id}: kept {len(survivors)} conformers")

            except Exception as e:
                self.log(f"[ERROR] Pruning failed for {mol_id}: {e}")
                self.update_progress(mol_id, success=False)

        # Write outputs
        self._write_outputs(pruned_entries, summary_rows)

        self.log_header("Pruning Stage Complete")
        self.job.mark_complete()

    # ------------------------------------------------------------
    # Load entries
    # ------------------------------------------------------------
    def _load_entries(self, summary_file):
        try:
            with open(summary_file) as f:
                entries = json.load(f)
        except Exception as e:
            self.fail(f"Failed to load summary_file: {e}")

        if not isinstance(entries, list):
            self.fail("summary_file must contain a list of conformer entries.")

        return entries

    # ------------------------------------------------------------
    # Group conformers by molecule
    # ------------------------------------------------------------
    def _group_by_molecule(self, entries):
        groups = {}
        for e in entries:
            lookup = e.get("lookup_id")
            if not lookup:
                self.log("[WARNING] Entry missing lookup_id, skipping")
                continue

            # Extract molecule ID (everything before _confXXX)
            if "_conf" in lookup:
                mol_id = lookup.split("_conf")[0]
            else:
                mol_id = lookup  # fallback

            groups.setdefault(mol_id, []).append(e)

        self.log(f"Found {len(groups)} molecule groups")
        return groups

    # ------------------------------------------------------------
    # Pruning logic for a single molecule
    # ------------------------------------------------------------
    def _prune_group(self, mol_id, conformers, params):
        """
        Returns:
            survivors: list of conformer dicts
            summary_row: dict for pruning_summary.csv
        """

        # Remove conformers with missing energies
        valid = []
        removed_missing_energy = 0

        for e in conformers:
            if e.get("energy") is None or (isinstance(e["energy"], float) and np.isnan(e["energy"])):
                removed_missing_energy += 1
                self.log(f"[WARNING] {mol_id}: conformer {e.get('lookup_id')} removed (missing energy)")
            else:
                valid.append(e)

        if not valid:
            raise RuntimeError(f"All conformers for {mol_id} missing energies.")

        # Choose pruning strategy
        strategy = params.get("pruning_strategy", "keep_lowest_n").lower()

        if strategy == "keep_lowest_n":
            survivors = self._prune_keep_lowest_n(valid, params)

        elif strategy == "energy_window":
            survivors = self._prune_energy_window(valid, params)

        elif strategy == "all":
            survivors = valid

        elif strategy == "rmsd":
            survivors = self._prune_rmsd(valid, params)  # placeholder

        else:
            raise RuntimeError(f"Unknown pruning_strategy: {strategy}")

        # Build summary row
        summary_row = {
            "molecule_id": mol_id,
            "total_conformers": len(conformers),
            "valid_energy_conformers": len(valid),
            "removed_missing_energy": removed_missing_energy,
            "kept_after_pruning": len(survivors),
            "strategy": strategy,
        }

        return survivors, summary_row

    # ------------------------------------------------------------
    # Strategy: keep lowest N energies
    # ------------------------------------------------------------
    def _prune_keep_lowest_n(self, conformers, params):
        n = params.get("n", 5)
        sorted_conf = sorted(conformers, key=lambda e: e["energy"])
        return sorted_conf[:n]

    # ------------------------------------------------------------
    # Strategy: keep all within ΔE of minimum
    # ------------------------------------------------------------
    def _prune_energy_window(self, conformers, params):
        window = params.get("pruning_energy_window", 5.0)  # kcal/mol
        sorted_conf = sorted(conformers, key=lambda e: e["energy"])
        min_e = sorted_conf[0]["energy"]
        return [e for e in sorted_conf if (e["energy"] - min_e) <= window]

    # ------------------------------------------------------------
    # Strategy: RMSD clustering (placeholder)
    # ------------------------------------------------------------
    def _prune_rmsd(self, conformers, params):
        self.log("[WARNING] RMSD pruning not implemented; keeping all conformers")
        return conformers

    # ------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------
    def _write_outputs(self, pruned_entries, summary_rows):
        outputs_dir = self.outputs_dir

        # energies.json
        energies_path = os.path.join(outputs_dir, "energies.json")
        with AtomicWriter(energies_path) as f:
            json.dump(pruned_entries, f, indent=2)

        # pruning_summary.csv
        summary_path = os.path.join(outputs_dir, "pruning_summary.csv")
        df = pd.DataFrame(summary_rows)

        df = df[
            [
                "molecule_id",
                "total_conformers",
                "valid_energy_conformers",
                "removed_missing_energy",
                "kept_after_pruning",
                "strategy",
            ]
        ]

        with AtomicWriter(summary_path) as f:
            df.to_csv(f, index=False)

        self.log(f"Pruned energies written to: {energies_path}")
        self.log(f"Pruning summary written to: {summary_path}")
