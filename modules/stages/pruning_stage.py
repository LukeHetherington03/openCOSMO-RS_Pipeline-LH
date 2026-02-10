import os
import pandas as pd
import numpy as np

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter
from modules.utils.conformers import ConformerRecord, ConformerSet


class PruningStage(BaseStage):
    """
    Prunes conformers for each molecule.

    Input:
        stage_input = energies.json from Generation or Optimisation

    Output:
        energies.json (canonical pruned conformer set)
        pruning_summary.csv
    """

    def execute(self):
        # Strict mode
        self.strict_mode = self.strict("pruning")

        # Declare canonical output
        self.set_stage_output("energies.json")

        # Retrieve stage input
        energies_file = self.require_file(
            self.get_stage_input(),
            "stage_input energies.json"
        )

        self.log(f"Stage input: {energies_file}")
        self.log(f"Strict mode: {self.strict_mode}")

        # Load conformers
        conformer_set = ConformerSet.load(energies_file)
        groups = conformer_set.group_by_molecule()

        molecule_ids = list(groups.keys())
        self.set_items(molecule_ids)

        self.warnings = {
            "molecules_all_missing_energy": 0,
        }

        pruned_records = []
        summary_rows = []

        for inchi_key, confs in groups.items():
            try:
                survivors, summary_row = self._prune_group(inchi_key, confs, self.parameters)

                for rec in survivors:
                    pruned_records.append(rec)

                summary_rows.append(summary_row)
                self.update_progress(inchi_key, success=True)

                self.log(
                    f"[INFO] Pruning complete for {inchi_key}: "
                    f"kept {len(survivors)} conformers"
                )

            except Exception as e:
                self.log(f"[ERROR] Pruning failed for {inchi_key}: {e}")
                self.update_progress(inchi_key, success=False)
                if self.strict_mode:
                    self.fail(f"Pruning failed for {inchi_key}: {e}")

        pruned_set = ConformerSet(records=pruned_records)
        self._write_outputs(pruned_set, summary_rows)

        self._log_warning_summary()

    # ------------------------------------------------------------
    # Pruning logic for a single molecule
    # ------------------------------------------------------------
    def _prune_group(self, inchi_key, conformers, params):
        total_confs = len(conformers)

        # 1) Remove missing-energy conformers
        valid = []
        removed_missing_energy = 0

        for rec in conformers:
            e = rec.energy
            if e is None or (isinstance(e, float) and np.isnan(e)):
                removed_missing_energy += 1
                self.log(
                    f"[WARNING] {inchi_key}: conformer {rec.lookup_id} removed (missing energy)"
                )
            else:
                valid.append(rec)

        if not valid:
            msg = f"All conformers for {inchi_key} missing energies."
            self.log(f"[ERROR] {msg}")
            self.warnings["molecules_all_missing_energy"] += 1
            if self.strict_mode:
                raise RuntimeError(msg)
            return [], {
                "inchi_key": inchi_key,
                "total_conformers": total_confs,
                "valid_energy_conformers": 0,
                "removed_missing_energy": removed_missing_energy,
                "kept_after_pruning": 0,
                **self._extract_params(params),
            }

        survivors = list(valid)

        # Deterministic pruning pipeline
        p = self._extract_params(params)

        if p["rmsd_threshold"] is not None:
            survivors = self._prune_rmsd(survivors, p["rmsd_threshold"])

        if p["energy_window"] is not None:
            survivors = self._prune_energy_window(survivors, p["energy_window"])

        if p["max_energy"] is not None:
            survivors = self._prune_max_energy(survivors, p["max_energy"])

        if p["percentile"] is not None:
            survivors = self._prune_percentile(survivors, p["percentile"])

        if p["n"] is not None:
            survivors = self._prune_keep_lowest_n(survivors, p["n"])

        if p["n_high"] is not None:
            survivors = self._prune_keep_highest_n(survivors, p["n_high"])

        if p["n_start"] is not None and p["n"] is not None:
            survivors = self._prune_slice(survivors, p["n_start"], p["n"])

        summary_row = {
            "inchi_key": inchi_key,
            "total_conformers": total_confs,
            "valid_energy_conformers": len(valid),
            "removed_missing_energy": removed_missing_energy,
            "kept_after_pruning": len(survivors),
            **p,
        }

        return survivors, summary_row

    # ------------------------------------------------------------
    # Parameter extraction
    # ------------------------------------------------------------
    def _extract_params(self, params):
        return {
            "rmsd_threshold": self._get_float(params, "rmsd_threshold"),
            "energy_window": self._get_float(params, "energy_window"),
            "max_energy": self._get_float(params, "max_energy"),
            "percentile": self._get_float(params, "percentile"),
            "n": self._get_int(params, "n"),
            "n_high": self._get_int(params, "n_high"),
            "n_start": self._get_int(params, "n_start"),
        }

    # ------------------------------------------------------------
    # Strategy helpers (unchanged)
    # ------------------------------------------------------------
    def _get_float(self, params, key):
        if key not in params:
            return None
        try:
            return float(params[key])
        except Exception:
            self.log(f"[WARNING] Invalid float for {key}: {params[key]} (ignored)")
            return None

    def _get_int(self, params, key):
        if key not in params:
            return None
        try:
            return int(params[key])
        except Exception:
            self.log(f"[WARNING] Invalid int for {key}: {params[key]} (ignored)")
            return None

    # ------------------------------------------------------------
    # Pruning strategies (unchanged)
    # ------------------------------------------------------------
    def _prune_rmsd(self, conformers, rmsd_threshold):
        self.log(
            f"[WARNING] RMSD pruning (rmsd_threshold={rmsd_threshold}) not implemented; keeping all conformers"
        )
        return conformers

    def _prune_energy_window(self, conformers, energy_window):
        if not conformers:
            return conformers
        energies = np.array([c.energy for c in conformers], dtype=float)
        min_e = float(energies.min())
        survivors = [c for c in conformers if (c.energy - min_e) <= energy_window]
        self.log(
            f"[INFO] Energy window pruning (ΔE <= {energy_window}): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_max_energy(self, conformers, max_energy):
        survivors = [c for c in conformers if c.energy <= max_energy]
        self.log(
            f"[INFO] Max energy pruning (E <= {max_energy}): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_percentile(self, conformers, percentile):
        if not conformers:
            return conformers
        energies = np.array([c.energy for c in conformers], dtype=float)
        cutoff = float(np.percentile(energies, percentile))
        survivors = [c for c in conformers if c.energy <= cutoff]
        self.log(
            f"[INFO] Percentile pruning (<= {percentile}th, cutoff={cutoff:.4f}): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_keep_lowest_n(self, conformers, n):
        if n <= 0 or not conformers:
            return []
        sorted_conf = sorted(conformers, key=lambda c: c.energy)
        survivors = sorted_conf[:n]
        self.log(
            f"[INFO] Keep lowest N pruning (N={n}): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_keep_highest_n(self, conformers, n_high):
        if n_high <= 0 or not conformers:
            return []
        sorted_conf = sorted(conformers, key=lambda c: c.energy, reverse=True)
        survivors = sorted_conf[:n_high]
        self.log(
            f"[INFO] Keep highest N pruning (N={n_high}): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_slice(self, conformers, n_start, n):
        if n <= 0 or not conformers:
            return []
        total = len(conformers)

        if n_start < 0:
            start = total + n_start
        else:
            start = n_start

        end = start + n
        survivors = conformers[max(0, start):max(0, end)]
        self.log(
            f"[INFO] Slice pruning (start={n_start}, n={n}): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    # ------------------------------------------------------------
    # Output writing
    # ------------------------------------------------------------
    def _write_outputs(self, pruned_set: ConformerSet, summary_rows):
        outputs_dir = self.outputs_dir

        # Canonical output
        energies_path = self.get_stage_output()
        pruned_set.save(energies_path)

        # Summary CSV
        summary_path = os.path.join(outputs_dir, "pruning_summary.csv")
        df = pd.DataFrame(summary_rows)

        cols = [
            "inchi_key",
            "total_conformers",
            "valid_energy_conformers",
            "removed_missing_energy",
            "kept_after_pruning",
            "rmsd_threshold",
            "energy_window",
            "max_energy",
            "percentile",
            "n",
            "n_high",
            "n_start",
        ]
        df = df[cols]

        with AtomicWriter(summary_path) as f:
            df.to_csv(f, index=False)

        self.log(f"Pruned energies written to: {energies_path}")
        self.log(f"Pruning summary written to: {summary_path}")

    # ------------------------------------------------------------
    # Warning summary
    # ------------------------------------------------------------
    def _log_warning_summary(self):
        total = sum(self.warnings.values())
        if total == 0:
            self.log("Pruning completed with no warnings.")
            return

        self.log("Pruning completed with warnings:")
        if self.warnings["molecules_all_missing_energy"] > 0:
            self.log(
                " - Molecules with all conformers missing energy: "
                f"{self.warnings['molecules_all_missing_energy']}"
            )
