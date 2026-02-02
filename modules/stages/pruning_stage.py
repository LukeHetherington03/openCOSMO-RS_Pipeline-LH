import os
import pandas as pd
import numpy as np

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter
from modules.utils.conformers import ConformerRecord, ConformerSet


class PruningStage(BaseStage):
    """
    Prunes conformers for each molecule.

    Features:
      - Uses ConformerSet / ConformerRecord
      - Reads and writes energies.json (canonical conformer store)
      - Groups by inchi_key (molecule-level)
      - Deterministic pruning pipeline driven by args:
            * rmsd_threshold
            * energy_window
            * max_energy
            * percentile
            * n        (keep lowest N)
            * n_high   (keep highest N)
            * n_start  (slice start, supports negative)
      - Always removes conformers with missing energies
      - Produces:
            * pruned energies.json
            * pruning_summary.csv
    """

    # ------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------
    def execute(self):
        self.log_header("Starting Pruning Stage")

        cfg = self.config or {}
        self.strict_mode = bool(cfg.get("pruning", {}).get("strict", False))

        params = self.parameters or {}

        energies_file = params.get("energies_file")
        if energies_file is None:
            energies_file = self._auto_detect_energies()

        if not os.path.isfile(energies_file):
            self.fail(f"PruningStage: energies_file does not exist: {energies_file}")

        self.log(f"Energies file: {energies_file}")
        self.log(f"Strict mode: {self.strict_mode}")

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
                survivors, summary_row = self._prune_group(inchi_key, confs, params)
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
        self.log_header("Pruning Stage Complete")
        self.job.mark_complete()

    # ------------------------------------------------------------
    # Input handling
    # ------------------------------------------------------------
    def _auto_detect_energies(self) -> str:
        candidates = [
            os.path.join(self.inputs_dir, "energies.json"),
            os.path.join(self.outputs_dir, "energies.json"),
        ]
        for path in candidates:
            if os.path.isfile(path):
                self.log(f"Auto-detected energies.json at: {path}")
                return path
        self.fail(
            "PruningStage: energies_file not provided and energies.json "
            "could not be auto-detected."
        )

    # ------------------------------------------------------------
    # Pruning logic for a single molecule
    # ------------------------------------------------------------
    def _prune_group(self, inchi_key, conformers, params):
        """
        Returns:
            survivors: list[ConformerRecord]
            summary_row: dict for pruning_summary.csv
        """

        total_confs = len(conformers)

        # 1) Remove conformers with missing or NaN energies
        valid = []
        removed_missing_energy = 0

        for rec in conformers:
            e = rec.energy
            if e is None or (isinstance(e, float) and np.isnan(e)):
                removed_missing_energy += 1
                self.log(
                    f"[WARNING] {inchi_key}: conformer {rec.lookup_id} "
                    "removed (missing energy)"
                )
            else:
                valid.append(rec)

        if not valid:
            msg = f"All conformers for {inchi_key} missing energies."
            self.log(f"[ERROR] {msg}")
            self.warnings["molecules_all_missing_energy"] += 1
            if self.strict_mode:
                raise RuntimeError(msg)
            # Non-strict: keep original list (no pruning possible)
            valid = []

        survivors = list(valid)

        # If nothing valid, build summary and return early
        if not survivors:
            summary_row = {
                "inchi_key": inchi_key,
                "total_conformers": total_confs,
                "valid_energy_conformers": 0,
                "removed_missing_energy": removed_missing_energy,
                "kept_after_pruning": 0,
                "rmsd_threshold": params.get("rmsd_threshold"),
                "energy_window": params.get("energy_window"),
                "max_energy": params.get("max_energy"),
                "percentile": params.get("percentile"),
                "n": params.get("n"),
                "n_high": params.get("n_high"),
                "n_start": params.get("n_start"),
            }
            return [], summary_row

        # Deterministic pruning pipeline
        rmsd_threshold = self._get_float(params, "rmsd_threshold")
        energy_window = self._get_float(params, "energy_window")
        max_energy = self._get_float(params, "max_energy")
        percentile = self._get_float(params, "percentile")
        n = self._get_int(params, "n")
        n_high = self._get_int(params, "n_high")
        n_start = self._get_int(params, "n_start")

        # 2) RMSD clustering (placeholder: currently no-op)
        if rmsd_threshold is not None:
            survivors = self._prune_rmsd(survivors, rmsd_threshold)

        # 3) Energy window
        if energy_window is not None:
            survivors = self._prune_energy_window(survivors, energy_window)

        # 4) Max energy cutoff
        if max_energy is not None:
            survivors = self._prune_max_energy(survivors, max_energy)

        # 5) Percentile filter
        if percentile is not None:
            survivors = self._prune_percentile(survivors, percentile)

        # 6) Keep lowest N
        if n is not None:
            survivors = self._prune_keep_lowest_n(survivors, n)

        # 7) Keep highest N
        if n_high is not None:
            survivors = self._prune_keep_highest_n(survivors, n_high)

        # 8) Slice (index-based)
        if n_start is not None and n is not None:
            survivors = self._prune_slice(survivors, n_start, n)

        summary_row = {
            "inchi_key": inchi_key,
            "total_conformers": total_confs,
            "valid_energy_conformers": len(valid),
            "removed_missing_energy": removed_missing_energy,
            "kept_after_pruning": len(survivors),
            "rmsd_threshold": rmsd_threshold,
            "energy_window": energy_window,
            "max_energy": max_energy,
            "percentile": percentile,
            "n": n,
            "n_high": n_high,
            "n_start": n_start,
        }

        return survivors, summary_row

    # ------------------------------------------------------------
    # Strategy helpers: parsing
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
    # Strategy: RMSD clustering (placeholder)
    # ------------------------------------------------------------
    def _prune_rmsd(self, conformers, rmsd_threshold):
        self.log(
            f"[WARNING] RMSD pruning (rmsd_threshold={rmsd_threshold}) "
            "not implemented; keeping all conformers"
        )
        return conformers

    # ------------------------------------------------------------
    # Strategy: keep all within ΔE of minimum
    # ------------------------------------------------------------
    def _prune_energy_window(self, conformers, energy_window):
        if not conformers:
            return conformers
        energies = np.array([c.energy for c in conformers], dtype=float)
        min_e = float(energies.min())
        survivors = [
            c for c in conformers
            if (c.energy - min_e) <= energy_window
        ]
        self.log(
            f"[INFO] Energy window pruning (ΔE <= {energy_window}): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    # ------------------------------------------------------------
    # Strategy: max energy cutoff
    # ------------------------------------------------------------
    def _prune_max_energy(self, conformers, max_energy):
        survivors = [c for c in conformers if c.energy <= max_energy]
        self.log(
            f"[INFO] Max energy pruning (E <= {max_energy}): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    # ------------------------------------------------------------
    # Strategy: percentile filter
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # Strategy: keep lowest N energies
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # Strategy: keep highest N energies
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # Strategy: slice by index (supports negative n_start)
    # ------------------------------------------------------------
    def _prune_slice(self, conformers, n_start, n):
        if n <= 0 or not conformers:
            return []
        total = len(conformers)

        if n_start < 0:
            start = total + n_start  # e.g. -1 → last index
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

        energies_path = os.path.join(outputs_dir, "energies.json")
        pruned_set.save(energies_path)

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
