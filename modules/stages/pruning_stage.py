import os
import shutil
import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdMolAlign

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter
from modules.utils.conformers import ConformerRecord, ConformerSet


# Energy conversion constants
HARTREE_TO_KJ = 2625.49962
HARTREE_TO_KCAL = 627.509474


class PruningStage(BaseStage):
    """
    Prunes conformers for each molecule.

    Input:
        stage_input = energies.json from Generation or Optimisation

    Output:
        energies.json (canonical pruned conformer set)
        pruning_summary.csv
    """

    # -------------------------------------------------------------------------
    # Dynamic pruning registry
    # -------------------------------------------------------------------------
    PRUNING_METHODS = {
        "rmsd_threshold": "_prune_rmsd",
        "energy_window": "_prune_energy_window",
        "max_energy": "_prune_max_energy",
        "percentile": "_prune_percentile",
        "n": "_prune_keep_lowest_n",
        "n_high": "_prune_keep_highest_n",
        "n_start": "_prune_slice",
    }

    # -------------------------------------------------------------------------
    # Execution
    # -------------------------------------------------------------------------
    def execute(self):
        self.strict_mode = self.strict("pruning")
        self.set_stage_output("energies.json")

        energies_file = self.require_file(
            self.get_stage_input(),
            "stage_input energies.json"
        )

        # Copy input energies.json → inputs/
        inputs_energies = os.path.join(self.inputs_dir, "energies.json")
        if os.path.abspath(energies_file) != os.path.abspath(inputs_energies):
            shutil.copy(energies_file, inputs_energies)
        self.log(f"Copied energies.json → {inputs_energies}")

        self.log(f"Stage input: {energies_file}")
        self.log(f"Strict mode: {self.strict_mode}")

        conformer_set = ConformerSet.load(energies_file)
        groups = conformer_set.group_by_molecule()

        molecule_ids = list(groups.keys())
        self.set_items(molecule_ids)

        self.warnings = {
            "molecules_all_missing_energy": 0,
            "molecules_all_pruned": 0,
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
                    f"[INFO] Molecule {inchi_key}: {len(confs)} → {len(survivors)} conformers after pruning"
                )

            except Exception as e:
                self.log(f"[ERROR] Pruning failed for {inchi_key}: {e}")
                self.update_progress(inchi_key, success=False)
                if self.strict_mode:
                    self.fail(f"Pruning failed for {inchi_key}: {e}")

        pruned_set = ConformerSet(records=pruned_records)
        self._write_outputs(pruned_set, summary_rows)
        self._log_warning_summary()

    # -------------------------------------------------------------------------
    # Parameter extraction helpers
    # -------------------------------------------------------------------------
    def _get_float(self, params, key):
        val = params.get(key, None)
        if val is None:
            return None
        try:
            return float(val)
        except Exception:
            self.log(f"[WARNING] Could not parse float for '{key}': {val}")
            return None

    def _get_int(self, params, key):
        val = params.get(key, None)
        if val is None:
            return None
        try:
            return int(val)
        except Exception:
            self.log(f"[WARNING] Could not parse int for '{key}': {val}")
            return None

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

    # -------------------------------------------------------------------------
    # Pruning logic for a single molecule
    # -------------------------------------------------------------------------
    def _prune_group(self, inchi_key, conformers, params):
        total_confs = len(conformers)
        p = self._extract_params(params)
        keep_all = params.get("keep_all", False)

        # Fast path
        if keep_all:
            self.log(f"[INFO] keep_all=True → skipping pruning for {inchi_key}")
            return conformers, {
                "inchi_key": inchi_key,
                "total_conformers": total_confs,
                "valid_energy_conformers": total_confs,
                "removed_missing_energy": 0,
                "kept_after_pruning": total_confs,
                "methods_used": "none",
                **p,
            }

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
                "methods_used": "none",
                **p,
            }

        survivors = list(valid)
        methods_used = []

        # ---------------------------------------------------------------------
        # Dynamic pruning dispatch
        # ---------------------------------------------------------------------
        for param_name, method_name in self.PRUNING_METHODS.items():
            value = p.get(param_name)
            if value is None:
                continue

            prune_fn = getattr(self, method_name)
            before = len(survivors)

            # Slice requires two parameters
            if param_name == "n_start":
                survivors = prune_fn(survivors, p["n_start"], p["n"])
            else:
                survivors = prune_fn(survivors, value)

            after = len(survivors)
            methods_used.append(f"{param_name} ({before}→{after})")

        # Ensure survivors exist
        if not survivors:
            self.warnings["molecules_all_pruned"] += 1
            msg = f"All conformers pruned for {inchi_key}"
            self.log(f"[WARNING] {msg}")
            if self.strict_mode:
                raise RuntimeError(msg)

        survivors = sorted(survivors, key=lambda c: (c.energy, c.lookup_id))

        summary_row = {
            "inchi_key": inchi_key,
            "total_conformers": total_confs,
            "valid_energy_conformers": len(valid),
            "removed_missing_energy": removed_missing_energy,
            "kept_after_pruning": len(survivors),
            "methods_used": "; ".join(methods_used) if methods_used else "none",
            **p,
        }

        return survivors, summary_row

    # -------------------------------------------------------------------------
    # Energy conversion
    # -------------------------------------------------------------------------
    def _convert_energy(self, energy_hartree, units):
        if energy_hartree is None:
            return None

        units = units.lower()

        if units == "hartree":
            return energy_hartree
        if units == "kcal":
            return energy_hartree * HARTREE_TO_KCAL
        if units == "kJ" or units == "kj":
            return energy_hartree * HARTREE_TO_KJ

        self.log(f"[WARNING] Unknown energy units '{units}', assuming hartree")
        return energy_hartree

    # -------------------------------------------------------------------------
    # Pruning methods
    # -------------------------------------------------------------------------
    def _prune_rmsd(self, conformers, rmsd_threshold):
        if not conformers or rmsd_threshold <= 0:
            return conformers

        sorted_conf = sorted(conformers, key=lambda c: c.energy)
        survivors = []
        survivor_mols = []

        for rec in sorted_conf:
            try:
                mol = Chem.MolFromXYZFile(rec.xyz_path)
                if mol is None:
                    self.log(f"[WARNING] Could not load XYZ for RMSD pruning: {rec.lookup_id}")
                    continue
            except Exception as e:
                self.log(f"[WARNING] Failed to parse XYZ for {rec.lookup_id}: {e}")
                continue

            too_close = False
            for sm in survivor_mols:
                rms = rdMolAlign.GetBestRMS(mol, sm)
                if rms < rmsd_threshold:
                    too_close = True
                    break

            if not too_close:
                survivors.append(rec)
                survivor_mols.append(mol)

        self.log(
            f"[INFO] RMSD pruning (threshold={rmsd_threshold} Å): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_energy_window(self, conformers, energy_window, units="kcal"):
        if not conformers:
            return conformers

        converted = np.array([
            self._convert_energy(c.energy, units) for c in conformers
        ], dtype=float)

        min_e = float(converted.min())
        survivors = [
            c for c, e in zip(conformers, converted)
            if (e - min_e) <= energy_window
        ]

        self.log(
            f"[INFO] Energy window pruning (ΔE <= {energy_window} {units}): "
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

    # -------------------------------------------------------------------------
    # Output writing
    # -------------------------------------------------------------------------
    def _write_outputs(self, pruned_set: ConformerSet, summary_rows):
        outputs_dir = self.outputs_dir

        # Canonical output
        energies_path = self.get_stage_output()
        pruned_set.save(energies_path)

        # Summary CSV
        summary_path = os.path.join(outputs_dir, "pruning_summary.csv")
        df = pd.DataFrame(summary_rows)

        # Human-readable column names
        df = df.rename(columns={
            "inchi_key": "Molecule",
            "total_conformers": "Total Conformers",
            "valid_energy_conformers": "With Valid Energy",
            "removed_missing_energy": "Removed (Missing Energy)",
            "kept_after_pruning": "Final Count",
            "methods_used": "Pruning Steps Applied",
            "rmsd_threshold": "RMSD Threshold (Å)",
            "energy_window": "Energy Window",
            "max_energy": "Max Energy",
            "percentile": "Percentile Cutoff",
            "n": "Keep Lowest N",
            "n_high": "Keep Highest N",
            "n_start": "Slice Start",
        })

        cols = [
            "Molecule",
            "Total Conformers",
            "With Valid Energy",
            "Removed (Missing Energy)",
            "Final Count",
            "Pruning Steps Applied",
            "RMSD Threshold (Å)",
            "Energy Window",
            "Max Energy",
            "Percentile Cutoff",
            "Keep Lowest N",
            "Keep Highest N",
            "Slice Start",
        ]

        df = df.reindex(columns=cols)

        with AtomicWriter(summary_path) as f:
            df.to_csv(f, index=False)

        self.log(f"Pruned energies written to: {energies_path}")
        self.log(f"Pruning summary written to: {summary_path}")

    # -------------------------------------------------------------------------
    # Warning summary
    # -------------------------------------------------------------------------
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
        if self.warnings["molecules_all_pruned"] > 0:
            self.log(
                " - Molecules with all conformers pruned: "
                f"{self.warnings['molecules_all_pruned']}"
            )
