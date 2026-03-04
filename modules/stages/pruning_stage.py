#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modules/stages/pruning_stage.py

Conformer pruning stage.

=========================================================================
INPUT
=========================================================================

  stage_input : energies.json from GenerationStage or OptimisationStage
    ConformerSet JSON — one entry per conformer.

=========================================================================
OUTPUT  (canonical stage output)
=========================================================================

  energies.json   — pruned ConformerSet JSON
  pruning_summary.csv

=========================================================================
PRUNING METHODS  (applied in order if configured)
=========================================================================

  rmsd_threshold  float (Å)    Remove structurally similar conformers.
                                Sorted by energy first; lower-energy
                                survivor kept.
  energy_window   float        Remove conformers whose energy is more than
                                this value above the minimum.  Units from
                                config (default kcal/mol if hartree input).
  max_energy      float        Hard absolute energy ceiling.
  percentile      float (0–100) Remove conformers above this energy percentile.
  n               int          Keep only the N lowest-energy conformers.
  n_high          int          Keep only the N highest-energy conformers.
  n_start + n     int, int     Slice: keep n conformers starting at n_start.

  keep_all        bool         Skip all pruning — pass through unchanged.

Methods are applied sequentially in the order listed in PRUNING_METHODS.

=========================================================================
ITEM TRACKING
=========================================================================

  Items are inchi_keys — one per unique molecule.
  update_progress(inchi_key) called after each molecule.
  Serial execution — no run_parallel().

=========================================================================
STRICT MODE
=========================================================================

  self.strict("pruning") reads config["pruning"]["strict"].
  When True, any per-molecule failure raises and aborts the stage.

"""

import os
import shutil

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdMolAlign

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter
from modules.utils.conformers import ConformerRecord, ConformerSet

HARTREE_TO_KJ   = 2625.49962
HARTREE_TO_KCAL = 627.509474


class PruningStage(BaseStage):
    """
    Conformer pruning stage.

    Applies one or more pruning methods sequentially to each molecule's
    conformer pool.  Serial execution.

    Stages never call mark_complete() — BaseStage.run() owns lifecycle.
    """

    # Ordered pruning dispatch table
    PRUNING_METHODS = {
        "rmsd_threshold": "_prune_rmsd",
        "energy_window":  "_prune_energy_window",
        "max_energy":     "_prune_max_energy",
        "percentile":     "_prune_percentile",
        "n":              "_prune_keep_lowest_n",
        "n_high":         "_prune_keep_highest_n",
        "n_start":        "_prune_slice",
    }

    # =========================================================================
    # Entry point
    # =========================================================================

    def execute(self):
        self.strict_mode = self.strict("pruning")
        self.set_stage_output("energies.json")

        energies_file = self.require_file(
            self.get_stage_input(), "stage_input energies.json"
        )

        # Copy input into inputs/ for reproducibility
        inputs_energies = os.path.join(self.inputs_dir, "energies.json")
        os.makedirs(self.inputs_dir, exist_ok=True)
        if os.path.abspath(energies_file) != os.path.abspath(inputs_energies):
            shutil.copy(energies_file, inputs_energies)
        self.log_info(f"Stage input: {energies_file}")
        self.log_config(f"strict={self.strict_mode}")

        conformer_set = ConformerSet.load(energies_file)
        groups        = conformer_set.group_by_molecule()

        molecule_ids = list(groups.keys())
        self.set_items(molecule_ids)

        self.warnings = {
            "molecules_all_missing_energy": 0,
            "molecules_all_pruned":         0,
        }

        pruned_records = []
        summary_rows   = []

        for inchi_key, confs in groups.items():
            try:
                survivors, summary_row = self._prune_group(
                    inchi_key, confs, self.parameters
                )
                for rec in survivors:
                    pruned_records.append(rec)
                summary_rows.append(summary_row)
                self.log_info(
                    f"{inchi_key}: {len(confs)} → {len(survivors)} "
                    f"conformers after pruning"
                )
                self.update_progress(inchi_key, success=True)

            except Exception as e:
                self.log_error(f"Pruning failed for {inchi_key}: {e}")
                self.update_progress(inchi_key, success=False)
                if self.strict_mode:
                    self.fail(f"Pruning failed for {inchi_key}: {e}")

        self._write_outputs(ConformerSet(records=pruned_records), summary_rows)
        self._log_warning_summary()

    # =========================================================================
    # Pruning dispatch for one molecule
    # =========================================================================

    def _prune_group(
        self, inchi_key: str, conformers: list, params: dict
    ) -> tuple[list, dict]:
        total    = len(conformers)
        p        = self._extract_params(params)
        keep_all = params.get("keep_all", False)

        if keep_all:
            self.log_info(f"{inchi_key}: keep_all=True — pruning skipped")
            return conformers, self._summary_row(
                inchi_key, total, total, 0, total, "none", p
            )

        # ── Remove conformers with missing energy ─────────────────────────────
        valid                 = []
        removed_missing_energy = 0

        for rec in conformers:
            e = rec.energy
            if e is None or (isinstance(e, float) and np.isnan(e)):
                removed_missing_energy += 1
                self.log_warning(
                    f"{inchi_key}: {rec.lookup_id} removed — missing energy"
                )
            else:
                valid.append(rec)

        if not valid:
            self.log_error(f"{inchi_key}: all conformers have missing energy")
            self.warnings["molecules_all_missing_energy"] += 1
            if self.strict_mode:
                raise RuntimeError(
                    f"All conformers for {inchi_key} have missing energy"
                )
            return [], self._summary_row(
                inchi_key, total, 0, removed_missing_energy, 0, "none", p
            )

        # ── Apply pruning methods in order ────────────────────────────────────
        survivors    = list(valid)
        methods_used = []

        for param_name, method_name in self.PRUNING_METHODS.items():
            value = p.get(param_name)
            if value is None:
                continue

            prune_fn = getattr(self, method_name)
            before   = len(survivors)

            if param_name == "n_start":
                survivors = prune_fn(survivors, p["n_start"], p.get("n", 0))
            else:
                survivors = prune_fn(survivors, value)

            after = len(survivors)
            methods_used.append(f"{param_name} ({before}→{after})")

        if not survivors:
            self.warnings["molecules_all_pruned"] += 1
            msg = f"All conformers pruned for {inchi_key}"
            self.log_warning(msg)
            if self.strict_mode:
                raise RuntimeError(msg)

        survivors = sorted(survivors, key=lambda c: (c.energy, c.lookup_id))

        return survivors, self._summary_row(
            inchi_key, total, len(valid), removed_missing_energy,
            len(survivors),
            "; ".join(methods_used) if methods_used else "none",
            p,
        )

    @staticmethod
    def _summary_row(
        inchi_key, total, n_valid, removed_missing, kept, methods, p
    ) -> dict:
        return {
            "inchi_key":                 inchi_key,
            "total_conformers":          total,
            "valid_energy_conformers":   n_valid,
            "removed_missing_energy":    removed_missing,
            "kept_after_pruning":        kept,
            "methods_used":              methods,
            **p,
        }

    # =========================================================================
    # Parameter extraction
    # =========================================================================

    def _extract_params(self, params: dict) -> dict:
        def _f(k):
            v = params.get(k)
            if v is None:
                return None
            try:
                return float(v)
            except Exception:
                self.log_warning(f"Could not parse float for '{k}': {v}")
                return None

        def _i(k):
            v = params.get(k)
            if v is None:
                return None
            try:
                return int(v)
            except Exception:
                self.log_warning(f"Could not parse int for '{k}': {v}")
                return None

        return {
            "rmsd_threshold": _f("rmsd_threshold"),
            "energy_window":  _f("energy_window"),
            "max_energy":     _f("max_energy"),
            "percentile":     _f("percentile"),
            "n":              _i("n"),
            "n_high":         _i("n_high"),
            "n_start":        _i("n_start"),
        }

    # =========================================================================
    # Energy conversion
    # =========================================================================

    def _convert_energy(self, energy_hartree: float, units: str) -> float | None:
        if energy_hartree is None:
            return None
        u = str(units).lower()
        if u == "hartree":
            return energy_hartree
        if u == "kcal":
            return energy_hartree * HARTREE_TO_KCAL
        if u in ("kj", "kJ"):
            return energy_hartree * HARTREE_TO_KJ
        self.log_warning(f"Unknown energy units '{units}' — assuming hartree")
        return energy_hartree

    # =========================================================================
    # Pruning methods
    # =========================================================================

    def _prune_rmsd(self, conformers: list, rmsd_threshold: float) -> list:
        if not conformers or rmsd_threshold <= 0:
            return conformers

        sorted_conf  = sorted(conformers, key=lambda c: c.energy)
        survivors    = []
        survivor_mols = []

        for rec in sorted_conf:
            try:
                mol = Chem.MolFromXYZFile(rec.xyz_path)
                if mol is None:
                    self.log_warning(
                        f"RMSD pruning: could not load XYZ for {rec.lookup_id}"
                    )
                    continue
            except Exception as e:
                self.log_warning(
                    f"RMSD pruning: XYZ parse failed for {rec.lookup_id} — {e}"
                )
                continue

            too_close = any(
                rdMolAlign.GetBestRMS(mol, sm) < rmsd_threshold
                for sm in survivor_mols
            )
            if not too_close:
                survivors.append(rec)
                survivor_mols.append(mol)

        self.log_info(
            f"RMSD pruning (threshold={rmsd_threshold} Å): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_energy_window(
        self, conformers: list, energy_window: float, units: str = "kcal"
    ) -> list:
        if not conformers:
            return conformers

        converted = np.array(
            [self._convert_energy(c.energy, units) for c in conformers],
            dtype=float,
        )
        min_e     = float(converted.min())
        survivors = [
            c for c, e in zip(conformers, converted)
            if (e - min_e) <= energy_window
        ]
        self.log_info(
            f"Energy window (ΔE ≤ {energy_window} {units}): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_max_energy(self, conformers: list, max_energy: float) -> list:
        survivors = [c for c in conformers if c.energy <= max_energy]
        self.log_info(
            f"Max energy (E ≤ {max_energy}): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_percentile(self, conformers: list, percentile: float) -> list:
        if not conformers:
            return conformers
        energies = np.array([c.energy for c in conformers], dtype=float)
        cutoff   = float(np.percentile(energies, percentile))
        survivors = [c for c in conformers if c.energy <= cutoff]
        self.log_info(
            f"Percentile (≤ {percentile}th, cutoff={cutoff:.4f}): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_keep_lowest_n(self, conformers: list, n: int) -> list:
        if n <= 0 or not conformers:
            return []
        survivors = sorted(conformers, key=lambda c: c.energy)[:n]
        self.log_info(
            f"Keep lowest N={n}: {len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_keep_highest_n(self, conformers: list, n_high: int) -> list:
        if n_high <= 0 or not conformers:
            return []
        survivors = sorted(conformers, key=lambda c: c.energy, reverse=True)[:n_high]
        self.log_info(
            f"Keep highest N={n_high}: {len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_slice(self, conformers: list, n_start: int, n: int) -> list:
        if n <= 0 or not conformers:
            return []
        total = len(conformers)
        start = (total + n_start) if n_start < 0 else n_start
        survivors = conformers[max(0, start): max(0, start + n)]
        self.log_info(
            f"Slice (start={n_start}, n={n}): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    # =========================================================================
    # Output writing
    # =========================================================================

    def _write_outputs(self, pruned_set: ConformerSet, summary_rows: list):
        # energies.json
        energies_path = self.get_stage_output()
        pruned_set.save(energies_path)

        # pruning_summary.csv
        summary_path = os.path.join(self.outputs_dir, "pruning_summary.csv")
        df = pd.DataFrame(summary_rows).rename(columns={
            "inchi_key":               "Molecule",
            "total_conformers":        "Total Conformers",
            "valid_energy_conformers": "With Valid Energy",
            "removed_missing_energy":  "Removed (Missing Energy)",
            "kept_after_pruning":      "Final Count",
            "methods_used":            "Pruning Steps Applied",
            "rmsd_threshold":          "RMSD Threshold (Å)",
            "energy_window":           "Energy Window",
            "max_energy":              "Max Energy",
            "percentile":              "Percentile Cutoff",
            "n":                       "Keep Lowest N",
            "n_high":                  "Keep Highest N",
            "n_start":                 "Slice Start",
        })

        col_order = [
            "Molecule", "Total Conformers", "With Valid Energy",
            "Removed (Missing Energy)", "Final Count",
            "Pruning Steps Applied", "RMSD Threshold (Å)",
            "Energy Window", "Max Energy", "Percentile Cutoff",
            "Keep Lowest N", "Keep Highest N", "Slice Start",
        ]
        df = df.reindex(columns=col_order)

        with AtomicWriter(summary_path) as f:
            df.to_csv(f, index=False)

        self.log_info(
            f"Outputs: energies.json ({len(pruned_set)} conformers), "
            f"pruning_summary.csv ({len(summary_rows)} molecules)"
        )

    # =========================================================================
    # Warning summary
    # =========================================================================

    def _log_warning_summary(self):
        total = sum(self.warnings.values())
        if total == 0:
            self.log_info("Pruning completed with no warnings.")
            return

        self.log_info("Pruning completed with warnings:")
        if self.warnings["molecules_all_missing_energy"]:
            self.log_warning(
                f"Molecules with all conformers missing energy: "
                f"{self.warnings['molecules_all_missing_energy']}"
            )
        if self.warnings["molecules_all_pruned"]:
            self.log_warning(
                f"Molecules with all conformers pruned: "
                f"{self.warnings['molecules_all_pruned']}"
            )