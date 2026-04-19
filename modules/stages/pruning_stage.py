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
  rdkit_post_opt_legacy  bool  Replicate the 4 internal filters from the IC
                                colleague pipeline's calculate_rdkit(), applied
                                in order: convergence → overlapping atoms →
                                Boltzmann (T=298.15 K, cutoff 1%) → RMSD 1.0 Å.
                                Use after forcefield_mmff optimisation.

  window_pair     float (kcal/mol)
                               Testing/paper tool — keeps exactly two conformers:
                                (1) the global energy minimum, and (2) the
                                lowest-energy conformer whose energy lies above
                                the given window.  Degrades to one conformer if
                                nothing exists outside the window.  Not intended
                                for normal pipeline use.

  select_conf     str|list     Keep all conformers for molecules whose inchi_key
                                is in the given list; drop all conformers for
                                every other molecule.  Accepts a single inchi_key
                                string or a list of strings.

  remove_conf     str|list     Drop all conformers for molecules whose inchi_key
                                is in the given list; keep all others unchanged.
                                Accepts a single inchi_key string or a list.

  conf_select     list[int]    Pick specific conformers by 0-based index from
                                the current pool, e.g. [0, 2, 4].  Out-of-range
                                indices are logged and skipped (not an error).

  keep_all        bool         Skip all pruning — pass through unchanged.

All energy thresholds are in kcal/mol.  All stored conformer energies are
normalised to kcal/mol by OptimisationStage.

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

class PruningStage(BaseStage):
    """
    Conformer pruning stage.

    Applies one or more pruning methods sequentially to each molecule's
    conformer pool.  Serial execution.

    Stages never call mark_complete() — BaseStage.run() owns lifecycle.
    """

    # Ordered pruning dispatch table
    PRUNING_METHODS = {
        "select_conf":           "_prune_select_conf",
        "remove_conf":           "_prune_remove_conf",
        "rdkit_post_opt_legacy": "_prune_rdkit_post_opt_legacy",
        "rmsd_threshold":        "_prune_rmsd",
        "energy_window":         "_prune_energy_window",
        "max_energy":            "_prune_max_energy",
        "percentile":            "_prune_percentile",
        "n":                     "_prune_keep_lowest_n",
        "n_high":                "_prune_keep_highest_n",
        "conf_select":           "_prune_conf_select",
        "window_pair":           "_prune_window_pair",
        "lowest_outside_window": "_prune_lowest_outside_window",
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
        self._current_inchi_key = inchi_key

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

        def _b(k):
            v = params.get(k)
            if v is None:
                return None
            if isinstance(v, bool):
                return v
            if isinstance(v, str):
                return v.lower() not in ("false", "0", "no", "")
            return bool(v)

        def _list_str(k):
            v = params.get(k)
            if v is None:
                return None
            if isinstance(v, str):
                return [v]
            if isinstance(v, list):
                return [str(item) for item in v]
            self.log_warning(f"Expected string or list for '{k}': {v!r}")
            return None

        def _list_int(k):
            v = params.get(k)
            if v is None:
                return None
            if not isinstance(v, list):
                self.log_warning(f"Expected list for '{k}': {v!r}")
                return None
            result = []
            for item in v:
                try:
                    result.append(int(item))
                except Exception:
                    self.log_warning(f"Could not parse int in '{k}': {item!r}")
            return result if result else None

        return {
            "select_conf":           _list_str("select_conf"),
            "remove_conf":           _list_str("remove_conf"),
            "rdkit_post_opt_legacy": _b("rdkit_post_opt_legacy"),
            "rmsd_threshold":        _f("rmsd_threshold"),
            "energy_window":         _f("energy_window"),
            "max_energy":            _f("max_energy"),
            "percentile":            _f("percentile"),
            "n":                     _i("n"),
            "n_high":                _i("n_high"),
            "conf_select":           _list_int("conf_select"),
            "window_pair":           _f("window_pair"),
            "lowest_outside_window": _f("lowest_outside_window"),
        }

    # =========================================================================
    # Pruning methods  (all energies in kcal/mol)
    # =========================================================================

    def _prune_select_conf(self, conformers: list, inchi_keys: list) -> list:
        """Keep all conformers if current molecule is in inchi_keys; drop all otherwise."""
        if self._current_inchi_key in inchi_keys:
            self.log_info(
                f"select_conf: {self._current_inchi_key} matched — "
                f"keeping all {len(conformers)} conformers"
            )
            return conformers
        self.log_info(
            f"select_conf: {self._current_inchi_key} not in list — "
            f"dropping all {len(conformers)} conformers"
        )
        return []

    def _prune_remove_conf(self, conformers: list, inchi_keys: list) -> list:
        """Drop all conformers if current molecule is in inchi_keys; keep all otherwise."""
        if self._current_inchi_key in inchi_keys:
            self.log_info(
                f"remove_conf: {self._current_inchi_key} matched — "
                f"dropping all {len(conformers)} conformers"
            )
            return []
        self.log_info(
            f"remove_conf: {self._current_inchi_key} not in list — "
            f"keeping all {len(conformers)} conformers"
        )
        return conformers

    def _prune_rdkit_post_opt_legacy(
        self, conformers: list, enabled: bool
    ) -> list:
        """
        Replicates the 4 post-MMFF filters applied inside the IC colleague
        pipeline's calculate_rdkit() (ConformerGenerator_IC.py lines 1097–1131).

        Applied in order:
          1. Convergence  — remove if latest optimisation_history status != "converged"
          2. Overlap      — remove if any atom pair < 0.1 Å (IC threshold)
          3. Boltzmann    — remove if rel_prob < 1% at T=298.15 K (~2.73 kcal/mol)
          4. RMSD 1.0 Å   — cluster, keep lower-energy representative (IC rms_threshold)
        """
        if not enabled:
            return conformers

        import math
        import numpy as np

        # ── 1. Convergence filter ─────────────────────────────────────────────
        survivors = []
        for rec in conformers:
            hist   = rec.optimisation_history or []
            status = hist[-1].get("status", "unknown") if hist else "unknown"
            if status == "converged":
                survivors.append(rec)
            else:
                self.log_info(
                    f"rdkit_post_opt_legacy: {rec.lookup_id} removed — "
                    f"MMFF status={status!r}"
                )
        self.log_info(
            f"rdkit_post_opt_legacy convergence: {len(conformers)} → {len(survivors)}"
        )
        if not survivors:
            return survivors

        # ── 2. Overlapping atoms filter ───────────────────────────────────────
        OVERLAP_THRESHOLD = 0.1  # Å — matches IC pipeline

        def _has_overlapping_atoms(xyz_path: str) -> bool:
            mol = Chem.MolFromXYZFile(xyz_path)
            if mol is None:
                return False
            pos  = mol.GetConformer().GetPositions()     # (N, 3)
            diff = pos[:, None, :] - pos[None, :, :]    # (N, N, 3)
            dist = np.linalg.norm(diff, axis=2)         # (N, N)
            np.fill_diagonal(dist, 1.0)                 # ignore self
            return bool(np.any(dist < OVERLAP_THRESHOLD))

        before    = len(survivors)
        survivors = [r for r in survivors if not _has_overlapping_atoms(r.xyz_path)]
        self.log_info(
            f"rdkit_post_opt_legacy overlap (<{OVERLAP_THRESHOLD} Å): "
            f"{before} → {len(survivors)}"
        )
        if not survivors:
            return survivors

        # ── 3. Boltzmann filter ───────────────────────────────────────────────
        # IC: rel_prob = exp(-E/RT) / exp(-E_min/RT); remove if rel_prob < 1e-2
        # Energies in kcal/mol.  RT = R·T / 4184 ≈ 0.5927 kcal/mol at 298.15 K.
        RT_KCAL = 8.314 * 298.15 / 4184.0
        e_min   = min(r.energy for r in survivors)
        before  = len(survivors)
        survivors = [
            r for r in survivors
            if math.exp(-(r.energy - e_min) / RT_KCAL) >= 0.01
        ]
        self.log_info(
            f"rdkit_post_opt_legacy Boltzmann (T=298.15 K, cutoff=1e-2): "
            f"{before} → {len(survivors)}"
        )
        if not survivors:
            return survivors

        # ── 4. RMSD 1.0 Å filter ─────────────────────────────────────────────
        # IC: _get_idx_to_keep_by_rms_window(rms_threshold=1.0)
        before    = len(survivors)
        survivors = self._prune_rmsd(survivors, rmsd_threshold=1.0)
        self.log_info(
            f"rdkit_post_opt_legacy RMSD (1.0 Å): {before} → {len(survivors)}"
        )
        return survivors

    @staticmethod
    def _heavy_atom_mol(mol: Chem.Mol) -> Chem.Mol:
        """
        Return a new mol containing only non-hydrogen atoms with their 3D
        positions.  Used to compute RMSD over heavy atoms only, matching the
        legacy IC pipeline's rms_only_heavy_atoms=True behaviour (lines 248-250
        of ConformerGenerator_IC.py).

        Works on mols loaded from XYZ (no bond information required).
        """
        rw    = Chem.RWMol()
        conf  = mol.GetConformer()
        heavy = [(i, atom) for i, atom in enumerate(mol.GetAtoms())
                 if atom.GetAtomicNum() != 1]
        new_conf = Chem.Conformer(len(heavy))
        for j, (i, atom) in enumerate(heavy):
            rw.AddAtom(Chem.Atom(atom.GetAtomicNum()))
            new_conf.SetAtomPosition(j, conf.GetAtomPosition(i))
        rw.AddConformer(new_conf, assignId=True)
        return rw.GetMol()

    def _prune_rmsd(self, conformers: list, rmsd_threshold: float) -> list:
        if not conformers or rmsd_threshold <= 0:
            return conformers

        sorted_conf   = sorted(conformers, key=lambda c: c.energy)
        survivors     = []
        survivor_mols = []

        for rec in sorted_conf:
            try:
                mol = Chem.MolFromXYZFile(rec.xyz_path)
                if mol is None:
                    self.log_warning(
                        f"RMSD pruning: could not load XYZ for {rec.lookup_id}"
                    )
                    continue
                mol = self._heavy_atom_mol(mol)
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
            f"RMSD pruning (heavy atoms, threshold={rmsd_threshold} Å): "
            f"{len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_energy_window(self, conformers: list, energy_window: float) -> list:
        if not conformers:
            return conformers

        energies  = np.array([c.energy for c in conformers], dtype=float)
        min_e     = float(energies.min())
        survivors = [
            c for c, e in zip(conformers, energies)
            if (e - min_e) <= energy_window
        ]
        self.log_info(
            f"Energy window (ΔE ≤ {energy_window} kcal/mol): "
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

    def _prune_conf_select(self, conformers: list, indices: list) -> list:
        survivors = []
        for idx in indices:
            if idx < 0 or idx >= len(conformers):
                self.log_warning(
                    f"conf_select: index {idx} out of range "
                    f"(0–{len(conformers) - 1}) — skipped"
                )
            else:
                survivors.append(conformers[idx])
        self.log_info(
            f"conf_select {indices}: {len(conformers)} → {len(survivors)}"
        )
        return survivors

    def _prune_window_pair(self, conformers: list, energy_window: float) -> list:
        """
        Testing/paper method — keeps exactly two conformers separated by an
        energy window:
          1. The global energy minimum (lowest-energy conformer).
          2. The lowest-energy conformer whose energy exceeds
             min_energy + energy_window (kcal/mol).

        If no conformer exists outside the window, returns just the minimum
        with a warning.  Not intended for normal pipeline runs.
        """
        if not conformers:
            return conformers

        sorted_conf = sorted(conformers, key=lambda c: c.energy)
        minimum     = sorted_conf[0]
        e_min       = minimum.energy

        outside = [c for c in sorted_conf if (c.energy - e_min) > energy_window]

        if outside:
            pair = [minimum, outside[0]]
            self.log_info(
                f"window_pair (window={energy_window} kcal/mol): "
                f"{len(conformers)} → 2 "
                f"(min={e_min:.3f}, outlier={outside[0].energy:.3f} kcal/mol)"
            )
            return pair

        self.log_warning(
            f"window_pair: no conformer found above window={energy_window} kcal/mol "
            f"(min={e_min:.3f}, max={sorted_conf[-1].energy:.3f}) — returning minimum only"
        )
        return [minimum]

    def _prune_lowest_outside_window(
        self, conformers: list, energy_window: float
    ) -> list:
        """
        Testing/paper method — keeps only the lowest-energy conformer whose
        energy lies strictly above min_energy + energy_window (kcal/mol).

        Pair with a separate n=1 job to compare best vs worst-of-window
        conformer solubilities independently.  Returns an empty list (with a
        warning) if no conformer exists outside the window.
        """
        if not conformers:
            return conformers

        sorted_conf = sorted(conformers, key=lambda c: c.energy)
        e_min       = sorted_conf[0].energy

        outside = [c for c in sorted_conf if (c.energy - e_min) > energy_window]

        if outside:
            self.log_info(
                f"lowest_outside_window (window={energy_window} kcal/mol): "
                f"{len(conformers)} → 1 "
                f"(selected={outside[0].energy:.3f}, min={e_min:.3f} kcal/mol)"
            )
            return [outside[0]]

        self.log_warning(
            f"lowest_outside_window: no conformer found above "
            f"window={energy_window} kcal/mol "
            f"(min={e_min:.3f}, max={sorted_conf[-1].energy:.3f}) — returning empty"
        )
        return []

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
            "select_conf":             "Select Molecules (inchi keys)",
            "remove_conf":             "Remove Molecules (inchi keys)",
            "total_conformers":        "Total Conformers",
            "valid_energy_conformers": "With Valid Energy",
            "removed_missing_energy":  "Removed (Missing Energy)",
            "kept_after_pruning":      "Final Count",
            "methods_used":            "Pruning Steps Applied",
            "rmsd_threshold":          "RMSD Threshold (Å)",
            "energy_window":           "Energy Window (kcal/mol)",
            "max_energy":              "Max Energy (kcal/mol)",
            "percentile":              "Percentile Cutoff",
            "n":                       "Keep Lowest N",
            "n_high":                  "Keep Highest N",
            "conf_select":             "Conformer Select (indices)",
            "window_pair":             "Window Pair (kcal/mol)",
            "lowest_outside_window":   "Lowest Outside Window (kcal/mol)",
        })

        col_order = [
            "Molecule", "Total Conformers", "With Valid Energy",
            "Removed (Missing Energy)", "Final Count",
            "Pruning Steps Applied",
            "Select Molecules (inchi keys)", "Remove Molecules (inchi keys)",
            "RMSD Threshold (Å)",
            "Energy Window (kcal/mol)", "Max Energy (kcal/mol)",
            "Percentile Cutoff", "Keep Lowest N", "Keep Highest N",
            "Conformer Select (indices)",
            "Window Pair (kcal/mol)", "Lowest Outside Window (kcal/mol)",
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