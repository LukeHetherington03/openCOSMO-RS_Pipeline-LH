#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modules/stages/solubility_stage.py

COSMO-RS solubility stage with solvent list support.

=========================================================================
INPUT
=========================================================================

  stage_input : orcacosmo_summary.json from OrcacosmoStage
    List of conformer entries, one per conformer:
      lookup_id, inchi_key, conf_num, xyz_path, energy, smiles,
      provenance, optimisation_history, orcacosmo_history, orcacosmo_path

    Multiple conformers per molecule are grouped by inchi_key.

=========================================================================
OUTPUT  (canonical stage output)
=========================================================================

  solubility_results.json
    List of per-molecule result dicts — see SCHEMA section below.

=========================================================================
SCHEMA  (solubility_results.json entry)
=========================================================================

  {
    "inchi_key":                           str,
    "mol_name":                            str,
    "smiles":                              str,
    "melting_temp":                        float | "liquid" | "N/A",
    "melting_temp_source":                 str,
    "experimental_solubility_mol_frac":    float | null,
    "aqsol_predicted_solubility_mol_frac": float | null,
    "conformer_history": [
      {
        "lookup_id":            str,
        "conf_num":             int | null,
        "energy":               float | null,
        "provenance":           dict,
        "optimisation_history": list,
        "orcacosmo_history":    dict,
        "orcacosmo_path":       str
      }
    ],
    "solvent_results": [
      {
        "combo_id":             str,
        "combo_label":          str,
        "solvents":             [{"name": str, "mol_frac": float}],
        "temperature_K":        float,
        "predicted_solubility": float | null,
        "n_solute_confs":       int,
        "n_solvent_confs":      {str: int},
        "mixture_inputs_path":  str,
        "raw_output_path":      str,
        "error":                str | null
      }
    ]
  }

=========================================================================
AUXILIARY OUTPUTS  (per molecule)
=========================================================================

  results/<inchi_key>/
      inputs/
          solute_raw/           original .orcacosmo (provenance, unmodified)
          solute/               renamed <inchi_key>_cNNN.orcacosmo
          solvents/<name>/      conformers per solvent (shared across combos)
      <combo_id>/
          mixture_inputs.txt
          raw_output.txt
          solubility.json       combo-level result (supports mid-molecule resume)
  results/<inchi_key>.json      per-molecule assembled result

  solubility_results.csv        flat: one row per (molecule, combo)
  solubility_human_summary.txt  aligned fixed-width table
  checkpoints/<inchi_key>.json  written atomically after all combos complete

=========================================================================
STAGE ARGUMENTS  (pipeline_spec args)
=========================================================================

  solvent_list : str, optional
    Filename of a solvent list JSON in CONSTANT_FILES/solvent_lists/.
    Can be a stem ("water_only"), filename ("water_only.json"), or
    absolute path.  If omitted, falls back to default_list.json, then
    to a built-in water-only emergency combo.

  temperature : float, optional
    Stage-level temperature override (Kelvin).
    Hierarchy: pipeline_spec > combo "temperature" field > defaults JSON.

  verbose : bool, optional  [default: false]

  strict : bool, optional  [default: false]
    If true, any missing solvent directory aborts the stage immediately.
    If false, combos using missing solvents are skipped with a warning;
    the stage fails only if zero valid combos remain.

=========================================================================
SOLVENT LIST FILE FORMAT
=========================================================================

  Location: CONSTANT_FILES/solvent_lists/<filename>

  {
    "description": "...",
    "combos": {
      "water_only": {
        "label": "Pure water",
        "solvents": [{"name": "water", "mol_frac": 1.0}]
      },
      "60_40_water_ethanol": {
        "label": "60:40 water/ethanol",
        "temperature": 310.15,
        "solvents": [
          {"name": "water",   "mol_frac": 0.6},
          {"name": "ethanol", "mol_frac": 0.4}
        ]
      }
    }
  }

  Mol fracs are treated as ratios and normalised automatically.
  Solvent conformer dirs resolved as: CONSTANT_FILES/solvents/<name>/

=========================================================================
SOLVENT LIST FALLBACK CHAIN
=========================================================================

  1. parameters["solvent_list"]           user-specified
  2. solvent_lists_dir/default_list.json  stage default
  3. Built-in water-only emergency combo  always works

=========================================================================
PARALLELISM
=========================================================================

  Uses BaseStage.run_parallel().
  Item granularity: one molecule (inchi_key) per worker.
  All valid combos for a molecule run sequentially inside the worker.

  Solubility calculations are single-core (COSMO-RS engine).
  cores_per_item will typically be 1 for this stage.

=========================================================================
COMBO-LEVEL RESUME
=========================================================================

  If a molecule's worker crashes mid-combo (e.g. OOM), on the next run:
  - BaseStage skips molecules whose checkpoint already exists.
  - For molecules without a checkpoint, the worker checks whether
    results/<inchi_key>/<combo_id>/solubility.json already exists and
    is valid before running that combo.  Completed combos are skipped.

=========================================================================
METADATA FIELDS READ  (written by CleaningStage v3+)
=========================================================================

  inchi_key, smiles, mol_name, mol_name_iupac
  melting_temp, melting_temp_source
  experimental_solubility_mol_frac
  aqsol_predicted_solubility_mol_frac
  Gfus_mode, Hfus

  Legacy fallbacks (pre-v3):
    key_inchi      -> inchi_key
    Tm             -> melting_temp
    melting_point  -> melting_temp

=========================================================================
STRICT MODE HIERARCHY
=========================================================================

  1. parameters["strict"]               runtime override
  2. solubility_defaults.json "strict"  stage default
  3. False                              fallback

"""

import csv
import json
import math
import multiprocessing
import os
import shutil
import traceback
from datetime import datetime
from pathlib import Path

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter
from modules.utils.mixture_inputs_builder import (
    build_mixture_inputs,
    normalise_solvent_ratios,
)
from modules.solubility_engine.legacy_cpp_wrapper import run_legacy_cosmors


_WATER_MOLARITY = 55.51  # mol/L at 25 °C


def _mol_frac_to_logS(x):
    try:
        v = float(x)
        if v <= 0.0 or v >= 1.0:
            return None
        return math.log10(v * _WATER_MOLARITY / (1.0 - v))
    except (TypeError, ValueError):
        return None


# =============================================================================
# Module-level worker — must be at module scope to be picklable
# =============================================================================

def _solubility_worker(args: dict) -> dict:
    """
    Pure worker function.  Processes a single molecule (inchi_key):
      - stages solute and solvent .orcacosmo files
      - for each valid combo: builds mixture_inputs.txt, runs COSMO-RS
      - writes per-molecule outputs and returns checkpoint dict

    Combo-level resume: if results/<inchi_key>/<combo_id>/solubility.json
    already exists and is valid, that combo is skipped.
    """
    worker_name = multiprocessing.current_process().name
    worker_pid  = os.getpid()

    item_id      = args["item_id"]
    inchi_key    = args["inchi_key"]
    conformers   = args["conformers"]        # list of conformer entry dicts
    outputs_dir  = args["outputs_dir"]
    valid_combos = args["valid_combos"]      # {combo_id: cdef} already validated
    solvent_bank = args["solvent_bank"]      # path to solvents/ root dir
    metadata_dir = args["metadata_dir"]
    defaults     = args["defaults"]
    opencosmo    = args["opencosmo"]
    verbose      = args["verbose"]

    results_dir = os.path.join(outputs_dir, "results")

    try:
        import time
        t_start = time.perf_counter()

        # ── Per-molecule directory structure ──────────────────────────────────
        mol_dir        = os.path.join(results_dir, inchi_key)
        inputs_dir     = os.path.join(mol_dir, "inputs")
        solute_raw_dir = os.path.join(inputs_dir, "solute_raw")
        solute_dir     = os.path.join(inputs_dir, "solute")
        solvents_dir   = os.path.join(inputs_dir, "solvents")

        for d in (mol_dir, inputs_dir, solute_raw_dir, solute_dir, solvents_dir):
            os.makedirs(d, exist_ok=True)

        # ── Stage solute conformers ───────────────────────────────────────────
        all_orca_paths = [Path(c["orcacosmo_path"]) for c in conformers]

        for p in all_orca_paths:
            shutil.copy(str(p), os.path.join(solute_raw_dir, p.name))

        for i, p in enumerate(sorted(all_orca_paths)):
            shutil.copy(str(p), os.path.join(solute_dir, f"{inchi_key}_c{i:03d}.orcacosmo"))

        n_solute = len(list(Path(solute_dir).glob("*.orcacosmo")))

        # ── Stage all required solvents once (shared across combos) ──────────
        all_solvent_names = {
            s["name"]
            for cdef in valid_combos.values()
            for s in cdef["solvents"]
        }

        copied_solvent_dirs = {}   # name -> str path inside results tree
        n_solvent_confs_map = {}   # name -> int

        for solvent_name in all_solvent_names:
            src_dir  = os.path.join(solvent_bank, solvent_name)
            dest_dir = os.path.join(solvents_dir, solvent_name)
            os.makedirs(dest_dir, exist_ok=True)

            files = sorted(Path(src_dir).glob("*.orcacosmo"))
            for p in files:
                shutil.copy(str(p), os.path.join(dest_dir, p.name))

            copied_solvent_dirs[solvent_name] = dest_dir
            n_solvent_confs_map[solvent_name] = len(files)

        # ── Load metadata ─────────────────────────────────────────────────────
        meta = _load_metadata(inchi_key, metadata_dir)

        # ── Run each valid combo ──────────────────────────────────────────────
        solvent_results = []
        for combo_id, cdef in valid_combos.items():
            result = _run_combo(
                inchi_key           = inchi_key,
                mol_dir             = mol_dir,
                combo_id            = combo_id,
                cdef                = cdef,
                meta                = meta,
                solute_dir          = solute_dir,
                n_solute            = n_solute,
                copied_solvent_dirs = copied_solvent_dirs,
                n_solvent_confs_map = n_solvent_confs_map,
                defaults            = defaults,
                opencosmo           = opencosmo,
                verbose             = verbose,
            )
            solvent_results.append(result)

        elapsed = time.perf_counter() - t_start

        # ── Build conformer_history ───────────────────────────────────────────
        conformer_history = [
            {
                "lookup_id":            c["lookup_id"],
                "conf_num":             c.get("conf_num"),
                "energy":               c.get("energy"),
                "provenance":           c.get("provenance", {}),
                "optimisation_history": c.get("optimisation_history", []),
                "orcacosmo_history":    c.get("orcacosmo_history", {}),
                "orcacosmo_path":       c.get("orcacosmo_path"),
            }
            for c in conformers
        ]

        # ── Assemble per-molecule result ──────────────────────────────────────
        mol_result = {
            "inchi_key":                           inchi_key,
            "mol_name":                            meta["mol_name"],
            "smiles":                              meta["smiles"],
            "melting_temp":                        meta["Tm"],
            "melting_temp_source":                 meta["Tm_source"],
            "experimental_solubility_mol_frac":    meta["experimental_solubility_mol_frac"],
            "log_solubility_experimental":         _mol_frac_to_logS(meta["experimental_solubility_mol_frac"]),
            "aqsol_predicted_solubility_mol_frac": meta["aqsol_predicted_solubility_mol_frac"],
            "conformer_history":                   conformer_history,
            "solvent_results":                     solvent_results,
        }

        mol_json_path = os.path.join(results_dir, f"{inchi_key}.json")
        with AtomicWriter(mol_json_path) as f:
            json.dump(mol_result, f, indent=2)

        # ── Determine overall status ──────────────────────────────────────────
        n_ok   = sum(1 for r in solvent_results if r["predicted_solubility"] is not None)
        n_fail = len(solvent_results) - n_ok

        if n_ok == 0:
            raise RuntimeError(
                f"All {len(solvent_results)} combo(s) failed for {inchi_key}"
            )

        return {
            # BaseStage contract
            "item_id":     item_id,
            "status":      "ok",
            "worker_name": worker_name,
            "worker_pid":  worker_pid,
            "log_details": [
                ("combos",   f"{n_ok} ok / {n_fail} failed"),
                ("confs",    f"{n_solute} solute"),
                ("elapsed",  f"{elapsed:.1f}s"),
            ],
            # Stage payload — full mol_result inlined for checkpoint
            **mol_result,
        }

    except Exception as e:
        tb = traceback.format_exc()
        return {
            "item_id":     item_id,
            "status":      "failed",
            "error":       str(e),
            "traceback":   tb,
            "worker_name": worker_name,
            "worker_pid":  worker_pid,
            "log_details": [],
        }


# =============================================================================
# Combo runner
# =============================================================================

def _run_combo(
    inchi_key:           str,
    mol_dir:             str,
    combo_id:            str,
    cdef:                dict,
    meta:                dict,
    solute_dir:          str,
    n_solute:            int,
    copied_solvent_dirs: dict,
    n_solvent_confs_map: dict,
    defaults:            dict,
    opencosmo:           dict,
    verbose:             bool,
) -> dict:
    """Run a single solvent combo for one molecule."""
    combo_dir = os.path.join(mol_dir, combo_id)
    os.makedirs(combo_dir, exist_ok=True)

    mix_path = os.path.join(combo_dir, "mixture_inputs.txt")
    raw_path = os.path.join(combo_dir, "raw_output.txt")
    sol_path = os.path.join(combo_dir, "solubility.json")
    temperature = cdef["temperature"]

    # ── Combo-level resume ────────────────────────────────────────────────────
    if os.path.exists(sol_path):
        try:
            with open(sol_path) as f:
                existing = json.load(f)
            if existing.get("predicted_solubility") is not None or existing.get("error"):
                return existing   # already done
        except Exception:
            pass  # corrupt — rerun

    # ── Normalise solvent mol fracs ───────────────────────────────────────────
    try:
        solvents_normalised = normalise_solvent_ratios(
            cdef["solvents"],
            defaults.get("initial_solute_x", 0),
        )
    except Exception as e:
        return _failed_combo(
            combo_id, cdef, temperature, n_solute,
            n_solvent_confs_map, sol_path, raw_path, mix_path,
            f"Solvent normalisation failed: {e}",
        )

    # ── Build mixture_inputs.txt and run COSMO-RS (with SORcf fallback) ───────
    combo_solvent_dirs = {
        s["name"]: Path(copied_solvent_dirs[s["name"]])
        for s in cdef["solvents"]
    }
    combo_n_solvent_confs = {
        s["name"]: n_solvent_confs_map[s["name"]]
        for s in cdef["solvents"]
    }

    sorcf_schedule = defaults["saturation"].get(
        "sorcf_fallback", [defaults["saturation"].get("SORcf", 1.0)]
    )
    timeout = defaults.get("cosmors_timeout", 1500)

    result     = None
    used_sorcf = sorcf_schedule[0]
    last_error = None

    for sorcf in sorcf_schedule:
        attempt_defaults = {
            **defaults,
            "saturation": {**defaults["saturation"], "SORcf": sorcf},
        }

        try:
            mixture_text = build_mixture_inputs(
                solute_meta = {
                    "inchi_key":    inchi_key,
                    "smiles":       meta["smiles"],
                    "melting_temp": meta["Tm"],
                    "Gfus_mode":    meta["Gfus_mode"],
                    "Hfus":         meta["Hfus"],
                },
                solute_dir           = Path(solute_dir),
                solvent_dirs         = combo_solvent_dirs,
                n_solute_confs       = n_solute,
                n_solvent_confs      = combo_n_solvent_confs,
                defaults             = attempt_defaults,
                temperature          = temperature,
                solvents             = solvents_normalised,
                solute_multiplicity  = int(meta.get("multiplicity", 1)),
            )
        except Exception as e:
            return _failed_combo(
                combo_id, cdef, temperature, n_solute,
                n_solvent_confs_map, sol_path, raw_path, mix_path,
                f"mixture_inputs build failed: {e}",
            )

        with open(mix_path, "w") as f:
            f.write(mixture_text)

        try:
            result = run_legacy_cosmors(
                mixture_text,
                python_src    = opencosmo["python_src"],
                cpp_bindings  = opencosmo["cpp_bindings"],
                driver_script = opencosmo["python_driver"],
                timeout       = timeout,
            )
        except Exception as e:
            with open(raw_path, "w") as f:
                f.write("")
            last_error = f"COSMO-RS engine error: {e}"
            result = None
            break  # non-timeout error — don't retry with different SORcf

        if result.get("timed_out"):
            last_error = f"timed out after {timeout}s with SORcf={sorcf}"
            result = None
            continue  # try next SORcf value

        # Engine responded (converged or returned no solubility) — stop retrying
        used_sorcf = sorcf
        break

    with open(raw_path, "w") as f:
        f.write((result or {}).get("raw_stdout") or "")

    predicted = (result or {}).get("solubility")
    if predicted is None:
        return _failed_combo(
            combo_id, cdef, temperature, n_solute,
            n_solvent_confs_map, sol_path, raw_path, mix_path,
            last_error or "COSMO-RS engine returned no solubility value",
        )

    combo_result = {
        "combo_id":                   combo_id,
        "combo_label":                cdef["label"],
        "solvents":                   solvents_normalised,
        "temperature_K":              temperature,
        "predicted_solubility":       predicted,
        "log_solubility_predicted":   _mol_frac_to_logS(predicted),
        "saturation_mole_fractions":  result.get("saturation_mole_fractions", []),
        "n_solute_confs":             n_solute,
        "n_solvent_confs":            combo_n_solvent_confs,
        "mixture_inputs_path":        mix_path,
        "raw_output_path":            raw_path,
        "sorcf_used":                 used_sorcf,
        "error":                      None,
    }

    with AtomicWriter(sol_path) as f:
        json.dump(combo_result, f, indent=2)

    return combo_result


def _failed_combo(
    combo_id:            str,
    cdef:                dict,
    temperature:         float,
    n_solute:            int,
    n_solvent_confs_map: dict,
    sol_path:            str,
    raw_path:            str,
    mix_path:            str,
    error_msg:           str,
) -> dict:
    combo_n_solvent_confs = {
        s["name"]: n_solvent_confs_map.get(s["name"], 0)
        for s in cdef["solvents"]
    }
    result = {
        "combo_id":             combo_id,
        "combo_label":          cdef.get("label", combo_id),
        "solvents":             cdef["solvents"],
        "temperature_K":        temperature,
        "predicted_solubility": None,
        "n_solute_confs":       n_solute,
        "n_solvent_confs":      combo_n_solvent_confs,
        "mixture_inputs_path":  mix_path,
        "raw_output_path":      raw_path,
        "error":                error_msg,
    }
    with AtomicWriter(sol_path) as f:
        json.dump(result, f, indent=2)
    return result


# =============================================================================
# Metadata loader  (module-level — called from worker)
# =============================================================================

def _load_metadata(inchi_key: str, metadata_dir: str) -> dict:
    """
    Load per-molecule metadata JSON.

    Melting temp resolution order:
      1. meta["melting_temp"]    canonical (CleaningStage v3+)
      2. meta["Tm"]              legacy field
      3. meta["melting_point"]   legacy field
      4. PubChem live query      written back to disk on success
    Raises RuntimeError if no melting point can be resolved.
    """
    meta_path = os.path.join(metadata_dir, f"{inchi_key}.json")
    if not os.path.exists(meta_path):
        raise RuntimeError(
            f"Metadata file not found for {inchi_key}: {meta_path}"
        )
    with open(meta_path) as f:
        meta = json.load(f)

    def _try_float(val):
        if val is None or val == "" or val == "N/A":
            return None
        try:
            v = float(val)
            return None if v != v else v   # guard NaN
        except (TypeError, ValueError):
            return None

    # Legacy inchi_key fallback
    ik = meta.get("inchi_key") or meta.get("key_inchi", inchi_key)

    # ── Melting temperature resolution ───────────────────────────────────────
    Tm     = None
    Tm_src = None

    # "liquid" string is a valid special case (no melting point needed)
    raw_tm = meta.get("melting_temp")
    if isinstance(raw_tm, str) and raw_tm.strip().lower() == "liquid":
        Tm     = "liquid"
        Tm_src = meta.get("melting_temp_source", "provided")

    if Tm is None:
        Tm = _try_float(meta.get("melting_temp"))
        if Tm is not None:
            Tm_src = meta.get("melting_temp_source", "metadata")

    if Tm is None:
        Tm = _try_float(meta.get("Tm"))
        if Tm is not None:
            Tm_src = "legacy_Tm_field"

    if Tm is None:
        Tm = _try_float(meta.get("melting_point"))
        if Tm is not None:
            Tm_src = "legacy_melting_point_field"

    if Tm is None:
        # PubChem live query as last resort
        try:
            from modules.utils.molecule_utils import MoleculeUtils
            mp_info = MoleculeUtils.get_melting_point(inchikey=inchi_key)
            Tm      = _try_float(mp_info.get("melting_temp"))
            if Tm is not None:
                Tm_src = "pubchem_runtime"
                _write_back_melting_temp(
                    meta_path, meta, Tm, Tm_src,
                    mp_info.get("melting_temp_detail", ""),
                )
        except Exception:
            pass

    if Tm is None:
        raise ValueError(
            f"No melting point available for {inchi_key}. "
            "Provide melting_temp in the input CSV or metadata JSON."
        )

    mol_name = (
        meta.get("mol_name")
        or meta.get("mol_name_iupac")
        or ik
    )

    return {
        "inchi_key":                           ik,
        "smiles":                              meta.get("smiles", ""),
        "mol_name":                            mol_name,
        "Tm":                                  Tm,
        "Tm_source":                           Tm_src,
        "Gfus_mode":                           meta.get("Gfus_mode", "MyrdalYalkowsky"),
        "Hfus":                                meta.get("Hfus", "N/A"),
        "experimental_solubility_mol_frac":    meta.get("experimental_solubility_mol_frac"),
        "aqsol_predicted_solubility_mol_frac": meta.get("aqsol_predicted_solubility_mol_frac"),
    }


def _write_back_melting_temp(
    meta_path: str,
    meta:      dict,
    Tm:        float,
    Tm_src:    str,
    detail:    str,
):
    """Write PubChem-resolved melting temp back to metadata JSON."""
    try:
        updated = {
            **meta,
            "melting_temp":        Tm,
            "melting_temp_c":      round(Tm - 273.15, 4),
            "melting_temp_source": Tm_src,
            "melting_temp_detail": detail,
        }
        with AtomicWriter(meta_path) as f:
            json.dump(updated, f, indent=2)
    except Exception:
        pass   # non-fatal — just couldn't cache the result


# =============================================================================
# Stage class
# =============================================================================

class SolubilityStage(BaseStage):
    """
    COSMO-RS solubility stage with solvent list support.

    Item granularity: one molecule (inchi_key) per worker.
    All valid combos for a molecule run sequentially inside the worker.

    execute() loads config, resolves and validates the solvent list,
    groups conformers by molecule, then delegates to run_parallel().
    Final outputs are assembled from checkpoints.

    Stages never call mark_complete() — BaseStage.run() owns lifecycle.
    """

    # =========================================================================
    # Entry point
    # =========================================================================

    def execute(self):
        self.set_stage_output("solubility_results.json")

        self._prepare_directories()
        self._load_stage_config()
        self._load_solvent_list()
        self._validate_solvent_dirs()
        self._load_orcacosmo_summary()

        inchi_keys = [e["inchi_key"] for e in self._molecules]

        if self.job.pending_items:
            self.log_info(
                f"Resuming — {len(self.job.pending_items)} pending / "
                f"{len(inchi_keys) - len(self.job.pending_items)} "
                f"already complete"
            )
        else:
            self.set_items(inchi_keys)
            self._clear_checkpoints()
            self.log_info(
                f"Fresh run — {len(inchi_keys)} molecules, "
                f"{len(self.valid_combos)} combo(s) each"
            )

        self.run_parallel(_solubility_worker, self._build_args)
        self._assemble_output()

    # =========================================================================
    # Directories
    # =========================================================================

    def _prepare_directories(self):
        self.results_dir = os.path.join(self.outputs_dir, "results")
        os.makedirs(self.results_dir, exist_ok=True)

    # =========================================================================
    # Config
    # =========================================================================

    def _load_stage_config(self):
        cfg = self.config
        if not cfg:
            self.fail("SolubilityStage: missing config")

        # ── Defaults file ─────────────────────────────────────────────────────
        defaults_path = cfg.get("solubility", {}).get("defaults")
        self._sol_defaults = {}
        if defaults_path and os.path.exists(defaults_path):
            with open(defaults_path) as f:
                self._sol_defaults = json.load(f)
        else:
            self.log_warning(
                "solubility_defaults.json not found — using built-in defaults"
            )

        if "n_workers" in self._sol_defaults:
            self.log_warning(
                "solubility_defaults.json contains 'n_workers' — "
                "this is owned by ResourceAllocator and will be ignored"
            )

        # ── Constant file dirs ────────────────────────────────────────────────
        const_cfg = cfg.get("constant_files", {})
        self.metadata_dir      = const_cfg.get("metadata_dir")
        self.solvent_bank      = Path(const_cfg.get("solvent_dir", ""))
        self.solvent_lists_dir = Path(const_cfg.get("solvent_lists", ""))

        if not self.metadata_dir:
            self.fail("Missing metadata_dir in config['constant_files']")
        if not self.solvent_bank or not self.solvent_bank.is_dir():
            self.fail(f"solvent_dir not found: {self.solvent_bank}")

        # ── opencosmo paths ───────────────────────────────────────────────────
        opencosmo_cfg = cfg.get("opencosmo", {})
        self.opencosmo = {
            "python_src":    opencosmo_cfg.get("python_src"),
            "cpp_bindings":  opencosmo_cfg.get("cpp_bindings"),
            "python_driver": opencosmo_cfg.get("python_driver"),
        }
        missing = [k for k, v in self.opencosmo.items() if not v]
        if missing:
            self.fail(f"Missing opencosmo config keys: {missing}")

        # ── Stage-level temperature (combo field overrides this) ──────────────
        self.stage_temperature = float(
            self.parameters.get("temperature")
            or self._sol_defaults.get("default_temperature")
            or 298.15
        )

        # ── Verbose ───────────────────────────────────────────────────────────
        self.verbose = bool(self.parameters.get("verbose", False))

        # ── Strict mode hierarchy ─────────────────────────────────────────────
        params_strict = self.parameters.get("strict")
        if params_strict is not None:
            self.strict_mode = bool(params_strict)
        else:
            self.strict_mode = bool(self._sol_defaults.get("strict", False))

        self.log_config(
            f"stage_temperature={self.stage_temperature} K  "
            f"strict={self.strict_mode}"
        )
        self.log_resources(
            f"n_workers={self.n_workers}  "
            f"cores_per_item={self.cores_per_item}"
        )

    # =========================================================================
    # Solvent list loading  (3-level fallback)
    # =========================================================================

    def _resolve_solvent_list_path(self, value: str) -> Path | None:
        """
        Locate a solvent list JSON from a user-supplied string.

        Resolution order:
          1. Absolute path or contains path separators — use directly.
          2. solvent_lists_dir / value           (exact filename)
          3. solvent_lists_dir / value + ".json" (stem only)
        Returns None if nothing matched.
        """
        p = Path(value)
        if p.is_absolute() or "/" in value or "\\" in value:
            return p if p.exists() else None

        candidate = self.solvent_lists_dir / value
        if candidate.exists():
            return candidate

        candidate_json = self.solvent_lists_dir / (value + ".json")
        if candidate_json.exists():
            return candidate_json

        return None

    def _load_solvent_list(self):
        """
        Load solvent combos via the 3-level fallback chain:

          1. parameters["solvent_list"]           user-specified
          2. solvent_lists_dir/default_list.json  stage default
          3. Built-in water-only emergency combo  always available

        Populates self.combos: {combo_id: cdef}
        """
        list_value  = self.parameters.get("solvent_list")
        list_path   = None
        source_note = ""

        # ── Level 1: user-specified ───────────────────────────────────────────
        if list_value:
            list_path = self._resolve_solvent_list_path(str(list_value))
            if list_path is None:
                self.log_warning(
                    f"solvent_list='{list_value}' could not be resolved "
                    f"(searched {self.solvent_lists_dir}) — "
                    "falling back to default_list.json"
                )
            else:
                source_note = f"user-specified ('{list_value}')"

        # ── Level 2: default_list.json ────────────────────────────────────────
        if list_path is None:
            default_path = self.solvent_lists_dir / "default_list.json"
            if default_path.exists():
                list_path   = default_path
                source_note = "default_list.json"
                if not list_value:
                    self.log_warning(
                        "No solvent_list specified — using default_list.json"
                    )

        # ── Level 3: built-in water-only emergency fallback ───────────────────
        if list_path is None:
            self.log_warning(
                "No solvent list available — "
                "using built-in water-only emergency fallback"
            )
            self._apply_combos(
                raw_combos  = {
                    "water_only": {
                        "label":    "Pure water (emergency fallback)",
                        "solvents": [{"name": "water", "mol_frac": 1.0}],
                    }
                },
                source_note = "built-in water-only emergency fallback",
            )
            return

        # ── Parse the resolved JSON file ──────────────────────────────────────
        try:
            with open(list_path) as f:
                raw = json.load(f)
        except Exception as e:
            self.fail(f"Failed to parse solvent list '{list_path}': {e}")

        raw_combos = raw.get("combos", {})
        if not raw_combos:
            self.fail(f"Solvent list '{list_path}' contains no combos")

        self._apply_combos(
            raw_combos  = raw_combos,
            description = raw.get("description", ""),
            source_note = source_note,
            list_path   = list_path,
        )

    def _apply_combos(
        self,
        raw_combos:  dict,
        description: str  = "",
        source_note: str  = "",
        list_path         = None,
    ):
        """
        Populate self.combos from a raw combos dict.

        Temperature hierarchy per combo:
          1. parameters["temperature"]   runtime override (applies to all combos)
          2. combo_def["temperature"]    combo-specific
          3. self.stage_temperature      from defaults / parameters
        """
        self.combos = {}
        for combo_id, combo_def in raw_combos.items():
            if "temperature" in self.parameters:
                temperature = float(self.parameters["temperature"])
            elif "temperature" in combo_def:
                temperature = float(combo_def["temperature"])
            else:
                temperature = self.stage_temperature

            self.combos[combo_id] = {
                "label":       combo_def.get("label", combo_id),
                "temperature": temperature,
                "solvents": [
                    {"name": s["name"], "ratio": float(s.get("mol_frac", 1.0))}
                    for s in combo_def["solvents"]
                ],
            }

        path_str = str(list_path) if list_path else "built-in"
        self.log_config(f"Loaded solvent list ({source_note}): {path_str}")
        if description:
            self.log_config(f"Description: {description}")
        self.log_config(
            f"{len(self.combos)} combo(s): {list(self.combos.keys())}"
        )

    # =========================================================================
    # Solvent dir validation
    # =========================================================================

    def _validate_solvent_dirs(self):
        """
        Check every solvent referenced in self.combos has a conformer dir.
        Populates self.valid_combos (subset of self.combos).
        Missing solvents: warn and skip combos unless strict_mode.
        """
        required = {
            s["name"]
            for cdef in self.combos.values()
            for s in cdef["solvents"]
        }

        valid_solvents   = set()
        invalid_solvents = set()

        for name in required:
            solvent_dir = self.solvent_bank / name
            files = (
                list(solvent_dir.glob("*.orcacosmo"))
                if solvent_dir.exists() else []
            )
            if files:
                valid_solvents.add(name)
            else:
                invalid_solvents.add(name)
                msg = (
                    f"Solvent '{name}' has no conformer directory at "
                    f"{solvent_dir}"
                )
                if self.strict_mode:
                    self.fail(msg)
                else:
                    self.log_warning(
                        f"{msg} — combos using this solvent will be skipped"
                    )

        self.valid_combos = {
            cid: cdef
            for cid, cdef in self.combos.items()
            if all(s["name"] in valid_solvents for s in cdef["solvents"])
        }

        skipped = set(self.combos) - set(self.valid_combos)
        if skipped:
            self.log_warning(
                f"Combos skipped (missing solvents): {sorted(skipped)}"
            )

        if not self.valid_combos:
            self.fail("No valid combos after solvent validation — cannot proceed")

        self.log_config(
            f"Valid combos: {list(self.valid_combos.keys())}"
        )

    # =========================================================================
    # Load and group orcacosmo_summary.json
    # =========================================================================

    def _load_orcacosmo_summary(self):
        summary_file = self.require_file(
            self.get_stage_input(), "orcacosmo_summary.json"
        )
        self.log_info(f"Stage input: {summary_file}")

        os.makedirs(self.inputs_dir, exist_ok=True)
        shutil.copy(
            summary_file,
            os.path.join(self.inputs_dir, "orcacosmo_summary.json")
        )

        with open(summary_file) as f:
            entries = json.load(f)

        # Group conformer entries by inchi_key
        grouped: dict = {}
        for entry in entries:
            ik = entry["inchi_key"]
            grouped.setdefault(ik, []).append(entry)

        self._molecules = [
            {"inchi_key": ik, "conformers": confs}
            for ik, confs in grouped.items()
        ]
        self._molecule_map = {m["inchi_key"]: m for m in self._molecules}

        total_confs = sum(len(m["conformers"]) for m in self._molecules)
        self.log_info(
            f"Loaded {len(self._molecules)} molecules / "
            f"{total_confs} conformers"
        )

    # =========================================================================
    # Build args for worker
    # =========================================================================

    def _build_args(self, item: str) -> dict:
        """item = inchi_key"""
        mol = self._molecule_map[item]

        # Log multiplicity at point of use so any default is visible
        _meta_path = os.path.join(self.metadata_dir, f"{item}.json")
        if os.path.exists(_meta_path):
            try:
                with open(_meta_path) as _mf:
                    _m = json.load(_mf)
                _mult     = int(_m.get("multiplicity", 1))
                _mult_src = _m.get("multiplicity_source", "default")
                if _mult_src == "default":
                    self.log_warning(
                        f"[COSMO-RS] {item}: multiplicity defaulted to 1 "
                        f"— not provided and cannot be reliably determined"
                    )
                else:
                    self.log_info(
                        f"[COSMO-RS] {item}: multiplicity={_mult} "
                        f"(source: {_mult_src})"
                    )
            except Exception:
                pass

        return {
            **self._base_args(item),
            "inchi_key":    item,
            "conformers":   mol["conformers"],
            "valid_combos": self.valid_combos,
            "solvent_bank": str(self.solvent_bank),
            "metadata_dir": self.metadata_dir,
            "defaults":     self._sol_defaults,
            "opencosmo":    self.opencosmo,
            "verbose":      self.verbose,
        }

    # =========================================================================
    # Assemble output from checkpoints
    # =========================================================================

    def _assemble_output(self):
        # Load per-molecule result JSONs (written by worker into results/)
        # rather than the raw checkpoints — these have the full mol_result
        # structure including solvent_results.
        all_results = []
        for mol_json in sorted(Path(self.results_dir).glob("*.json")):
            try:
                with open(mol_json) as f:
                    all_results.append(json.load(f))
            except Exception as e:
                self.log_warning(f"Could not read {mol_json.name}: {e}")

        if not all_results:
            self.log_warning("No results to assemble")
            return

        # ── solubility_results.json ───────────────────────────────────────────
        results_path = self.get_stage_output()
        with AtomicWriter(results_path) as f:
            json.dump(all_results, f, indent=2)
        self.log_info(
            f"Wrote solubility_results.json — {len(all_results)} molecules"
        )

        self._write_csv(all_results)
        self._write_human_summary(all_results)

    # =========================================================================
    # CSV output
    # =========================================================================

    def _write_csv(self, all_results: list):
        csv_path = os.path.join(self.outputs_dir, "solubility_results.csv")
        fieldnames = [
            "inchi_key", "mol_name", "smiles",
            "melting_temp", "melting_temp_source",
            "experimental_solubility_mol_frac", "log_solubility_experimental",
            "aqsol_predicted_solubility_mol_frac",
            "combo_id", "combo_label",
            "temperature_K", "predicted_solubility", "log_solubility_predicted",
            "abs_error", "logS_abs_error",
            "n_solute_confs", "n_solvent_confs",
            "sorcf_used",
            "error", "mol_json_path",
        ]
        rows = []
        for mol in all_results:
            ik       = mol["inchi_key"]
            mol_json = os.path.join(self.results_dir, f"{ik}.json")
            exp      = mol.get("experimental_solubility_mol_frac")

            for sr in mol.get("solvent_results", []):
                pred      = sr.get("predicted_solubility")
                abs_error = (
                    round(abs(pred - exp), 8)
                    if pred is not None and exp is not None
                    else None
                )
                n_solv_str = ",".join(
                    f"{k}:{v}" for k, v in sr.get("n_solvent_confs", {}).items()
                )
                logS_exp  = mol.get("log_solubility_experimental")
                logS_pred = sr.get("log_solubility_predicted")
                logS_err  = (
                    round(abs(logS_pred - logS_exp), 6)
                    if logS_pred is not None and logS_exp is not None
                    else None
                )
                rows.append({
                    "inchi_key":                          ik,
                    "mol_name":                           mol.get("mol_name"),
                    "smiles":                             mol.get("smiles"),
                    "melting_temp":                       mol.get("melting_temp"),
                    "melting_temp_source":                mol.get("melting_temp_source"),
                    "experimental_solubility_mol_frac":   exp,
                    "log_solubility_experimental":        logS_exp,
                    "aqsol_predicted_solubility_mol_frac":mol.get("aqsol_predicted_solubility_mol_frac"),
                    "combo_id":                           sr.get("combo_id"),
                    "combo_label":                        sr.get("combo_label"),
                    "temperature_K":                      sr.get("temperature_K"),
                    "predicted_solubility":               pred,
                    "log_solubility_predicted":           logS_pred,
                    "abs_error":                          abs_error,
                    "logS_abs_error":                     logS_err,
                    "n_solute_confs":                     sr.get("n_solute_confs"),
                    "n_solvent_confs":                    n_solv_str,
                    "sorcf_used":                         sr.get("sorcf_used"),
                    "error":                              sr.get("error"),
                    "mol_json_path":                      mol_json,
                })

        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        self.log_info(f"Wrote solubility_results.csv — {len(rows)} rows")

    # =========================================================================
    # Human summary
    # =========================================================================

    def _write_human_summary(self, all_results: list):
        summary_path = os.path.join(
            self.outputs_dir, "solubility_human_summary.txt"
        )
        sep  = "  "
        cols = [
            ("inchi_key",    27),
            ("mol_name",     22),
            ("smiles",       30),
            ("exp_sol",      12),
            ("combo_id",     22),
            ("combo_label",  25),
            ("pred_sol",     14),
            ("abs_error",    12),
            ("Tm",            9),
            ("Tm_source",    18),
            ("n_sol_confs",  11),
            ("n_solv_confs", 14),
            ("sorcf_used",    9),
        ]

        def _cell(val, w):
            s = str(val) if val is not None else "N/A"
            return s[:w].ljust(w)

        header  = sep.join(h.ljust(w) for h, w in cols) + sep + "json_path"
        divider = "─" * len(header)
        lines   = [header, divider]
        n_total = n_ok = 0

        for mol in all_results:
            ik       = mol["inchi_key"]
            mol_json = os.path.join(self.results_dir, f"{ik}.json")
            exp      = mol.get("experimental_solubility_mol_frac")

            for sr in mol.get("solvent_results", []):
                pred  = sr.get("predicted_solubility")
                n_total += 1
                if pred is not None:
                    n_ok += 1

                abs_err_str = (
                    f"{abs(pred - exp):.4e}"
                    if pred is not None and exp is not None
                    else "N/A"
                )
                n_solv_str  = ",".join(
                    f"{k}:{v}" for k, v in sr.get("n_solvent_confs", {}).items()
                )
                pred_str    = f"{pred:.6e}" if pred is not None else "FAILED"

                values = [
                    (ik,                                27),
                    (mol.get("mol_name", ""),           22),
                    (mol.get("smiles", ""),             30),
                    (exp,                               12),
                    (sr.get("combo_id", ""),            22),
                    (sr.get("combo_label", ""),         25),
                    (pred_str,                          14),
                    (abs_err_str,                       12),
                    (mol.get("melting_temp", ""),        9),
                    (mol.get("melting_temp_source", ""), 18),
                    (sr.get("n_solute_confs", ""),      11),
                    (n_solv_str,                        14),
                    (sr.get("sorcf_used", ""),           9),
                ]
                lines.append(
                    sep.join(_cell(v, w) for v, w in values) + sep + mol_json
                )

        lines += [
            "",
            f"Total:      {n_total}",
            f"Successful: {n_ok}",
            f"Failed:     {n_total - n_ok}",
        ]

        with open(summary_path, "w") as f:
            f.write("\n".join(lines) + "\n")

        self.log_info(
            f"Wrote solubility_human_summary.txt — {n_total} rows"
        )