#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modules/stages/optimisation_stage.py

Geometry optimisation stage.

=========================================================================
INPUT
=========================================================================

  stage_input : energies.json from GenerationStage (or PruningStage)
    Each entry is a conformer record dict:
      lookup_id, inchi_key, conf_num, xyz_path, energy, smiles,
      provenance, optimisation_history

=========================================================================
OUTPUT  (canonical stage output)
=========================================================================

  energies.json  (ConformerSet JSON)
    One entry per conformer — the input record with optimisation_history
    extended by one entry:
    {
        "lookup_id":            str,           {inchi_key}_conf{conf_num:03d}
        "inchi_key":            str,
        "conf_num":             int | null,
        "smiles":               str,
        "energy":               float,         final optimised energy
        "xyz_path":             str,           path to optimised XYZ
        "provenance":           dict,
        "optimisation_history": [
            {
                "stage":           "optimisation",
                "engine":          str,
                "level":           str,
                "status":          "converged" | "partial" | "failed",
                "energy":          float | null,
                "xyz_path":        str | null,
                "log_path":        str,
                "elapsed_seconds": float,
                "timing":          { wall_seconds, cpu_seconds, ... },
                "backend_meta":    { run_status, version, command_str },
                "timestamp":       str  ISO-8601
            },
            ...
        ]
    }

=========================================================================
AUXILIARY OUTPUTS
=========================================================================

  checkpoints/<lookup_id>.json     written atomically after each success
  xyz/<lookup_id>_opt<N>.xyz       optimised geometry per iteration
  log/<lookup_id>_opt<N>.log       backend log per iteration
  optimisation_summary.csv         human-readable summary

=========================================================================
DEFAULTS  (optimisation_defaults.json)
=========================================================================

  {
    "engine":                 "gxtb_opt_normal",
    "level":                  "normal",
    "max_iter":               250,
    "keep_scratch":           false,
    "global_fail_threshold":  0.8,
    "strict":                 false
  }

=========================================================================
ENGINES  (optimisation_engines.json)
=========================================================================

  See optimisation_engines.json — entries cover:
    orca_opt_*, orca_sp_*, orca_xtb2_alpb_opt
    xtb_opt_{loose,normal,tight,vtight}
    gxtb_opt_{loose,normal,tight,vtight}
    forcefield_{mmff,uff}

=========================================================================
PARALLELISM
=========================================================================

  Uses BaseStage.run_parallel().
  n_workers and cores_per_item injected by PipelineRunner.

  ORCA backend:  cores_per_item → %pal nprocs
  XTB backend:   cores_per_item → OMP_NUM_THREADS
  gXTB backend:  cores_per_item → OMP_NUM_THREADS

=========================================================================
STRICT MODE HIERARCHY
=========================================================================

  1. parameters["strict"]                   runtime override
  2. optimisation_defaults.json "strict"    stage default
  3. False                                  fallback

=========================================================================
GLOBAL FAILURE THRESHOLD
=========================================================================

  After the pool drains, if the fraction of failed items exceeds
  global_fail_threshold the stage raises an error.  This catches
  systemic issues that would otherwise produce a near-empty output.

  Distinct from strict_mode (which aborts mid-pool on first failure).

"""

import json
import math
import multiprocessing
import os
import shutil
import subprocess
import time
from datetime import datetime

import pandas as pd

try:
    import rdkit
    from rdkit.Chem import AllChem
    from modules.utils.rdkit_helper import load_mol_from_xyz as _ff_load_mol
    from modules.utils.rdkit_helper import write_xyz as _ff_write_xyz
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter
from modules.utils.conformers import ConformerRecord, ConformerSet

HARTREE_TO_KCAL = 627.509474   # all stored energies are in kcal/mol

# ── Parser registry (imported lazily to keep worker picklable) ────────────────
# Parsers are called inside the worker — imports happen at module load time
# which is fine since this module is imported in the main process before
# the pool is created.

from modules.parsers.opt.forcefield_log_parser import ForcefieldLogParser
from modules.parsers.opt.gxtb_log_parser import GxTBLogParser
from modules.parsers.opt.orca_log_parser import ORCALogParser
from modules.parsers.opt.xtb_log_parser import XTBLogParser

_PARSERS = {
    "orca":        ORCALogParser,
    "xtb":         XTBLogParser,
    "gxtb":        GxTBLogParser,
    "forcefield":  ForcefieldLogParser,
}

_ORCA_OPT_KEYWORDS = {
    "loose":  "LooseOpt",
    "normal": "Opt",
    "tight":  "TightOpt",
    "vtight": "VeryTightOpt",
}

VALID_LEVELS = ("loose", "normal", "tight", "vtight")


# =============================================================================
# Module-level worker — must be at module scope to be picklable
# =============================================================================

def _optimisation_worker(args: dict) -> dict:
    """
    Pure worker function.  Receives a fully self-contained args dict,
    runs the appropriate backend, parses the log, updates optimisation
    history, and returns a checkpoint dict.

    Contract fields in return dict:
      status        "ok" | "failed"
      log_details   list of (key, value) tuples  → [COMPLETE] log rows
      error         str                           → [FAILED] log row
      worker_name   str
      worker_pid    int
    """
    worker_name = multiprocessing.current_process().name
    worker_pid  = os.getpid()

    item_id        = args["item_id"]
    charge         = args.get("charge", 0)
    multiplicity   = args.get("multiplicity", 1)
    lookup_id      = args["lookup_id"]
    workdir        = args["workdir"]
    outputs_dir    = args["outputs_dir"]
    cores_per_item = args["cores_per_item"]
    engine_name    = args["engine_name"]
    engine_spec    = args["engine_spec"]
    level          = args["level"]
    max_iter       = args["max_iter"]
    keep_scratch   = args["keep_scratch"]
    config         = args["config"]
    record         = args["record"]      # conformer record as plain dict

    xyz_dir     = os.path.join(outputs_dir, "xyz")
    log_dir     = os.path.join(outputs_dir, "log")
    scratch_dir = os.path.join(workdir, "scratch")

    for d in (workdir, xyz_dir, log_dir, scratch_dir):
        os.makedirs(d, exist_ok=True)

    # Iteration number = length of existing history
    iteration  = len(record.get("optimisation_history", []))
    input_xyz  = record["xyz_path"]
    output_xyz = os.path.join(xyz_dir, f"{lookup_id}_opt{iteration}.xyz")
    output_log = os.path.join(log_dir, f"{lookup_id}_opt{iteration}.log")

    try:
        t_wall_start = time.perf_counter()
        t_cpu_start  = time.process_time()

        backend_result = _run_backend(
            input_xyz      = input_xyz,
            output_xyz     = output_xyz,
            output_log     = output_log,
            scratch_dir    = scratch_dir,
            engine_name    = engine_name,
            engine_spec    = engine_spec,
            level          = level,
            max_iter       = max_iter,
            cores_per_item = cores_per_item,
            config         = config,
            charge         = charge,
            multiplicity   = multiplicity,
        )

        t_wall_end = time.perf_counter()
        t_cpu_end  = time.process_time()
        elapsed    = t_wall_end - t_wall_start

        # Parse log, determine status, validate energy
        parsed = _parse_backend_log(output_log, engine_spec)
        status = _determine_status(parsed, output_xyz)
        energy = _validate_energy(parsed.get("energy"))

        # Normalise to kcal/mol — quantum backends (xtb, gxtb, orca) output
        # Hartree; forcefield backends (MMFF94/UFF) already output kcal/mol.
        if energy is not None and engine_spec["family"] != "forcefield":
            energy = energy * HARTREE_TO_KCAL

        # Build history entry — schema unchanged from original stage
        history_entry = {
            "stage":           "optimisation",
            "engine":          engine_name,
            "level":           level,
            "status":          status,
            "energy":          energy,
            "xyz_path":        output_xyz if os.path.exists(output_xyz) else None,
            "log_path":        output_log,
            "elapsed_seconds": elapsed,
            "timing": {
                "start_wall":   t_wall_start,
                "end_wall":     t_wall_end,
                "start_cpu":    t_cpu_start,
                "end_cpu":      t_cpu_end,
                "wall_seconds": t_wall_end - t_wall_start,
                "cpu_seconds":  t_cpu_end  - t_cpu_start,
            },
            "backend_meta": {
                "run_status":  backend_result.get("run_status"),
                "version":     backend_result.get("version", "unknown"),
                "command_str": backend_result.get("command_str", ""),
            },
            "timestamp": datetime.utcnow().isoformat() + "Z",
        }

        # Build updated record dict
        history = list(record.get("optimisation_history", []))
        history.append(history_entry)

        updated = dict(record)
        updated["optimisation_history"] = history
        updated["energy"]               = energy
        updated["xyz_path"]             = (
            output_xyz if os.path.exists(output_xyz) else None
        )

        if not keep_scratch:
            shutil.rmtree(scratch_dir, ignore_errors=True)

        energy_str = f"{energy:.4f} kcal/mol" if energy is not None else "n/a"

        return {
            # BaseStage contract
            "item_id":     item_id,
            "status":      "ok",
            "worker_name": worker_name,
            "worker_pid":  worker_pid,
            "log_details": [
                ("engine",  engine_name),
                ("level",   level),
                ("result",  status),
                ("energy",  energy_str),
                ("elapsed", f"{elapsed:.1f}s"),
            ],
            # Stage payload — full updated conformer record
            **updated,
        }

    except Exception as e:
        if not keep_scratch:
            shutil.rmtree(scratch_dir, ignore_errors=True)
        return {
            "item_id":     item_id,
            "status":      "failed",
            "error":       str(e),
            "worker_name": worker_name,
            "worker_pid":  worker_pid,
            "log_details": [],
        }


# =============================================================================
# Backend dispatcher
# =============================================================================

def _run_backend(
    input_xyz:      str,
    output_xyz:     str,
    output_log:     str,
    scratch_dir:    str,
    engine_name:    str,
    engine_spec:    dict,
    level:          str,
    max_iter:       int,
    cores_per_item: int,
    config:         dict,
    charge:         int = 0,
    multiplicity:   int = 1,
) -> dict:
    family = engine_spec["family"]

    if family == "orca":
        return _backend_orca(
            input_xyz      = input_xyz,
            output_xyz     = output_xyz,
            output_log     = output_log,
            scratch_dir    = scratch_dir,
            engine_spec    = engine_spec,
            level          = level,
            cores_per_item = cores_per_item,
            config         = config,
            charge         = charge,
            multiplicity   = multiplicity,
        )
    if family == "gxtb":
        return _backend_gxtb(
            input_xyz      = input_xyz,
            output_xyz     = output_xyz,
            output_log     = output_log,
            scratch_dir    = scratch_dir,
            engine_spec    = engine_spec,
            level          = level,
            max_iter       = max_iter,
            cores_per_item = cores_per_item,
            config         = config,
            charge         = charge,
            multiplicity   = multiplicity,
        )
    if family == "xtb":
        return _backend_xtb(
            input_xyz      = input_xyz,
            output_xyz     = output_xyz,
            output_log     = output_log,
            scratch_dir    = scratch_dir,
            engine_spec    = engine_spec,
            level          = level,
            max_iter       = max_iter,
            cores_per_item = cores_per_item,
            config         = config,
            charge         = charge,
            multiplicity   = multiplicity,
        )
    if family == "forcefield":
        return _backend_forcefield(
            input_xyz   = input_xyz,
            output_xyz  = output_xyz,
            output_log  = output_log,
            engine_spec = engine_spec,
            max_iter    = max_iter,
            charge      = charge,
        )

    raise RuntimeError(f"Unknown engine family: {family}")


# =============================================================================
# XYZ helpers
# =============================================================================

def _extract_last_xyz_frame(trj_path, out_path):
    """
    Extract the last frame from an XYZ trajectory file into out_path.
    Returns True on success, False if the file is empty or malformed.
    Used as a fallback when ORCA writes only a trajectory (_trj.xyz).
    """
    try:
        with open(trj_path) as fh:
            raw = fh.read()
        lines  = raw.splitlines()
        frames = []
        i = 0
        while i < len(lines):
            stripped = lines[i].strip()
            if not stripped:
                i += 1
                continue
            try:
                n_atoms = int(stripped)
            except ValueError:
                i += 1
                continue
            end = i + n_atoms + 2
            if end > len(lines):
                break
            frames.append(lines[i:end])
            i = end
        if not frames:
            return False
        with open(out_path, "w") as fh:
            fh.write("\n".join(frames[-1]) + "\n")
        return True
    except Exception:
        return False


# =============================================================================
# Backend: ORCA
# =============================================================================

def _backend_orca(
    input_xyz:      str,
    output_xyz:     str,
    output_log:     str,
    scratch_dir:    str,
    engine_spec:    dict,
    level:          str,
    cores_per_item: int,
    config:         dict,
    charge:         int = 0,
    multiplicity:   int = 1,
) -> dict:
    # Short stable base name avoids filesystem path-length limits when
    # InChIKeys are used.  The per-item scratch_dir provides uniqueness.
    lookup   = "mol"
    inp_file = os.path.join(scratch_dir, lookup + ".inp")
    orca_exe = config["orca"]["executable"]

    opt       = engine_spec.get("opt",    True)
    sp        = engine_spec.get("sp",     False)
    cpcm_mode = engine_spec.get("cpcm",   "none")
    alpb      = engine_spec.get("alpb",   False)
    method    = engine_spec.get("method")
    basis     = engine_spec.get("basis")
    solvent   = engine_spec.get("solvent", "water")

    lines = ["%MaxCore 2000", ""]

    if cores_per_item > 1:
        lines += [f"%pal nprocs {cores_per_item}", "end", ""]

    # CPCM block
    if cpcm_mode == "default":
        lines += ["%cpcm", "end", ""]
    elif cpcm_mode == "custom":
        chem_dir  = config["constant_files"]["chemistry_dir"]
        cpcm_path = os.path.join(chem_dir, "cpcm_radii.json")
        with open(cpcm_path) as f:
            cpcm_cfg = json.load(f)
        lines.append("%cpcm")
        for Z, r in cpcm_cfg["radii"].items():
            lines.append(f"  radius[{Z}] {r}")
        lines += [f"  cut_area {cpcm_cfg['cut_area']}", "end", ""]

    # Method line
    if alpb:
        opt_kw      = _ORCA_OPT_KEYWORDS.get(level, "Opt") if (opt and not sp) else "SP"
        method_line = f"! XTB2 {opt_kw} ALPB({solvent})"
    else:
        opt_kw      = _ORCA_OPT_KEYWORDS.get(level, "Opt") if (opt and not sp) else "SP"
        method_line = (
            f"! {method} {basis} {opt_kw}"
            if basis else
            f"! {method} {opt_kw}"
        )

    lines += [method_line, "", f'%base "{lookup}"', "", f"* xyz {charge} {multiplicity}"]

    with open(input_xyz) as xyz_f:
        xyz_lines = xyz_f.readlines()
    lines.extend(l.rstrip() for l in xyz_lines[2:])
    lines += ["*", ""]

    with open(inp_file, "w") as f:
        f.write("\n".join(lines))

    cmd        = [orca_exe, inp_file]
    run_status = "ok"

    try:
        with open(output_log, "w") as log:
            subprocess.run(
                cmd, cwd=scratch_dir,
                stdout=log, stderr=log, check=True,
            )
        final_xyz = os.path.join(scratch_dir, lookup + ".xyz")
        trj_xyz   = os.path.join(scratch_dir, lookup + "_trj.xyz")
        if os.path.exists(final_xyz):
            shutil.copy(final_xyz, output_xyz)
        elif os.path.exists(trj_xyz):
            # Some ORCA jobs write only a trajectory; extract last frame
            if not _extract_last_xyz_frame(trj_xyz, output_xyz):
                run_status = "failed"
        else:
            run_status = "failed"
    except subprocess.CalledProcessError:
        run_status = "failed"

    return {
        "run_status":  run_status,
        "version":     "unknown",
        "command_str": " ".join(cmd),
    }


# =============================================================================
# Backend: gXTB
# =============================================================================

def _backend_gxtb(
    input_xyz:      str,
    output_xyz:     str,
    output_log:     str,
    scratch_dir:    str,
    engine_spec:    dict,
    level:          str,
    max_iter:       int,
    cores_per_item: int,
    config:         dict,
    charge:         int = 0,
    multiplicity:   int = 1,
) -> dict:
    xtb_bin  = config["xtb"]["executable"]
    gxtb_bin = config["gxtb"]["executable"]

    tmp_xyz = os.path.join(scratch_dir, "molecule.xyz")
    opt_tmp = os.path.join(scratch_dir, "xtbopt.xyz")

    for p in (tmp_xyz, opt_tmp):
        if os.path.exists(p):
            os.remove(p)

    shutil.copy(input_xyz, tmp_xyz)

    # Use opt_flag from engine_spec if present (e.g. "--opt loose" → ["--opt", "loose"])
    # Fall back to level-based dict if engine_spec doesn't have it
    _opt_flag_str = engine_spec.get("opt_flag")
    if _opt_flag_str:
        opt_flag = _opt_flag_str.split()
    else:
        opt_flag = {
            "loose":  ["--opt", "loose"],
            "normal": ["--opt"],
            "tight":  ["--opt", "tight"],
            "vtight": ["--opt", "vtight"],
        }.get(level, ["--opt"])

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(cores_per_item)

    driver_string = f"{gxtb_bin} -grad -c xtbdriver.xyz"

    cmd = [
        xtb_bin, "molecule.xyz",
        "--driver", driver_string,
        *opt_flag,
        "--iterations", str(max_iter),
        "--chrg", str(charge),
        "--uhf",  str(multiplicity - 1),
    ]

    run_status = "ok"
    try:
        with open(output_log, "w") as log:
            subprocess.run(
                cmd, cwd=scratch_dir,
                stdout=log, stderr=log,
                text=True, check=False, env=env,
            )
        if os.path.exists(opt_tmp):
            shutil.copy(opt_tmp, output_xyz)
        else:
            run_status = "failed"
    except Exception:
        run_status = "failed"

    return {
        "run_status":  run_status,
        "version":     "unknown",
        "command_str": " ".join(cmd),
    }


# =============================================================================
# Backend: XTB
# =============================================================================

def _backend_xtb(
    input_xyz:      str,
    output_xyz:     str,
    output_log:     str,
    scratch_dir:    str,
    engine_spec:    dict,
    level:          str,
    max_iter:       int,
    cores_per_item: int,
    config:         dict,
    charge:         int = 0,
    multiplicity:   int = 1,
) -> dict:
    xtb_bin = config["xtb"]["executable"]
    gfn     = engine_spec.get("gfn", 2)

    tmp_xyz = os.path.join(scratch_dir, "input.xyz")
    opt_tmp = os.path.join(scratch_dir, "xtbopt.xyz")

    for p in (tmp_xyz, opt_tmp):
        if os.path.exists(p):
            os.remove(p)

    shutil.copy(input_xyz, tmp_xyz)

    # Level-based iteration count — XTB converges better with more iterations
    # at tighter levels. max_iter is the fallback for unrecognised levels.
    level_iters = {
        "loose":  100,
        "normal": 250,
        "tight":  500,
        "vtight": 1000,
    }
    iters = level_iters.get(level, max_iter)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(cores_per_item)

    cmd = [
        xtb_bin, tmp_xyz,
        "--opt",
        "--gfn", str(gfn),
        "--iterations", str(iters),
        "--chrg", str(charge),
        "--uhf",  str(multiplicity - 1),
    ]

    run_status = "ok"
    try:
        with open(output_log, "w") as log:
            subprocess.run(
                cmd, cwd=scratch_dir,
                stdout=log, stderr=log,
                check=False, env=env,
            )
        if os.path.exists(opt_tmp):
            shutil.copy(opt_tmp, output_xyz)
        else:
            run_status = "failed"
    except Exception:
        run_status = "failed"

    return {
        "run_status":  run_status,
        "version":     "unknown",
        "command_str": " ".join(cmd),
    }


# =============================================================================
# Backend: Forcefield (RDKit MMFF94 / UFF)
# =============================================================================


def _backend_forcefield(
    input_xyz:   str,
    output_xyz:  str,
    output_log:  str,
    engine_spec: dict,
    max_iter:    int = 500,
    charge:      int = 0,
) -> dict:
    if not _RDKIT_AVAILABLE:
        raise RuntimeError(
            "RDKit is not installed — cannot run forcefield optimisation. "
            "Install it with: conda install -c conda-forge rdkit"
        )

    ff_name = engine_spec.get("forcefield", "MMFF94").upper()
    use_mmff = ff_name.startswith("MMFF")

    try:
        mol = _ff_load_mol(input_xyz, charge=charge)

        if use_mmff:
            props = AllChem.MMFFGetMoleculeProperties(mol)
            ff    = AllChem.MMFFGetMoleculeForceField(mol, props)
        else:
            ff = AllChem.UFFGetMoleculeForceField(mol)

        ff.Initialize()
        ret_code = ff.Minimize(maxIts=max_iter)
        # ret_code: 0 = converged, 1 = not converged within maxIts, -1 = error
        opt_status = "converged" if ret_code == 0 else "not_converged"
        energy     = ff.CalcEnergy()

        comment = (
            f"energy: {energy:.10f} kcal/mol  ff: {ff_name}"
            f"  status: {opt_status}  rdkit: {rdkit.__version__}"
        )
        _ff_write_xyz(mol, output_xyz, comment=comment)

        with open(output_log, "w") as f:
            f.write(f"FORCEFIELD_FAMILY  {ff_name}\n")
            f.write(f"OPTIMIZATION_STATUS  {opt_status}\n")
            f.write(f"FORCEFIELD_ENERGY  {energy:.10f}\n")
            f.write(f"ATOMS  {mol.GetNumAtoms()}\n")

        return {
            "run_status":  "ok",
            "version":     rdkit.__version__,
            "command_str": "",
        }

    except Exception as exc:
        with open(output_log, "w") as f:
            f.write(f"FORCEFIELD_FAMILY  {ff_name}\n")
            f.write(f"OPTIMIZATION_STATUS  failed\n")
            f.write(f"ERROR  {exc}\n")
        return {
            "run_status":  "failed",
            "version":     rdkit.__version__ if _RDKIT_AVAILABLE else "unknown",
            "command_str": "",
        }


# =============================================================================
# Log parsing / status helpers
# =============================================================================

def _parse_backend_log(log_path: str, engine_spec: dict) -> dict:
    family     = engine_spec["family"]
    parser_cls = _PARSERS.get(family)

    if not parser_cls or not log_path or not os.path.exists(log_path):
        return {"energy": None, "converged": None}
    try:
        return parser_cls.parse(log_path)
    except Exception:
        return {"energy": None, "converged": None}


def _determine_status(parsed: dict, xyz_path: str) -> str:
    conv       = parsed.get("converged")
    xyz_exists = bool(xyz_path and os.path.exists(xyz_path))

    if conv is True:
        return "converged" if xyz_exists else "partial"
    if conv is False:
        return "failed"
    return "partial"


def _validate_energy(energy) -> float | None:
    if energy is None:
        return None
    try:
        val = float(energy)
        if math.isnan(val) or abs(val) > 1e6:
            return None
        return val
    except Exception:
        return None


def _is_record_usable(record: dict) -> bool:
    history = record.get("optimisation_history", [])
    if not history:
        return False
    last   = history[-1]
    xyz    = last.get("xyz_path")
    energy = last.get("energy")
    if not xyz or not os.path.exists(xyz):
        return False
    if energy is None or math.isnan(energy) or abs(energy) > 1e6:
        return False
    return True


# =============================================================================
# Stage class
# =============================================================================

class OptimisationStage(BaseStage):
    """
    Geometry optimisation stage.

    execute() loads config and conformer records, delegates per-item work
    to _optimisation_worker via BaseStage.run_parallel(), then assembles
    energies.json from checkpoints.

    Stages never call mark_complete() — BaseStage.run() owns lifecycle.
    """

    # =========================================================================
    # Entry point
    # =========================================================================

    def execute(self):
        self.set_stage_output("energies.json")

        self._load_stage_config()
        self._load_entries()

        lookup_ids = [r["lookup_id"] for r in self._all_records]

        # Resume vs fresh start
        if self.job.pending_items:
            self.log_info(
                f"Resuming — {len(self.job.pending_items)} pending / "
                f"{len(lookup_ids) - len(self.job.pending_items)} "
                f"already complete"
            )
        else:
            self.set_items(lookup_ids, sort_by_complexity=True)
            self._clear_checkpoints()
            self.log_info(
                f"Fresh run — {len(lookup_ids)} conformers to optimise"
            )

        # Skip items with missing XYZ before entering pool
        for record in self._all_records:
            lid = record["lookup_id"]
            if lid not in self.job.pending_items:
                continue
            xyz = record.get("xyz_path")
            if not xyz or not os.path.exists(xyz):
                self.log_skip(lid, f"XYZ not found: {xyz}")
                self.update_progress(lid, success=False)

        self.run_parallel(_optimisation_worker, self._build_args)
        self._check_global_threshold()
        self._assemble_output()

    # =========================================================================
    # Config
    # =========================================================================

    def _load_stage_config(self):
        cfg = self.config
        if not cfg:
            self.fail("OptimisationStage: missing config")

        # ── Defaults file ─────────────────────────────────────────────────────
        defaults_path = cfg.get("optimisation", {}).get("defaults")
        self._opt_defaults = {}
        if defaults_path and os.path.exists(defaults_path):
            with open(defaults_path) as f:
                self._opt_defaults = json.load(f)
        else:
            self.log_warning(
                "optimisation_defaults.json not found — using built-in defaults"
            )

        # ── Engine registry ───────────────────────────────────────────────────
        engines_path = cfg.get("optimisation", {}).get("engines")
        if not engines_path or not os.path.exists(engines_path):
            self.fail("optimisation_engines.json not found in config")
        with open(engines_path) as f:
            self._engine_registry = json.load(f)

        # ── Merge defaults + runtime parameter overrides ──────────────────────
        # Runtime parameters can override engine, level, max_iter etc.
        merged = dict(self._opt_defaults)
        for key in ("engine", "level", "max_iter", "keep_scratch",
                    "global_fail_threshold"):
            if key in self.parameters:
                merged[key] = self.parameters[key]

        # ── Engine selection ──────────────────────────────────────────────────
        self.engine_name = merged.get("engine", "gxtb_opt_normal")
        if self.engine_name not in self._engine_registry:
            self.fail(
                f"Unknown optimisation engine: '{self.engine_name}'. "
                f"Available: {sorted(self._engine_registry.keys())}"
            )
        self.engine_spec = self._engine_registry[self.engine_name]

        # ── Level — engine spec takes priority over defaults ──────────────────
        self.level = self.engine_spec.get(
            "level", merged.get("level", "normal")
        )
        if self.level not in VALID_LEVELS:
            self.fail(
                f"Invalid optimisation level: '{self.level}'. "
                f"Must be one of: {VALID_LEVELS}"
            )

        # ── Other parameters ──────────────────────────────────────────────────
        self.max_iter              = int(merged.get("max_iter", 250))
        self.keep_scratch          = bool(merged.get("keep_scratch", False))
        self.global_fail_threshold = float(
            merged.get("global_fail_threshold", 0.8)
        )

        # ── Strict mode hierarchy ─────────────────────────────────────────────
        # 1. parameters["strict"]    runtime override
        # 2. defaults["strict"]      stage default
        # 3. False                   fallback
        params_strict = self.parameters.get("strict")
        if params_strict is not None:
            self.strict_mode = bool(params_strict)
        else:
            self.strict_mode = bool(self._opt_defaults.get("strict", False))

        # ── Log resolved config ───────────────────────────────────────────────
        self.log_config(
            f"engine={self.engine_name}  "
            f"family={self.engine_spec['family']}  "
            f"level={self.level}  "
            f"max_iter={self.max_iter}"
        )
        self.log_config(
            f"keep_scratch={self.keep_scratch}  "
            f"global_fail_threshold={self.global_fail_threshold}  "
            f"strict={self.strict_mode}"
        )
        self.log_resources(
            f"cores_per_item={self.cores_per_item}  "
            f"n_workers={self.n_workers}  "
            f"budget={self.budget_cores}"
        )

    # =========================================================================
    # Load conformer records
    # =========================================================================

    def _load_entries(self):
        energies_file = self.require_file(
            self.get_stage_input(), "stage_input energies.json"
        )
        self.log_info(f"Stage input: {energies_file}")

        # Copy into inputs/ for reproducibility
        os.makedirs(self.inputs_dir, exist_ok=True)
        inputs_copy = os.path.join(self.inputs_dir, "energies.json")
        if os.path.abspath(energies_file) != os.path.abspath(inputs_copy):
            shutil.copy(energies_file, inputs_copy)

        # Load via ConformerSet then convert to plain dicts
        conformer_set   = ConformerSet.load(energies_file)
        self._all_records = []

        for rec in conformer_set.records:
            # Ensure optimisation_history is present
            if not hasattr(rec, "optimisation_history"):
                rec.optimisation_history = []

            # Convert to plain dict for pickling
            self._all_records.append(_conformer_to_dict(rec))

        # Copy XYZs into inputs/xyz for reproducibility
        xyz_dir = os.path.join(self.inputs_dir, "xyz")
        os.makedirs(xyz_dir, exist_ok=True)
        copied = 0
        for record in self._all_records:
            src = record.get("xyz_path")
            if src and os.path.exists(src):
                dst = os.path.join(xyz_dir, os.path.basename(src))
                if os.path.abspath(src) != os.path.abspath(dst):
                    shutil.copy(src, dst)
                    record["xyz_path"] = dst
                copied += 1

        self.log_info(
            f"Loaded {len(self._all_records)} conformer records "
            f"({copied} XYZ files staged)"
        )

    # =========================================================================
    # Build args for worker
    # =========================================================================

    def _build_args(self, item: str) -> dict:
        """
        Build self-contained args dict for _optimisation_worker.
        Extends _base_args with all backend-specific data.
        Everything must be JSON-serialisable / picklable.
        """
        record = self._record_map[item]
        provenance   = record.get("provenance", {})
        charge       = int(provenance.get("charge", 0))
        multiplicity = int(provenance.get("multiplicity", 1))
        return {
            **self._base_args(item),
            "lookup_id":    item,
            "engine_name":  self.engine_name,
            "engine_spec":  self.engine_spec,
            "level":        self.level,
            "max_iter":     self.max_iter,
            "keep_scratch": self.keep_scratch,
            "config":       self.config,
            "record":       record,
            "charge":       charge,
            "multiplicity": multiplicity,
        }

    # =========================================================================
    # Global failure threshold check
    # =========================================================================

    def _check_global_threshold(self):
        """
        After the pool drains, raise if too many items failed.
        Distinct from strict_mode — this detects systemic failure.
        """
        n_total  = len(self.job.items)
        n_failed = len(self.job.failed_items)

        if n_total == 0:
            return

        fail_rate = n_failed / n_total
        if fail_rate >= self.global_fail_threshold:
            self.fail(
                f"Global failure threshold exceeded: "
                f"{n_failed}/{n_total} items failed "
                f"({fail_rate:.0%} >= threshold {self.global_fail_threshold:.0%})"
            )

    # =========================================================================
    # Assemble output from checkpoints
    # =========================================================================

    def _assemble_output(self):
        checkpoints = self._load_all_checkpoints()

        if not checkpoints:
            self.fail("No successful checkpoints — cannot write energies.json")

        # Filter to usable records only
        usable = [cp for cp in checkpoints if _is_record_usable(cp)]

        if not usable:
            self.fail(
                "Optimisation produced zero usable conformers "
                "(all lack valid XYZ or energy)"
            )

        self.log_info(
            f"Assembling output: {len(usable)} usable / "
            f"{len(checkpoints)} total checkpoints"
        )

        # Reconstruct ConformerSet from checkpoint dicts
        records = [_dict_to_conformer(cp) for cp in usable]
        final_set = ConformerSet(records)
        final_set.save(self.get_stage_output())

        self.log_info(
            f"Wrote energies.json — {len(records)} conformers"
        )

        # Write CSV summary
        self._write_summary_csv(usable)

    # =========================================================================
    # CSV summary
    # =========================================================================

    def _write_summary_csv(self, checkpoints: list):
        csv_path = os.path.join(self.outputs_dir, "optimisation_summary.csv")
        rows = []
        for cp in checkpoints:
            history = cp.get("optimisation_history", [])
            last    = history[-1] if history else {}
            rows.append({
                "lookup_id":       cp.get("lookup_id"),
                "inchi_key":       cp.get("inchi_key"),
                "conf_num":        cp.get("conf_num"),
                "smiles":          cp.get("smiles"),
                "energy":          last.get("energy"),
                "status":          last.get("status"),
                "engine":          last.get("engine"),
                "level":           last.get("level"),
                "elapsed_seconds": last.get("elapsed_seconds"),
                "xyz_path":        last.get("xyz_path"),
            })

        if rows:
            pd.DataFrame(rows).to_csv(csv_path, index=False)
            self.log_info(f"Wrote optimisation_summary.csv — {len(rows)} rows")

    # =========================================================================
    # Record map (lazy — built on first access)
    # =========================================================================

    @property
    def _record_map(self) -> dict:
        if not hasattr(self, "_record_map_cache"):
            self._record_map_cache = {
                r["lookup_id"]: r for r in self._all_records
            }
        return self._record_map_cache


# =============================================================================
# ConformerRecord ↔ dict helpers
# =============================================================================

def _conformer_to_dict(rec: ConformerRecord) -> dict:
    """
    Serialise a ConformerRecord to a plain dict for pickling.
    Delegates to rec.to_dict() which is the canonical serialisation method.
    lookup_id is computed by ConformerRecord as {inchi_key}_conf{conf_num:03d}.
    """
    return rec.to_dict()


def _dict_to_conformer(d: dict) -> ConformerRecord:
    """
    Reconstruct a ConformerRecord from a checkpoint dict.
    Delegates to ConformerRecord.from_dict() — the canonical deserialiser.
    """
    return ConformerRecord.from_dict(d)