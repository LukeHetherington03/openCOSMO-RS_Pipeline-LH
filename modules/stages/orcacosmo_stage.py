#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modules/stages/orcacosmo_stage.py

ORCA CPCM COSMO stage.

=========================================================================
INPUT
=========================================================================

  stage_input : energies.json from OptimisationStage (or PruningStage)
    Each entry:
      lookup_id, inchi_key, conf_num, xyz_path, energy, smiles,
      provenance, optimisation_history

=========================================================================
OUTPUT  (canonical stage output)
=========================================================================

  orcacosmo_summary.json
    One entry per conformer — full history preserved.
    Schema per entry:
    {
        "item_id":              str,
        "lookup_id":            str,
        "inchi_key":            str,
        "conf_num":             int | null,
        "smiles":               str,
        "energy":               float,
        "provenance":           dict,
        "optimisation_history": list,
        "orcacosmo_history": {
            "method_used":        str,       basis set that succeeded
            "fallback_triggered": bool,
            "elapsed_default":    float,     wall time for default_basis attempt
            "elapsed_fallback":   float | null,
            "elapsed_total":      float,
            "orca_version":       str,
            "raw_output_paths":   { "log": str, "cpcm": str, "cpcm_corr": str }
        },
        "orcacosmo_path":       str,
        "status":               "ok"
    }

=========================================================================
AUXILIARY OUTPUTS
=========================================================================

  orcacosmo_summary.csv
  raw_outputs/<lookup_id>.{log,cpcm,cpcm_corr}
  parsed_outputs/<lookup_id>.{log,cpcm,cpcm_corr}.json
  orcacosmo_outputs/<lookup_id>.orcacosmo
  orca_inputs/<lookup_id>_<basis>.inp
  checkpoints/<lookup_id>.json       written atomically after each success

=========================================================================
DEFAULTS  (orcacosmo_defaults.json)
=========================================================================

  {
    "default_basis":  "TZVP",       basis set to use — "TZVP" or "TZVPD"
    "functional":     "BP86",        DFT functional
    "solvation":      "CPCM",        solvation model keyword
    "maxcore_mb":     2000,          ORCA %MaxCore in MB
    "strict":         false          abort on first failure
  }

  Fallback basis sets are handled by adding a second orcacosmo stage in
  the pipeline spec with the desired basis and parameters.fallback_only=true.
  That stage will skip conformers that already have a valid .orcacosmo file.

=========================================================================
PARALLELISM
=========================================================================

  Uses BaseStage.run_parallel().
  n_workers and cores_per_item injected by PipelineRunner.
  cores_per_item written into ORCA input as %pal nprocs.

=========================================================================
STRICT MODE HIERARCHY
=========================================================================

  1. parameters["strict"]              runtime override
  2. orcacosmo_defaults.json "strict"  stage default
  3. False                             fallback

"""

import json
import multiprocessing
import os
import shutil
import subprocess
import time

import pandas as pd

from modules.stages.base_stage import BaseStage
from modules.parsers.cosmo.orca_log_parser import OrcaLogParser
from modules.parsers.cosmo.orca_cpcm_parser import OrcaCpcmParser
from modules.parsers.cosmo.orca_cpcm_corr_parser import OrcaCpcmCorrParser
from modules.orcacosmo_reconstructors.orcacosmo_orchestrator import OrcaCosmoOrchestrator
from modules.utils.atomic_write import AtomicWriter
from modules.utils.log_helper import LogHelper


# =============================================================================
# Module-level worker — must be at module scope to be picklable
# =============================================================================

def _orcacosmo_worker(args: dict) -> dict:
    """
    Pure worker function.  Receives a fully self-contained args dict,
    runs ORCA, parses outputs, builds the .orcacosmo file, and returns
    a checkpoint dict conforming to the BaseStage contract.

    Contract fields in return dict:
      status        "ok" | "failed"
      log_details   list of (key, value) tuples  → [COMPLETE] log rows
      error         str                           → [FAILED] log row
      worker_name   str
      worker_pid    int
    """
    worker_name = multiprocessing.current_process().name
    worker_pid  = os.getpid()

    # ── Unpack args ───────────────────────────────────────────────────────────
    item_id        = args["item_id"]
    lookup_id      = args["lookup_id"]
    orca_base      = args["orca_base"]       # short safe name for ORCA internals
    xyz_path       = args["xyz_path"]
    workdir        = args["workdir"]
    outputs_dir    = args["outputs_dir"]
    cores_per_item = args["cores_per_item"]
    orca_command   = args["orca_command"]
    orca_version   = args["orca_version"]
    cpcm_radii     = args["cpcm_radii"]
    cpcm_cut_area  = args["cpcm_cut_area"]
    cpcm_file      = args["cpcm_file"]
    default_basis  = args["default_basis"]
    functional     = args["functional"]
    solvation      = args["solvation"]
    maxcore_mb     = args["maxcore_mb"]
    entry          = args["entry"]
    charge         = args.get("charge", 0)
    multiplicity   = args.get("multiplicity", 1)

    # ── Derived paths ─────────────────────────────────────────────────────────
    raw_dir         = os.path.join(outputs_dir, "raw_outputs")
    parsed_dir      = os.path.join(outputs_dir, "parsed_outputs")
    orcacosmo_dir   = os.path.join(outputs_dir, "orcacosmo_outputs")
    orca_inputs_dir = os.path.join(outputs_dir, "orca_inputs")

    for d in (workdir, raw_dir, parsed_dir, orcacosmo_dir, orca_inputs_dir):
        os.makedirs(d, exist_ok=True)

    # Copy XYZ into workdir with the short ORCA base name.
    # _write_orca_input references this as "{orca_base}.xyz" in the
    # xyzfile directive so ORCA can locate it.
    shutil.copy(xyz_path, os.path.join(workdir, f"{orca_base}.xyz"))

    try:
        # ── Run ORCA ──────────────────────────────────────────────────────────
        run_result = _run_orca_with_basis(
            lookup_id       = lookup_id,
            orca_base       = orca_base,
            workdir         = workdir,
            orca_command    = orca_command,
            cpcm_radii      = cpcm_radii,
            cpcm_cut_area   = cpcm_cut_area,
            functional      = functional,
            solvation       = solvation,
            default_basis   = default_basis,
            cores_per_item  = cores_per_item,
            maxcore_mb      = maxcore_mb,
            orca_inputs_dir = orca_inputs_dir,
            charge          = charge,
            multiplicity    = multiplicity,
        )

        method_used        = run_result["method_used"]
        fallback_triggered = run_result["fallback_triggered"]
        elapsed_default    = run_result["elapsed_default"]
        elapsed_fallback   = run_result["elapsed_fallback"]
        elapsed_total      = run_result["elapsed_total"]

        # Copy raw ORCA outputs from short workdir names → external lookup_id names.
        # ORCA writes {orca_base}.log, {orca_base}.cpcm etc. inside workdir.
        # Everything downstream uses lookup_id.
        for ext in ("log", "cpcm", "cpcm_corr"):
            src = os.path.join(workdir, f"{orca_base}.{ext}")
            dst = os.path.join(raw_dir, f"{lookup_id}.{ext}")
            if os.path.exists(src):
                shutil.copy(src, dst)

        # ── Parse ─────────────────────────────────────────────────────────────
        log_json       = _parse_output(lookup_id, raw_dir, "log",       OrcaLogParser)
        cpcm_json      = _parse_output(lookup_id, raw_dir, "cpcm",      OrcaCpcmParser)
        cpcm_corr_json = _parse_output(lookup_id, raw_dir, "cpcm_corr", OrcaCpcmCorrParser)

        # ── Write parser JSONs ────────────────────────────────────────────────
        for kind, data in [
            ("log",       log_json),
            ("cpcm",      cpcm_json),
            ("cpcm_corr", cpcm_corr_json),
        ]:
            out = os.path.join(parsed_dir, f"{lookup_id}.{kind}.json")
            with AtomicWriter(out) as f:
                json.dump(data, f, indent=2)

        # ── Build bundle + reconstruct .orcacosmo ─────────────────────────────
        bundle = {
            "meta": {
                "lookup_id":          lookup_id,
                "inchi_key":          entry["inchi_key"],
                "smiles":             entry.get("smiles"),
                "method_used":        method_used,
                "fallback_triggered": fallback_triggered,
                "functional":         functional,
                "solvation":          solvation,
                "orca_input_path":    os.path.join(
                    orca_inputs_dir,
                    f"{lookup_id}_{method_used.lower()}.inp"
                ),
                "cpcm_radii_source":  cpcm_file,
                "orca_version":       orca_version,
            },
            "paths": {
                "log":       os.path.join(parsed_dir, f"{lookup_id}.log.json"),
                "cpcm":      os.path.join(parsed_dir, f"{lookup_id}.cpcm.json"),
                "cpcm_corr": os.path.join(parsed_dir, f"{lookup_id}.cpcm_corr.json"),
                "xyz":       xyz_path,
            },
            "optimisation_entry": entry,
        }

        bundle_path = os.path.join(parsed_dir, f"{lookup_id}_bundle.json")
        with AtomicWriter(bundle_path) as f:
            json.dump(bundle, f, indent=2)

        orcacosmo_text = OrcaCosmoOrchestrator(bundle).reconstruct()
        orcacosmo_path = os.path.join(orcacosmo_dir, f"{lookup_id}.orcacosmo")
        with open(orcacosmo_path, "w") as f:
            f.write(orcacosmo_text)

        # ── Build checkpoint ──────────────────────────────────────────────────
        return {
            # BaseStage contract
            "item_id":     item_id,
            "status":      "ok",
            "worker_name": worker_name,
            "worker_pid":  worker_pid,
            "log_details": [
                ("basis",    method_used),
                ("fallback", "yes" if fallback_triggered else "no"),
                ("elapsed",  f"{elapsed_total:.1f}s"),
                ("output",   orcacosmo_path),
            ],
            # Stage payload
            "lookup_id":            lookup_id,
            "inchi_key":            entry["inchi_key"],
            "conf_num":             entry.get("conf_num"),
            "smiles":               entry.get("smiles", ""),
            "energy":               entry.get("energy"),
            "provenance":           entry.get("provenance", {}),
            "optimisation_history": entry.get("optimisation_history", []),
            "orcacosmo_history": {
                "method_used":        method_used,
                "default_basis":      default_basis,
                "fallback_triggered": fallback_triggered,
                "elapsed_default":    elapsed_default,
                "elapsed_fallback":   elapsed_fallback,
                "elapsed_total":      elapsed_total,
                "functional":         functional,
                "solvation":          solvation,
                "orca_version":       orca_version,
                "raw_output_paths": {
                    "log":       os.path.join(raw_dir, f"{lookup_id}.log"),
                    "cpcm":      os.path.join(raw_dir, f"{lookup_id}.cpcm"),
                    "cpcm_corr": os.path.join(raw_dir, f"{lookup_id}.cpcm_corr"),
                },
            },
            "orcacosmo_path": orcacosmo_path,
        }

    except Exception as e:
        return {
            "item_id":     item_id,
            "status":      "failed",
            "error":       str(e),
            "worker_name": worker_name,
            "worker_pid":  worker_pid,
            "log_details": [],
        }


# =============================================================================
# Pure helper functions
# =============================================================================

def _write_orca_input(
    path:           str,
    xyz_name:       str,
    base_name:      str,
    basis:          str,
    functional:     str,
    solvation:      str,
    cpcm_radii:     dict,
    cpcm_cut_area:  float,
    cores_per_item: int,
    maxcore_mb:     int,
    charge:         int = 0,
    multiplicity:   int = 1,
):
    """
    Write a single ORCA input file.

    The method line is assembled from config values:
        ! {solvation} {functional} def2-{basis} SP
    e.g. ! CPCM BP86 def2-TZVP SP
    """
    with open(path, "w") as f:
        f.write(f"%MaxCore {maxcore_mb}\n\n")
        f.write(f'%base "{base_name}"\n\n')

        if cores_per_item > 1:
            f.write(f"%pal nprocs {cores_per_item}\nend\n\n")

        f.write("%cpcm\n")
        for Z, r in cpcm_radii.items():
            f.write(f"radius[{Z}]  {r}\n")
        f.write(f"cut_area {cpcm_cut_area}\n")
        f.write("end\n\n")

        f.write(f"! {solvation} {functional} def2-{basis} SP\n\n")

        f.write("%elprop\n  Polar 1\n  Polaratom 1\nend\n\n")

        f.write(f"* xyzfile {charge} {multiplicity} {xyz_name}\n")


def _run_orca(
    workdir:        str,
    inp_file:       str,
    log_path:       str,
    orca_command:   str,
    cores_per_item: int,
):
    """Run ORCA subprocess, capturing stdout+stderr to log_path."""
    env = os.environ.copy()
    orca_dir = os.path.dirname(orca_command)
    env["LD_LIBRARY_PATH"] = f"{orca_dir}:{env.get('LD_LIBRARY_PATH', '')}"
    env["OMP_NUM_THREADS"]  = str(cores_per_item)

    # Remove inherited MPI environment to avoid slot allocation errors
    for key in list(env.keys()):
        if any(k in key for k in ("OMPI", "MPI", "PRTE")):
            del env[key]

    with open(log_path, "w") as f:
        subprocess.run(
            [orca_command, inp_file],
            cwd    = workdir,
            stdout = f,
            stderr = subprocess.STDOUT,
            check  = False,
            env    = env,
        )


def _validate_cpcm(cpcm: str, cpcm_corr: str, log: str) -> bool:
    """Basic size check — returns True if all three output files look valid."""
    for path, min_bytes in [(cpcm, 100), (cpcm_corr, 100), (log, 1000)]:
        if not os.path.exists(path):
            return False
        if os.path.getsize(path) < min_bytes:
            return False
    return True


def _cleanup_failed_run(workdir: str, base: str):
    """Remove stale ORCA output files before a retry."""
    for ext in (".log", ".cpcm", ".cpcm_corr", ".gbw", ".property.txt"):
        p = os.path.join(workdir, f"{base}{ext}")
        if os.path.exists(p):
            os.remove(p)


def _run_orca_with_basis(
    lookup_id:       str,
    orca_base:       str,
    workdir:         str,
    orca_command:    str,
    cpcm_radii:      dict,
    cpcm_cut_area:   float,
    functional:      str,
    solvation:       str,
    default_basis:   str,
    cores_per_item:  int,
    maxcore_mb:      int,
    orca_inputs_dir: str,
    charge:          int = 0,
    multiplicity:    int = 1,
) -> dict:
    """
    Run ORCA with default_basis.

    orca_base is the short safe name used for ORCA %base and all workdir
    filenames.  lookup_id is used only for the archived .inp copies saved
    to orca_inputs_dir so they remain identifiable after the run.

    Returns dict with:
      method_used        str    basis set that produced valid output
      fallback_triggered bool   always False (fallback handled externally)
      elapsed_default    float  wall time for the attempt
      elapsed_fallback   None
      elapsed_total      float

    Raises RuntimeError on failure.  Fallback basis sets should be
    implemented as a separate pipeline stage with fallback_only=True.
    """
    log       = os.path.join(workdir, f"{orca_base}.log")
    cpcm      = os.path.join(workdir, f"{orca_base}.cpcm")
    cpcm_corr = os.path.join(workdir, f"{orca_base}.cpcm_corr")

    inp = os.path.join(workdir, f"{orca_base}_{default_basis.lower()}.inp")
    _write_orca_input(
        path           = inp,
        xyz_name       = f"{orca_base}.xyz",
        base_name      = orca_base,
        basis          = default_basis,
        functional     = functional,
        solvation      = solvation,
        cpcm_radii     = cpcm_radii,
        cpcm_cut_area  = cpcm_cut_area,
        cores_per_item = cores_per_item,
        maxcore_mb     = maxcore_mb,
        charge         = charge,
        multiplicity   = multiplicity,
    )
    # Archive .inp under lookup_id so it's identifiable after the run
    shutil.copy(
        inp,
        os.path.join(orca_inputs_dir, f"{lookup_id}_{default_basis.lower()}.inp")
    )

    t0 = time.perf_counter()
    _run_orca(workdir, inp, log, orca_command, cores_per_item)
    elapsed = time.perf_counter() - t0

    if not _validate_cpcm(cpcm, cpcm_corr, log):
        raise RuntimeError(
            f"{default_basis} ORCA calculation failed for {lookup_id}. "
            f"To retry with a different basis, add a second orcacosmo stage "
            f"with the desired basis and fallback_only=True."
        )

    return {
        "method_used":        default_basis,
        "fallback_triggered": False,
        "elapsed_default":    elapsed,
        "elapsed_fallback":   None,
        "elapsed_total":      elapsed,
    }


def _parse_output(
    lookup_id:  str,
    raw_dir:    str,
    kind:       str,
    parser_cls,
) -> dict:
    """Run a parser on a raw output file. Raises RuntimeError if missing."""
    path = os.path.join(raw_dir, f"{lookup_id}.{kind}")
    if not os.path.exists(path):
        raise RuntimeError(
            f"Cannot parse {kind} — file not found: {path}"
        )
    return parser_cls(path).parse()


# =============================================================================
# Stage class
# =============================================================================

class OrcacosmoStage(BaseStage):
    """
    ORCA COSMO stage.

    execute() loads config and input data, delegates all per-item work
    to _orcacosmo_worker via BaseStage.run_parallel(), then assembles
    the final summary from checkpoints.

    Stages never call mark_complete() — BaseStage.run() owns lifecycle.
    """

    # =========================================================================
    # Entry point
    # =========================================================================

    def execute(self):
        self.set_stage_output("orcacosmo_summary.json")

        self._load_stage_config()
        self._load_cpcm_data()
        self._prepare_directories()
        self._load_entries()

        all_lookup_ids = self._discover_lookup_ids()

        # fallback_only: pass through conformers that already have a valid
        # .orcacosmo file; only submit those without one to the worker pool.
        if self.fallback_only:
            orcacosmo_dir = os.path.join(self.outputs_dir, "orcacosmo_outputs")
            self._pass_through_checkpoints = []
            to_process_ids = []
            for lid in all_lookup_ids:
                p = os.path.join(orcacosmo_dir, f"{lid}.orcacosmo")
                if os.path.exists(p) and os.path.getsize(p) > 100:
                    entry = self.entry_map.get(lid, {})
                    self._pass_through_checkpoints.append(
                        self._make_passthrough_checkpoint(lid, entry, p)
                    )
                else:
                    to_process_ids.append(lid)
            self.log_info(
                f"fallback_only=True — "
                f"{len(self._pass_through_checkpoints)} pass-through / "
                f"{len(to_process_ids)} to process"
            )
        else:
            self._pass_through_checkpoints = []
            to_process_ids = all_lookup_ids

        lookup_ids = to_process_ids

        # Resume vs fresh start
        if self.job.pending_items:
            self.log_info(
                f"Resuming — {len(self.job.pending_items)} pending / "
                f"{len(lookup_ids) - len(self.job.pending_items)} already complete"
            )
        else:
            self.set_items(lookup_ids, sort_by_complexity=True)
            self._clear_checkpoints()
            self.log_info(
                f"Fresh run — {len(lookup_ids)} conformers to process"
            )

        # Build short ORCA base name mapping.
        # ORCA appends suffixes to %base when constructing output filenames
        # (.cpcm, .cpcm_corr, .gbw, .property.txt, etc.).  Long lookup_ids
        # (e.g. INCHIKEY_c001) push these over ORCA's internal path length
        # limit, causing silent failures or crashes.  We assign a stable short
        # name item{N:04d} for the ORCA-internal base, keyed on sorted order
        # so the mapping is deterministic and survives resume.
        # All external outputs (raw_outputs/, parsed_outputs/, orcacosmo_outputs/,
        # checkpoints/) continue to use lookup_id.
        self._orca_base_map = {
            lid: f"item{idx:04d}"
            for idx, lid in enumerate(sorted(lookup_ids))
        }

        # Skip items with missing XYZ before entering pool
        for lookup_id in list(self.job.pending_items):
            xyz = os.path.join(self.inputs_dir, f"{lookup_id}.xyz")
            if not os.path.exists(xyz):
                self.log_skip(lookup_id, "XYZ file not found")
                self.update_progress(lookup_id, success=False)

        self.run_parallel(_orcacosmo_worker, self._build_args)
        self._assemble_summary()

    # =========================================================================
    # Config
    # =========================================================================

    def _load_stage_config(self):
        cfg = self.config
        if not cfg:
            self.fail("OrcacosmoStage: missing config")

        # ── ORCA executable ───────────────────────────────────────────────────
        orca_cfg = (
            cfg.get("orca_6_1_1")
            or cfg.get("orca")
            or cfg.get("orca_6_0_0")
            or {}
        )
        self.orca_command = orca_cfg.get("executable")
        self.orca_home    = orca_cfg.get("home")
        if not self.orca_command:
            self.fail("Missing ORCA executable in config")

        # ── Constant file dirs ────────────────────────────────────────────────
        const_cfg = cfg.get("constant_files", {})
        self.chemistry_dir = const_cfg.get("chemistry_dir")
        self.metadata_dir  = const_cfg.get("metadata_dir")
        if not self.chemistry_dir:
            self.fail("Missing chemistry_dir in config['constant_files']")
        if not self.metadata_dir:
            self.fail("Missing metadata_dir in config['constant_files']")

        # ── Stage defaults ────────────────────────────────────────────────────
        defaults_path = cfg.get("orcacosmo", {}).get("defaults")
        self._orca_defaults = {}
        if defaults_path and os.path.exists(defaults_path):
            with open(defaults_path) as f:
                self._orca_defaults = json.load(f)
        else:
            self.log_warning(
                "orcacosmo_defaults.json not found — using built-in defaults"
            )

        # ── Basis set ─────────────────────────────────────────────────────────
        # Stage parameter overrides the defaults JSON (allows per-stage basis
        # set in the pipeline spec, e.g. for a TZVPD fallback stage).
        self.default_basis = self.parameters.get(
            "default_basis",
            self._orca_defaults.get("default_basis", "TZVPD")
        )

        valid_bases = {"TZVP", "TZVPD"}
        if self.default_basis not in valid_bases:
            self.fail(
                f"default_basis '{self.default_basis}' is not valid. "
                f"Must be one of: {valid_bases}"
            )

        # ── fallback_only ─────────────────────────────────────────────────────
        # When True, conformers that already have a valid .orcacosmo file are
        # passed through unchanged; only those without one are submitted to
        # the worker pool.  Intended for pipeline stages using a fallback
        # basis set (e.g. TZVP after a TZVPD stage).
        self.fallback_only = bool(self.parameters.get("fallback_only", False))

        # ── Chemistry keywords ────────────────────────────────────────────────
        self.functional = self._orca_defaults.get("functional", "BP86")
        self.solvation  = self._orca_defaults.get("solvation",  "CPCM")

        # ── Memory ───────────────────────────────────────────────────────────
        self.maxcore_mb = int(self._orca_defaults.get("maxcore_mb", 2000))

        # ── Strict mode hierarchy ─────────────────────────────────────────────
        # 1. parameters["strict"]    runtime override
        # 2. defaults["strict"]      stage default
        # 3. False                   fallback
        params_strict = self.parameters.get("strict")
        if params_strict is not None:
            self.strict_mode = bool(params_strict)
        else:
            self.strict_mode = bool(self._orca_defaults.get("strict", False))

        # ── Log resolved config ───────────────────────────────────────────────
        self.log_config(
            f"default_basis={self.default_basis}  "
            f"fallback_only={self.fallback_only}  "
            f"functional={self.functional}  "
            f"solvation={self.solvation}"
        )
        self.log_config(
            f"maxcore_mb={self.maxcore_mb}  "
            f"strict={self.strict_mode}"
        )
        self.log_resources(
            f"cores_per_item={self.cores_per_item}  "
            f"n_workers={self.n_workers}  "
            f"budget={self.budget_cores}"
        )

        # ── ORCA version (best-effort) ────────────────────────────────────────
        self.orca_version = "unknown"
        try:
            out = subprocess.check_output(
                [self.orca_command, "-h"], stderr=subprocess.STDOUT
            )
            self.orca_version = (
                out.decode(errors="ignore").splitlines()[0].strip()
            )
        except Exception:
            pass

        self.log_config(f"orca_version={self.orca_version}")

    # =========================================================================
    # CPCM radii
    # =========================================================================

    def _load_cpcm_data(self):
        cpcm_file = os.path.join(self.chemistry_dir, "cpcm_radii.json")
        if not os.path.exists(cpcm_file):
            self.fail(f"Missing CPCM radii JSON: {cpcm_file}")

        with open(cpcm_file) as f:
            cpcm = json.load(f)

        self.cpcm_radii    = cpcm["radii"]
        self.cpcm_cut_area = cpcm["cut_area"]
        self.cpcm_file     = cpcm_file
        self.log_info(f"Loaded CPCM radii from {cpcm_file}")

    # =========================================================================
    # Directories
    # =========================================================================

    def _prepare_directories(self):
        for subdir in (
            "raw_outputs", "parsed_outputs",
            "orcacosmo_outputs", "orca_inputs",
        ):
            os.makedirs(os.path.join(self.outputs_dir, subdir), exist_ok=True)
        self.log_info("Output directories prepared")

    # =========================================================================
    # Load entries from energies.json
    # =========================================================================

    def _load_entries(self):
        energies_file = self.require_file(
            self.get_stage_input(), "stage_input energies.json"
        )
        self.log_info(f"Stage input: {energies_file}")

        os.makedirs(self.inputs_dir, exist_ok=True)
        shutil.copy(
            energies_file,
            os.path.join(self.inputs_dir, "energies.json")
        )

        with open(energies_file) as f:
            entries = json.load(f)

        self.entry_map = {e["lookup_id"]: e for e in entries}

        copied = 0
        for entry in entries:
            lid = entry.get("lookup_id")
            src = entry.get("xyz_path")
            if lid and src and os.path.exists(src):
                shutil.copy(src, os.path.join(self.inputs_dir, f"{lid}.xyz"))
                copied += 1

        self.log_info(
            f"Loaded {len(entries)} entries from energies.json "
            f"({copied} XYZ files copied)"
        )

    # =========================================================================
    # Discover lookup IDs
    # =========================================================================

    def _discover_lookup_ids(self) -> list:
        ids = sorted(
            os.path.splitext(f)[0]
            for f in os.listdir(self.inputs_dir)
            if f.endswith(".xyz")
        )
        self.log_info(f"Discovered {len(ids)} XYZ structures")
        return ids

    # =========================================================================
    # Build args for worker
    # =========================================================================

    def _build_args(self, item: str) -> dict:
        """
        Build self-contained args dict for _orcacosmo_worker.
        Extends _base_args with all ORCA-specific data.
        Everything here must be JSON-serialisable.
        """
        entry      = self.entry_map.get(item, {})
        provenance = entry.get("provenance", {})
        charge       = int(provenance.get("charge", 0))
        multiplicity = int(provenance.get("multiplicity", 1))
        return {
            **self._base_args(item),
            "lookup_id":       item,
            "orca_base":       self._orca_base_map[item],   # short safe name for ORCA internal use
            "xyz_path":        os.path.join(self.inputs_dir, f"{item}.xyz"),
            "orca_command":    self.orca_command,
            "orca_version":    self.orca_version,
            "cpcm_radii":      self.cpcm_radii,
            "cpcm_cut_area":   self.cpcm_cut_area,
            "cpcm_file":       self.cpcm_file,
            "default_basis":   self.default_basis,
            "functional":      self.functional,
            "solvation":       self.solvation,
            "maxcore_mb":      self.maxcore_mb,
            "entry":           entry,
            "charge":          charge,
            "multiplicity":    multiplicity,
        }

    # =========================================================================
    # Assemble canonical output from checkpoints
    # =========================================================================

    def _make_passthrough_checkpoint(self, lookup_id: str, entry: dict, orcacosmo_path: str) -> dict:
        """Build a minimal checkpoint dict for a pass-through conformer."""
        return {
            "item_id":              lookup_id,
            "status":               "ok",
            "lookup_id":            lookup_id,
            "inchi_key":            entry.get("inchi_key", ""),
            "conf_num":             entry.get("conf_num"),
            "smiles":               entry.get("smiles", ""),
            "energy":               entry.get("energy"),
            "provenance":           entry.get("provenance", {}),
            "optimisation_history": entry.get("optimisation_history", []),
            "orcacosmo_history": {
                "method_used":        "pass-through",
                "default_basis":      self.default_basis,
                "fallback_triggered": False,
                "elapsed_default":    0.0,
                "elapsed_fallback":   None,
                "elapsed_total":      0.0,
                "functional":         self.functional,
                "solvation":          self.solvation,
                "orca_version":       self.orca_version,
                "raw_output_paths":   {},
            },
            "orcacosmo_path": orcacosmo_path,
        }

    def _assemble_summary(self):
        checkpoints = self._load_all_checkpoints()

        # Merge pass-through records (fallback_only mode) before new checkpoints
        if self._pass_through_checkpoints:
            checkpoints = self._pass_through_checkpoints + checkpoints

        if not checkpoints:
            self.log_warning("No successful checkpoints to assemble")
            return

        # ── orcacosmo_summary.json ────────────────────────────────────────────
        results_path = self.get_stage_output()
        with AtomicWriter(results_path) as f:
            json.dump(checkpoints, f, indent=2)
        self.log_info(
            f"Wrote orcacosmo_summary.json — {len(checkpoints)} entries"
        )

        # ── orcacosmo_summary.csv ─────────────────────────────────────────────
        csv_path = os.path.join(self.outputs_dir, "orcacosmo_summary.csv")
        rows = []
        for cp in checkpoints:
            oh = cp.get("orcacosmo_history", {})
            rows.append({
                "lookup_id":          cp.get("lookup_id"),
                "inchi_key":          cp.get("inchi_key"),
                "conf_num":           cp.get("conf_num"),
                "smiles":             cp.get("smiles"),
                "energy":             cp.get("energy"),
                "method_used":        oh.get("method_used"),
                "default_basis":      oh.get("default_basis"),
                "fallback_triggered": oh.get("fallback_triggered"),
                "functional":         oh.get("functional"),
                "solvation":          oh.get("solvation"),
                "elapsed_total":      oh.get("elapsed_total"),
                "elapsed_default":    oh.get("elapsed_default"),
                "elapsed_fallback":   oh.get("elapsed_fallback"),
                "orcacosmo_path":     cp.get("orcacosmo_path"),
            })

        if rows:
            pd.DataFrame(rows).to_csv(csv_path, index=False)
            self.log_info(
                f"Wrote orcacosmo_summary.csv — {len(rows)} rows"
            )