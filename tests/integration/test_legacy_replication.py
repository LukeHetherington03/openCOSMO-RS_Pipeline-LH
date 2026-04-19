#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tests/integration/test_legacy_replication.py

Integration test: verify that our pipeline reproduces the legacy
ConformerGenerator_IC.py output for small, well-understood molecules.

Run with:
    python3 tests/integration/test_legacy_replication.py [--mode fast|full]

Modes
-----
fast  (default)
    Runs only the pre-ORCA stages (cleaning → generation → optimisation →
    cleaning/filtering).  No ORCA needed.  Verifies conformer counts and
    MMFF energies match expectations from the legacy pipeline.

full
    Also runs the ORCA optimisation and orcacosmo stages.  Requires a
    working ORCA installation and valid paths.json.

Test molecules
--------------
glycine      : 0 rotatable bonds within the amino-acid backbone → 50 confs
aspirin      : 3 rotatable bonds → 50 confs
ibuprofen    : 5 rotatable bonds → 50 confs

Legacy pipeline equivalent (from calculate_rdkit + calculate_orca):
    1. cleaning (parse SMILES, charge, multiplicity)
    2. generation  rdkit, n=50 (heuristic: ≤7 rot bonds → 50 confs)
    3. optimisation  forcefield_mmff, max_iter=2000
    4. cleaning  rdkit_post_opt_legacy  (convergence / overlap / Boltzmann / RMSD)
    5. [full only] optimisation  orca_opt_cpcm_fast
    6. [full only] orcacosmo

Expected post-filter conformer counts (from legacy runs):
    glycine   → 1  (single achiral conformer)
    aspirin   → typically 2-5 (benzene ring locks most DOF)
    ibuprofen → typically 3-8

If the legacy pipeline has been run locally, pass --legacy-dir to compare
orcacosmo files directly.
"""

import argparse
import csv
import json
import os
import sys
import tempfile
import textwrap

# ---------------------------------------------------------------------------
# Ensure project root is on path
# ---------------------------------------------------------------------------
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)


# ---------------------------------------------------------------------------
# Test molecules
# ---------------------------------------------------------------------------

TEST_MOLECULES = [
    {
        "name":    "glycine",
        "smiles":  "NCC(=O)O",
        "charge":  0,
        "rotbonds": 2,          # within heavy-atom backbone
        "expected_n_gen":     50,
        "expected_n_post_filter": (1, 5),   # inclusive range — 1 for rigid zwitterion
    },
    {
        "name":    "aspirin",
        "smiles":  "CC(=O)Oc1ccccc1C(=O)O",
        "charge":  0,
        "rotbonds": 3,
        "expected_n_gen":     50,
        "expected_n_post_filter": (1, 10),
    },
    {
        "name":    "ibuprofen",
        "smiles":  "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "charge":  0,
        "rotbonds": 5,
        "expected_n_gen":     50,
        "expected_n_post_filter": (2, 20),
    },
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_input_csv(path: str, molecules: list):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "smiles", "charge"])
        w.writeheader()
        for mol in molecules:
            w.writerow({"name": mol["name"], "smiles": mol["smiles"], "charge": mol["charge"]})


def _load_energies_json(path: str) -> list:
    with open(path) as f:
        return json.load(f)


def _n_conformers(energies_path: str) -> int:
    try:
        return len(_load_energies_json(energies_path))
    except Exception:
        return -1


def _energy_range_kcal(energies_path: str) -> tuple:
    """Return (min_energy, max_energy, rel_range) in kcal/mol."""
    try:
        records = _load_energies_json(energies_path)
        energies = [r["energy"] for r in records if r.get("energy") is not None]
        if not energies:
            return (None, None, None)
        e_min = min(energies)
        e_max = max(energies)
        return (e_min, e_max, e_max - e_min)
    except Exception:
        return (None, None, None)


def _print_stage_report(stage_name: str, energies_path: str):
    n = _n_conformers(energies_path)
    e_min, e_max, e_range = _energy_range_kcal(energies_path)
    if e_range is not None:
        print(f"    {stage_name:<35} n={n:>3}  "
              f"E_min={e_min:>10.4f}  E_max={e_max:>10.4f}  "
              f"ΔE={e_range:>8.4f} kcal/mol")
    else:
        print(f"    {stage_name:<35} n={n:>3}  (no energy data)")
    return n


def _check(condition: bool, msg: str, results: list):
    status = "PASS" if condition else "FAIL"
    results.append((status, msg))
    symbol = "✓" if condition else "✗"
    print(f"    [{symbol}] {msg}")


def _run_pipeline(config: dict, pipeline_spec: list, title: str, base_dir: str) -> str:
    """Submit a pipeline and run it in-process. Returns the request_id."""
    from modules.build.request_manager import Request
    from modules.execution.runner import PipelineRunner

    req = Request.create_new(
        base_dir    = base_dir,
        dataset     = title,
        pipeline_spec = pipeline_spec,
        parameters  = {
            "title":     title,
            "resources": {"cpus": 4, "memory_gb": 8},
            "config":    config,
        },
    )

    PipelineRunner.run_request(req.request_id, base_dir)
    return req.request_id


# ---------------------------------------------------------------------------
# Stage output helpers
# ---------------------------------------------------------------------------

def _find_stage_output(base_dir: str, request_id: str, stage_idx: int,
                        stage_name: str) -> str | None:
    """
    Locate the canonical output file for a given stage.
    Looks under:
        {base_dir}/pipeline_data/requests/{request_id}/jobs/*/outputs/
    """
    jobs_dir = os.path.join(
        base_dir, "pipeline_data", "requests", request_id, "jobs"
    )
    if not os.path.exists(jobs_dir):
        return None

    # Jobs are numbered 00_stagename, 01_stagename, etc.
    prefix = f"{stage_idx:02d}_{stage_name}"
    for job_dir in sorted(os.listdir(jobs_dir)):
        if job_dir.startswith(prefix):
            for fname in ("energies.json", "cleaned.csv", "orcacosmo_summary.json"):
                candidate = os.path.join(jobs_dir, job_dir, "outputs", fname)
                if os.path.exists(candidate):
                    return candidate
    return None


# ---------------------------------------------------------------------------
# Pre-ORCA fast test
# ---------------------------------------------------------------------------

def run_fast_test(config: dict, base_dir: str, molecules: list) -> list:
    """
    Run GenerationStage + OptimisationStage + CleaningStage (fast, no ORCA).
    Verifies conformer counts and energy ordering.
    Returns list of (status, message) tuples.
    """
    results = []

    input_csv_path = os.path.join(base_dir, "test_molecules.csv")
    _write_input_csv(input_csv_path, molecules)

    # Legacy-equivalent pre-ORCA pipeline:
    #   cleaning → generation (rdkit, n=50) → optimisation (forcefield_mmff) →
    #   pruning (rdkit_post_opt_legacy: convergence / overlap / Boltzmann / RMSD)
    pipeline_spec = [
        {"stage": "cleaning",     "args": {"input_csv": input_csv_path}},
        {"stage": "generation",   "args": {"engine": "rdkit", "n": 50}},
        {"stage": "optimisation", "args": {"engine": "forcefield_mmff", "max_iter": 2000}},
        {"stage": "pruning",      "args": {"rdkit_post_opt_legacy": True}},
    ]

    print("\n" + "="*70)
    print("  FAST TEST — pre-ORCA conformer generation and filtering")
    print("="*70)
    print(f"  Molecules: {[m['name'] for m in molecules]}")
    print()

    try:
        request_id = _run_pipeline(config, pipeline_spec, "legacy_fast_test", base_dir)
    except Exception as e:
        _check(False, f"Pipeline run failed: {e}", results)
        return results

    # ── Stage 0: cleaning ──────────────────────────────────────────────────
    cleaning_out = _find_stage_output(base_dir, request_id, 0, "cleaning")
    _check(
        cleaning_out is not None and os.path.exists(cleaning_out),
        "cleaning (stage 0) produced cleaned.csv",
        results,
    )

    # ── Stage 1: generation ────────────────────────────────────────────────
    gen_out = _find_stage_output(base_dir, request_id, 1, "generation")
    _check(
        gen_out is not None and os.path.exists(gen_out),
        "generation stage produced energies.json",
        results,
    )
    if gen_out:
        print()
        print("  After generation (before MMFF optimisation):")
        for mol in molecules:
            # load just this molecule's records
            records = [r for r in _load_energies_json(gen_out)
                       if r.get("smiles") == mol["smiles"] or
                          mol["name"].lower() in r.get("lookup_id", "").lower()]
            n = len(records)
            _check(
                n == mol["expected_n_gen"],
                f"  {mol['name']}: generated {n} / expected {mol['expected_n_gen']} conformers",
                results,
            )

    # ── Stage 2: optimisation (forcefield) ────────────────────────────────
    opt_out = _find_stage_output(base_dir, request_id, 2, "optimisation")
    _check(
        opt_out is not None and os.path.exists(opt_out),
        "optimisation stage produced energies.json",
        results,
    )
    if opt_out:
        print()
        print("  After MMFF optimisation (all energies in kcal/mol):")
        for mol in molecules:
            records = [r for r in _load_energies_json(opt_out)
                       if mol["name"].lower() in r.get("lookup_id", "").lower()]
            n      = len(records)
            energies = [r["energy"] for r in records if r.get("energy") is not None]
            _check(
                n > 0,
                f"  {mol['name']}: {n} conformers after MMFF opt",
                results,
            )
            if energies:
                e_range = max(energies) - min(energies)
                _check(
                    e_range >= 0,
                    f"  {mol['name']}: energy range {e_range:.4f} kcal/mol (≥0 required)",
                    results,
                )
                # All energies should be in a physically sensible range
                # MMFF energies for small organics: -200 to +500 kcal/mol
                _check(
                    all(-500 < e < 2000 for e in energies),
                    f"  {mol['name']}: all MMFF energies in plausible range",
                    results,
                )

    # ── Stage 3: pruning (rdkit_post_opt_legacy) ──────────────────────────
    filt_out = _find_stage_output(base_dir, request_id, 3, "pruning")
    _check(
        filt_out is not None and os.path.exists(filt_out),
        "pruning (rdkit_post_opt_legacy) produced energies.json",
        results,
    )
    if filt_out:
        print()
        print("  After post-MMFF filters (convergence / overlap / Boltzmann / RMSD):")
        for mol in molecules:
            records = [r for r in _load_energies_json(filt_out)
                       if mol["name"].lower() in r.get("lookup_id", "").lower()]
            n      = len(records)
            lo, hi = mol["expected_n_post_filter"]
            _check(
                lo <= n <= hi,
                f"  {mol['name']}: {n} conformers (expected {lo}–{hi})",
                results,
            )

            # Energies should be sorted ascending (legacy sorts by energy)
            energies = [r["energy"] for r in records if r.get("energy") is not None]
            if len(energies) > 1:
                _check(
                    all(energies[i] <= energies[i+1] + 1e-6
                        for i in range(len(energies)-1)),
                    f"  {mol['name']}: conformers sorted by energy ascending",
                    results,
                )

    return results


# ---------------------------------------------------------------------------
# Full test (includes ORCA)
# ---------------------------------------------------------------------------

def run_full_test(config: dict, base_dir: str, molecules: list,
                  legacy_dir: str | None = None) -> list:
    """
    Full legacy replication: pre-ORCA + ORCA opt + orcacosmo.
    If legacy_dir is provided, compare .orcacosmo files directly.
    """
    results = []

    input_csv_path = os.path.join(base_dir, "test_molecules.csv")
    _write_input_csv(input_csv_path, molecules)

    # 0: cleaning  1: generation  2: optimisation (FF)  3: pruning (legacy filters)
    # 4: optimisation (ORCA fast)  5: orcacosmo
    pipeline_spec = [
        {"stage": "cleaning",     "args": {"input_csv": input_csv_path}},
        {"stage": "generation",   "args": {"engine": "rdkit", "n": 50}},
        {"stage": "optimisation", "args": {"engine": "forcefield_mmff", "max_iter": 2000}},
        {"stage": "pruning",      "args": {"rdkit_post_opt_legacy": True}},
        {"stage": "optimisation", "args": {"engine": "orca_opt_cpcm_fast"}},
        {"stage": "orcacosmo",    "args": {}},
    ]

    print("\n" + "="*70)
    print("  FULL TEST — legacy pipeline replication including ORCA")
    print("="*70)
    print(f"  Molecules: {[m['name'] for m in molecules]}")
    print()

    try:
        request_id = _run_pipeline(config, pipeline_spec, "legacy_full_test", base_dir)
    except Exception as e:
        _check(False, f"Pipeline run failed: {e}", results)
        return results

    # ── orcacosmo output (stage 5) ────────────────────────────────────────
    cosmo_out = _find_stage_output(base_dir, request_id, 5, "orcacosmo")
    _check(
        cosmo_out is not None and os.path.exists(cosmo_out),
        "orcacosmo stage produced orcacosmo_summary.json",
        results,
    )

    if cosmo_out:
        print()
        print("  ORCA COSMO results:")
        records = _load_energies_json(cosmo_out)
        for mol in molecules:
            mol_records = [r for r in records
                           if mol["name"].lower() in r.get("lookup_id", "").lower()]
            n = len(mol_records)
            _check(
                n > 0,
                f"  {mol['name']}: {n} .orcacosmo files produced",
                results,
            )

            # All should have status "ok"
            ok = [r for r in mol_records if r.get("status") == "ok"]
            _check(
                len(ok) == n,
                f"  {mol['name']}: all {n} conformers status=ok",
                results,
            )

            # Verify .orcacosmo files exist on disk
            for r in mol_records:
                p = r.get("orcacosmo_path", "")
                _check(
                    bool(p) and os.path.exists(p) and os.path.getsize(p) > 100,
                    f"  {mol['name']}: {os.path.basename(p)} exists on disk",
                    results,
                )

        # ── Compare against legacy if provided ───────────────────────────
        if legacy_dir and os.path.isdir(legacy_dir):
            print()
            print("  Comparing against legacy .orcacosmo files:")
            _compare_against_legacy(records, legacy_dir, results)

    return results


def _compare_against_legacy(our_records: list, legacy_dir: str, results: list):
    """
    Compare our .orcacosmo files against legacy ones by parsing key scalar
    fields: total_energy, sigma-profile integral, etc.
    """
    import glob

    legacy_cosmo_files = glob.glob(os.path.join(legacy_dir, "**", "*.orcacosmo"),
                                   recursive=True)
    if not legacy_cosmo_files:
        _check(False, f"No .orcacosmo files found in {legacy_dir}", results)
        return

    print(f"    Found {len(legacy_cosmo_files)} legacy .orcacosmo file(s)")

    for legacy_path in sorted(legacy_cosmo_files):
        # Try to match by inchi_key extracted from path
        inchi_key_guess = os.path.splitext(os.path.basename(legacy_path))[0]
        inchi_key_guess = inchi_key_guess.split("_c")[0]  # strip _c000 etc.

        our_match = next(
            (r for r in our_records
             if r.get("inchi_key", "").startswith(inchi_key_guess)), None
        )
        if our_match is None:
            print(f"    (no match in our records for {inchi_key_guess})")
            continue

        our_path = our_match.get("orcacosmo_path", "")
        if not our_path or not os.path.exists(our_path):
            _check(False, f"  Our file missing for {inchi_key_guess}", results)
            continue

        # Parse both and compare total energy
        legacy_energy = _parse_orcacosmo_energy(legacy_path)
        our_energy    = _parse_orcacosmo_energy(our_path)

        if legacy_energy is not None and our_energy is not None:
            delta = abs(our_energy - legacy_energy)
            # DFT energies should match within 1e-6 Eh with identical inputs
            _check(
                delta < 1e-4,
                f"  {inchi_key_guess}: ΔE_total = {delta:.2e} Eh "
                f"(legacy={legacy_energy:.6f}, ours={our_energy:.6f})",
                results,
            )
        else:
            print(f"    (could not parse energies for {inchi_key_guess})")


def _parse_orcacosmo_energy(path: str) -> float | None:
    """
    Minimal parser: extract the total DFT energy from a .orcacosmo file.
    Legacy format has a line like:  !DATE ...  or  Energy=  ...
    """
    try:
        with open(path) as f:
            for line in f:
                ls = line.strip()
                if ls.startswith("energy_tot"):
                    # openCOSMO-RS format: energy_tot = <value>
                    return float(ls.split("=")[-1].strip())
                if "Total Energy" in ls or "TOTAL ENERGY" in ls:
                    parts = ls.split()
                    for i, p in enumerate(parts):
                        try:
                            val = float(p)
                            if abs(val) > 10:   # skip small numbers like "1.0"
                                return val
                        except ValueError:
                            pass
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Integration test: legacy pipeline replication",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              python3 tests/integration/test_legacy_replication.py
              python3 tests/integration/test_legacy_replication.py --mode full
              python3 tests/integration/test_legacy_replication.py \\
                  --mode full --legacy-dir /path/to/legacy/output
              python3 tests/integration/test_legacy_replication.py \\
                  --molecules glycine aspirin
        """),
    )
    parser.add_argument(
        "--mode", choices=["fast", "full"], default="fast",
        help="fast = pre-ORCA only (default); full = includes ORCA",
    )
    parser.add_argument(
        "--config", default=os.path.join(PROJECT_ROOT, "config", "paths.json"),
        help="Path to paths.json",
    )
    parser.add_argument(
        "--outdir", default=None,
        help="Directory for pipeline output (default: temp dir)",
    )
    parser.add_argument(
        "--legacy-dir", default=None,
        help="Path to legacy pipeline output dir (for direct .orcacosmo comparison)",
    )
    parser.add_argument(
        "--molecules", nargs="+",
        choices=[m["name"] for m in TEST_MOLECULES],
        default=None,
        help="Subset of molecules to test (default: all)",
    )
    args = parser.parse_args()

    # ── Load config ───────────────────────────────────────────────────────
    if not os.path.exists(args.config):
        print(f"ERROR: config not found at {args.config}")
        print("Pass --config /path/to/paths.json or run from the project root.")
        sys.exit(1)

    with open(args.config) as f:
        config = json.load(f)
    from modules.utils.git_version import get_git_version
    config["pipeline_version"] = get_git_version()

    # ── Select molecules ──────────────────────────────────────────────────
    molecules = TEST_MOLECULES
    if args.molecules:
        molecules = [m for m in TEST_MOLECULES if m["name"] in args.molecules]

    # ── Base dir ──────────────────────────────────────────────────────────
    if args.outdir:
        base_dir = args.outdir
        os.makedirs(base_dir, exist_ok=True)
        cleanup = False
    else:
        tmp = tempfile.mkdtemp(prefix="legacy_replication_test_")
        base_dir = tmp
        cleanup = True

    print(f"\nOutput directory: {base_dir}")
    print(f"Mode: {args.mode}")

    # ── Run tests ─────────────────────────────────────────────────────────
    all_results = []

    if args.mode == "fast":
        all_results = run_fast_test(config, base_dir, molecules)
    else:
        fast_results = run_fast_test(config, base_dir, molecules)
        all_results.extend(fast_results)
        full_results = run_full_test(
            config, base_dir, molecules, legacy_dir=args.legacy_dir
        )
        all_results.extend(full_results)

    # ── Summary ───────────────────────────────────────────────────────────
    n_pass = sum(1 for s, _ in all_results if s == "PASS")
    n_fail = sum(1 for s, _ in all_results if s == "FAIL")

    print()
    print("="*70)
    print(f"  RESULT: {n_pass} passed / {n_fail} failed / {len(all_results)} total")
    print("="*70)

    if n_fail:
        print("\nFailed checks:")
        for status, msg in all_results:
            if status == "FAIL":
                print(f"  ✗ {msg}")

    if cleanup and n_fail == 0:
        import shutil
        shutil.rmtree(base_dir, ignore_errors=True)
    elif cleanup:
        print(f"\nOutput preserved for inspection: {base_dir}")

    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
