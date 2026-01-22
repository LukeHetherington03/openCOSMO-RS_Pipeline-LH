#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import time
import shutil
import subprocess
import math

import pandas as pd
import numpy as np

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter

from modules.parsers.opt.gxtb_log_parser import GxTBLogParser
from modules.parsers.opt.orca_log_parser import ORCALogParser
from modules.parsers.opt.xtb_log_parser import XTBLogParser


# -------------------------------------------------------------------------
# Engine dispatch table
# -------------------------------------------------------------------------
OPTIMISATION_BACKENDS = {
    "xtb": ("_run_xtb", None),
    "gxtb": ("_run_gxtb", None),
    "orca": ("_run_orca_opt", "auto"),
    "orca_fast": ("_run_orca_opt", "fast"),
    "orca_final": ("_run_orca_opt", "final"),
    "forcefield": ("_run_forcefield", None),
}

# Parser dispatch table
PARSER_BACKENDS = {
    "xtb": XTBLogParser,
    "gxtb": GxTBLogParser,
    "orca": ORCALogParser,
    "orca_fast": ORCALogParser,
    "orca_final": ORCALogParser,
}


class OptimisationRecord:
    def __init__(self, lookup_id, inchi_key, energy, xyz_path, log_path, method, status, metadata=None):
        self.lookup_id = lookup_id
        self.inchikey = inchi_key
        self.energy = energy
        self.xyz_path = xyz_path
        self.log_path = log_path
        self.method = method
        self.status = status
        self.metadata = metadata or {}


    def to_energy_entry(self):
        return {
            "lookup_id": self.lookup_id,
            "inchi_key": self.inchikey,   
            "energy": self.energy,
            "xyz_path": self.xyz_path,
            "log_path": self.log_path,
            "metadata": {
                "method": self.method,
                "status": self.status,
                **self.metadata,
            },
        }



class OptimisationStage(BaseStage):
    """
    Geometry optimisation stage.

    Responsibilities:
      - Load conformers and XYZs from previous stage's energies.json
      - Optimise each conformer with selected engine
      - Track progress per conformer (job items)
      - Enforce energy sanity checks
      - Handle partial/failed optimisations without crashing
      - Produce:
          - outputs/xyz/*.xyz
          - outputs/summary.csv
          - outputs/energies.json
          - outputs/optimisation_summary.csv
          - job_state.json
          - stage_summary.json
    """

    # -------------------------------------------------------------------------
    # Entry point
    # -------------------------------------------------------------------------
    def execute(self):
        stage_start = time.perf_counter()
        cpu_count = os.cpu_count() or 1

        self.log_header("Starting Optimisation Stage")

        params = self.parameters
        engine = params.get("engine", "gxtb").lower()
        max_iter = params.get("max_iter", 250)
        global_fail_threshold = params.get("global_fail_threshold", 0.8)  # fraction

        self.log(f"Engine: {engine}, max_iter={max_iter}, global_fail_threshold={global_fail_threshold}")

        summary_file = params.get("summary_file")
        if summary_file is None:
            self.fail("OptimisationStage requires summary_file (energies.json) from previous stage.")

        prep, entries = self._prepare_inputs(summary_file)

        # Conformer-level items
        conformer_ids = [e["lookup_id"] for e in entries]
        self.set_items(conformer_ids)

        self.log(f"Input: {len(entries)} conformers after filtering missing XYZ")

        # Group by molecule for summary only
        molecules = self._group_by_molecule(entries)

        # Optimise all conformers
        records = self._optimise_all(entries, engine, max_iter, global_fail_threshold)

        # Unified output writer
        self._write_outputs(records)

        # Human-readable molecule-level summary
        self._write_human_summary_csv(records, molecules, engine, self.outputs_dir)

        # ---------------------------------------------------------------------
        # Write job_state.json
        # ---------------------------------------------------------------------
        stage_end = time.perf_counter()
        job_state_path = os.path.join(self.job.job_dir, "job_state.json")

        with AtomicWriter(job_state_path) as f:
            json.dump(
                {
                    "stage": "optimisation",
                    "engine": engine,
                    "max_iter": max_iter,
                    "num_input": len(entries),
                    "num_output": len(records),
                    "missing_xyz": prep["missing"],
                    "cpu_count": cpu_count,
                    "elapsed_seconds": stage_end - stage_start,
                },
                f,
                indent=2,
            )

        # ---------------------------------------------------------------------
        # Write stage summary
        # ---------------------------------------------------------------------
        summary_path = os.path.join(self.outputs_dir, "stage_summary.json")
        with AtomicWriter(summary_path) as f:
            json.dump(
                {
                    "stage": "optimisation",
                    "engine": engine,
                    "total_conformers": len(entries),
                    "optimised_conformers": len(records),
                    "missing_xyz": prep["missing"],
                    "cpu_count": cpu_count,
                    "elapsed_seconds": stage_end - stage_start,
                },
                f,
                indent=2,
            )

        self.job.mark_complete()

    # -------------------------------------------------------------------------
    # Input preparation
    # -------------------------------------------------------------------------
    def _prepare_inputs(self, summary_file):
        inputs_dir = self.inputs_dir
        outputs_dir = self.outputs_dir

        inputs_xyz_dir = os.path.join(inputs_dir, "xyz")
        outputs_xyz_dir = os.path.join(outputs_dir, "xyz")
        outputs_log_dir = os.path.join(outputs_dir, "log")

        os.makedirs(inputs_xyz_dir, exist_ok=True)
        os.makedirs(outputs_xyz_dir, exist_ok=True)
        os.makedirs(outputs_log_dir, exist_ok=True)

        # Copy energies.json locally into inputs/
        local_energies = os.path.join(inputs_dir, "energies.json")
        shutil.copy(summary_file, local_energies)

        with open(local_energies) as f:
            entries = json.load(f)

        if not isinstance(entries, list):
            self.fail("summary_file must contain a list of conformer entries.")

        # Copy XYZs into inputs/xyz, track missing
        missing = []
        valid_entries = []

        for entry in entries:
            lookup_id = entry.get("lookup_id")
            src = entry.get("xyz_path")

            if not lookup_id or not src:
                self.log(f"[WARNING] Entry missing lookup_id or xyz_path; skipping: {entry}")
                missing.append(lookup_id or "<unknown>")
                continue

            dst = os.path.join(inputs_xyz_dir, os.path.basename(src))

            if os.path.exists(src):
                if os.path.abspath(src) != os.path.abspath(dst):
                    shutil.copy(src, dst)
                entry["xyz_path"] = dst
                valid_entries.append(entry)
            else:
                self.log(f"[WARNING] Missing XYZ for {lookup_id}: {src}")
                missing.append(lookup_id)

        if not valid_entries:
            self.fail("No valid conformers available for optimisation (all missing XYZ).")

        return {"inputs_dir": inputs_dir, "outputs_dir": outputs_dir, "missing": missing}, valid_entries

    # -------------------------------------------------------------------------
    # Group by molecule (for summary only)
    # -------------------------------------------------------------------------
    def _group_by_molecule(self, entries):
        groups = {}
        for e in entries:
            mol_id = e.get("inchi_key")  
            if mol_id is None:
                continue
            groups.setdefault(mol_id, []).append(e)
        return groups


    # -------------------------------------------------------------------------
    # Optimise all conformers (conformer-level items)
    # -------------------------------------------------------------------------
    def _optimise_all(self, entries, engine, max_iter, global_fail_threshold):
        records = []
        total_confs = len(entries)
        failures = 0

        for idx, entry in enumerate(entries, start=1):
            lookup_id = entry["lookup_id"]
            self.log(f"[{idx}/{total_confs}] Optimising {lookup_id}", indent=2)

            try:
                record = self._optimise_single(entry, engine, max_iter)
                records.append(record)

                if record.status == "converged":
                    self.update_progress(lookup_id, success=True)
                else:
                    self.update_progress(lookup_id, success=False)
                    failures += 1

            except Exception as e:
                failures += 1
                self.update_progress(lookup_id, success=False)
                self.log(f"[ERROR] Exception during optimisation of {lookup_id}: {e}")

            # Global failure threshold check
            if total_confs > 0 and failures / total_confs >= global_fail_threshold:
                self.log(f"[ERROR] Global failure threshold reached ({failures}/{total_confs}); aborting optimisation.")
                break

        # Filter out completely unusable records
        usable = [r for r in records if self._is_record_usable(r)]

        if not usable:
            self.fail("Optimisation produced zero usable conformers.")

        dropped = len(records) - len(usable)
        if dropped > 0:
            self.log(f"[WARNING] Dropped {dropped} records due to unusable outputs (missing xyz or invalid energy).")

        return usable

    # -------------------------------------------------------------------------
    # Check if a record is usable
    # -------------------------------------------------------------------------
    def _is_record_usable(self, record):
        # Require xyz_path and finite energy if status is converged
        if record.status == "converged":
            if record.xyz_path is None or not os.path.exists(record.xyz_path):
                return False
            if record.energy is None or (
                isinstance(record.energy, float)
                and (math.isnan(record.energy) or abs(record.energy) > 1e6)
            ):
                return False
            return True

        # For failed/partial, we can keep them in energies.json if xyz exists
        if record.xyz_path is None or not os.path.exists(record.xyz_path):
            return False
        return True

    # -------------------------------------------------------------------------
    # Optimise a single conformer
    # -------------------------------------------------------------------------
    def _optimise_single(self, entry, engine, max_iter):
        lookup_id = entry["lookup_id"]
        xyz_path = entry["xyz_path"]
        inchi_key=entry["inchi_key"]

        start = time.perf_counter()
        backend_result = self._run_engine(xyz_path, engine, max_iter)
        elapsed_wall = time.perf_counter() - start

        log_out = backend_result.get("log_out")
        xyz_out = backend_result.get("xyz_out")
        backend_meta = backend_result.get("backend_meta", {})

        # Select parser
        parser_cls = PARSER_BACKENDS.get(engine)
        parsed = {
            "energy": None,
            "iterations": None,
            "converged": None,
            "elapsed_seconds": None,
        }

        if parser_cls and log_out and os.path.exists(log_out):
            try:
                parsed = parser_cls.parse(log_out)
            except Exception as e:
                self.log(f"[WARNING] Parser failed for {lookup_id} ({engine}): {e}")

        energy = parsed.get("energy")
        converged_flag = parsed.get("converged")

        # Option B: parser decides convergence, XYZ decides success
        if converged_flag is True:
            if xyz_out and os.path.exists(xyz_out):
                status = "converged"
            else:
                status = "partial"
        elif converged_flag is False:
            status = "failed"
        else:
            status = "partial"

        # Energy sanity check
        if energy is not None:
            try:
                energy_val = float(energy)
                if math.isnan(energy_val) or abs(energy_val) > 1e6:
                    self.log(f"[WARNING] Discarding insane energy for {lookup_id}: {energy}")
                    energy = None
                    if status == "converged":
                        status = "partial"
            except Exception:
                self.log(f"[WARNING] Failed to parse energy for {lookup_id}: {energy}")
                energy = None
                if status == "converged":
                    status = "partial"

        # Normalised record
        record = OptimisationRecord(
            lookup_id=lookup_id,
            inchi_key=inchi_key,
            energy=energy,
            xyz_path=xyz_out,
            log_path=log_out,
            method=engine,
            status=status,
            metadata={
                "elapsed_seconds": elapsed_wall,
                "parser_iterations": parsed.get("iterations"),
                "parser_converged": parsed.get("converged"),
                "parser_elapsed_seconds": parsed.get("elapsed_seconds"),
                **backend_meta,
            },
        )

        return record

    # -------------------------------------------------------------------------
    # Engine dispatch
    # -------------------------------------------------------------------------
    def _run_engine(self, xyz_path, engine, max_iter):
        if engine not in OPTIMISATION_BACKENDS:
            self.fail(f"Unknown optimisation engine: {engine}")

        method_name, level = OPTIMISATION_BACKENDS[engine]
        runner = getattr(self, method_name)

        if method_name == "_run_orca_opt":
            if level == "auto":
                level = self.parameters.get("level", "final")
            if level not in ("fast", "final"):
                level = "final"
            return runner(xyz_path, max_iter, level)

        return runner(xyz_path, max_iter)

    # -------------------------------------------------------------------------
    # Forcefield placeholder
    # -------------------------------------------------------------------------
    def _run_forcefield(self, xyz_path, max_iter):
        """
        Placeholder forcefield optimisation.

        Currently:
          - Copies input XYZ to output
          - No real optimisation or energy
        """
        job = self.job
        base = os.path.splitext(os.path.basename(xyz_path))[0]

        log_out = os.path.join(job.outputs_dir, "log", f"{base}_ff.log")
        xyz_out = os.path.join(job.outputs_dir, "xyz", f"{base}_ff_opt.xyz")

        os.makedirs(os.path.dirname(log_out), exist_ok=True)
        os.makedirs(os.path.dirname(xyz_out), exist_ok=True)

        with open(log_out, "w") as f:
            f.write("Forcefield optimisation placeholder\n")

        shutil.copy(xyz_path, xyz_out)

        return {
            "log_out": log_out,
            "xyz_out": xyz_out,
            "backend_meta": {
                "mode": "forcefield",
            },
        }

    # -------------------------------------------------------------------------
    # ORCA backend
    # -------------------------------------------------------------------------
    def _run_orca_opt(self, xyz_path, max_iter, level):
        job = self.job
        base = os.path.splitext(os.path.basename(xyz_path))[0]

        # ------------------------------
        # Dedicated scratch directory
        # ------------------------------
        scratch_root = os.path.join(job.job_dir, "orca_opt_scratch")
        scratch_dir = os.path.join(scratch_root, base)
        os.makedirs(scratch_dir, exist_ok=True)

        # Input/output paths
        inp_file = os.path.join(scratch_dir, f"{base}.inp")
        log_out = os.path.join(job.outputs_dir, "log", f"{base}_{level}.log")
        xyz_out = os.path.join(job.outputs_dir, "xyz", f"{base}_{level}_opt.xyz")

        os.makedirs(os.path.dirname(log_out), exist_ok=True)
        os.makedirs(os.path.dirname(xyz_out), exist_ok=True)

        # ------------------------------
        # ORCA input
        # ------------------------------
        if level == "fast":
            method_line = "! BP86 def2-TZVP(-f) Opt"
        else:
            method_line = "! BP86 def2-TZVP Opt"

        with open(inp_file, "w") as f:
            f.write(method_line + "\n")
            f.write("* xyz 0 1\n")
            with open(xyz_path) as xyz:
                lines = xyz.readlines()
                f.writelines(lines[2:])  # skip XYZ header
            f.write("*\n")

        backend_meta = {
            "mode": "orca",
            "level": level,
            "scratch_dir": scratch_dir,
            "run_status": "unknown",
        }

        # ------------------------------
        # Run ORCA inside scratch
        # ------------------------------
        try:
            with open(log_out, "w") as log:
                subprocess.run(
                    ["orca", inp_file],
                    cwd=scratch_dir,      # <-- KEY CHANGE
                    stdout=log,
                    stderr=log,
                    check=True,
                )
            backend_meta["run_status"] = "ok"

            # ORCA writes final geometry as <base>.xyz in cwd
            final_xyz = os.path.join(scratch_dir, f"{base}.xyz")
            if os.path.exists(final_xyz):
                shutil.copy(final_xyz, xyz_out)
            else:
                self.log(f"[WARNING] ORCA did not produce final XYZ for {base}")
                xyz_out = None

        except subprocess.CalledProcessError:
            backend_meta["run_status"] = "failed"
            xyz_out = None

        return {
            "log_out": log_out if os.path.exists(log_out) else None,
            "xyz_out": xyz_out if os.path.exists(xyz_out) else None,
            "backend_meta": backend_meta,
        }

    # -------------------------------------------------------------------------
    # gXTB backend
    # -------------------------------------------------------------------------
    def _run_gxtb(self, xyz_path, max_iter):
        job = self.job
        base = os.path.splitext(os.path.basename(xyz_path))[0]

        tmp_exec = os.path.join(job.job_dir, "gxtb_scratch")
        os.makedirs(tmp_exec, exist_ok=True)

        tmp_xyz = os.path.join(tmp_exec, f"{base}.xyz")
        opt_tmp = os.path.join(tmp_exec, "xtbopt.xyz")

        for fpath in (tmp_xyz, opt_tmp):
            if os.path.exists(fpath):
                os.remove(fpath)

        shutil.copy(xyz_path, tmp_xyz)

        log_out = os.path.join(job.outputs_dir, "log", f"{base}_gxtb.log")
        xyz_out = os.path.join(job.outputs_dir, "xyz", f"{base}_opt.xyz")

        os.makedirs(os.path.dirname(log_out), exist_ok=True)
        os.makedirs(os.path.dirname(xyz_out), exist_ok=True)

        max_cores = os.cpu_count() or 1
        ncores = max(1, math.floor(max_cores * 0.8))

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(ncores)

        cmd = [
            "xtb", tmp_xyz,
            "--driver", "gxtb -grad -c xtbdriver.xyz",
            "--opt",
            "--iterations", str(max_iter),
        ]

        backend_meta = {
            "mode": "gxtb",
            "ncores_used": ncores,
            "scratch_dir": tmp_exec,
            "run_status": "unknown",
        }

        try:
            with open(log_out, "w") as log:
                subprocess.run(
                    cmd,
                    cwd=tmp_exec,
                    stdout=log,
                    stderr=log,
                    check=True,
                    env=env,
                )
            backend_meta["run_status"] = "ok"

            if os.path.exists(opt_tmp):
                shutil.move(opt_tmp, xyz_out)
            else:
                self.log(f"[WARNING] gXTB did not produce xtbopt.xyz for {base}")
                xyz_out = None

        except subprocess.CalledProcessError:
            backend_meta["run_status"] = "failed"

        # Cleanup scratch input
        if os.path.exists(tmp_xyz):
            os.remove(tmp_xyz)

        return {
            "log_out": log_out if os.path.exists(log_out) else None,
            "xyz_out": xyz_out if os.path.exists(xyz_out) else None,
            "backend_meta": backend_meta,
        }

    # -------------------------------------------------------------------------
    # Classic XTB backend
    # -------------------------------------------------------------------------
    def _run_xtb(self, xyz_path, max_iter):
        job = self.job
        base = os.path.splitext(os.path.basename(xyz_path))[0]

        scratch = os.path.join(job.job_dir, "xtb_scratch")
        os.makedirs(scratch, exist_ok=True)

        tmp_xyz = os.path.join(scratch, f"{base}.xyz")
        opt_tmp = os.path.join(scratch, "xtbopt.xyz")

        # Clean scratch
        for fpath in (tmp_xyz, opt_tmp):
            if os.path.exists(fpath):
                os.remove(fpath)

        shutil.copy(xyz_path, tmp_xyz)

        log_out = os.path.join(job.outputs_dir, "log", f"{base}_xtb.log")
        xyz_out = os.path.join(job.outputs_dir, "xyz", f"{base}_opt.xyz")

        os.makedirs(os.path.dirname(log_out), exist_ok=True)
        os.makedirs(os.path.dirname(xyz_out), exist_ok=True)

        cmd = [
            "xtb", tmp_xyz,
            "--opt",
            "--gfn", job.parameters.get("level", "2"),
            "--iterations", str(max_iter),
        ]

        backend_meta = {
            "mode": "xtb",
            "scratch_dir": scratch,
            "run_status": "unknown",
        }

        try:
            with open(log_out, "w") as log:
                subprocess.run(
                    cmd,
                    cwd=scratch,
                    stdout=log,
                    stderr=log,
                    check=True,
                )
            backend_meta["run_status"] = "ok"

            if os.path.exists(opt_tmp):
                shutil.move(opt_tmp, xyz_out)
            else:
                self.log(f"[WARNING] XTB did not produce xtbopt.xyz for {base}")
                xyz_out = None

        except subprocess.CalledProcessError:
            backend_meta["run_status"] = "failed"

        if os.path.exists(tmp_xyz):
            os.remove(tmp_xyz)

        return {
            "log_out": log_out if os.path.exists(log_out) else None,
            "xyz_out": xyz_out if os.path.exists(xyz_out) else None,
            "backend_meta": backend_meta,
        }

    # -------------------------------------------------------------------------
    # Unified output writer
    # -------------------------------------------------------------------------
    def _write_outputs(self, records):
        """
        Produces:
          - outputs/xyz/*.xyz
          - outputs/summary.csv
          - outputs/energies.json
        """
        outputs_dir = self.outputs_dir
        outputs_xyz_dir = os.path.join(outputs_dir, "xyz")
        os.makedirs(outputs_xyz_dir, exist_ok=True)

        # Copy XYZs into outputs/xyz, skip missing paths
        for rec in records:
            src = rec.xyz_path
            if not src or not os.path.exists(src):
                self.log(f"[WARNING] Missing xyz_path for {rec.lookup_id}, skipping copy.")
                rec.xyz_path = None
                continue

            dst = os.path.join(outputs_xyz_dir, os.path.basename(src))
            if os.path.abspath(src) != os.path.abspath(dst):
                try:
                    shutil.copy(src, dst)
                except Exception as e:
                    self.log(f"[ERROR] Failed to copy {src} → {dst}: {e}")
                    rec.xyz_path = None
                    continue

            rec.xyz_path = dst

        # summary.csv (human-readable per conformer)
        summary_csv = os.path.join(outputs_dir, "summary.csv")
        df = pd.DataFrame(
            [
                {
                    "lookup_id": r.lookup_id,
                    "energy": r.energy,
                    "xyz_path": r.xyz_path,
                    "log_path": r.log_path,
                    "status": r.status,
                    "method": r.method,
                    "elapsed_seconds": r.metadata.get("elapsed_seconds"),
                }
                for r in records
            ]
        )
        df.to_csv(summary_csv, index=False)

        # energies.json (canonical)
        energies_json = os.path.join(outputs_dir, "energies.json")
        with AtomicWriter(energies_json) as f:
            json.dump([r.to_energy_entry() for r in records], f, indent=2)

        self.log(f"Optimisation outputs written to: {outputs_dir}")

    def _write_human_summary_csv(self, records, molecules, engine, outputs_dir):
        summary_csv = os.path.join(outputs_dir, "optimisation_summary.csv")
        mol_stats = []

        # Pre-index records by molecule ID (inchi_key)
        by_mol = {}
        for r in records:
            mol_id = r.inchikey   # ← authoritative
            by_mol.setdefault(mol_id, []).append(r)

        for mol_id, mol_entries in molecules.items():
            mol_records = by_mol.get(mol_id, [])

            n_confs = len(mol_records)
            n_success = sum(1 for r in mol_records if r.status == "converged")
            n_failed = sum(1 for r in mol_records if r.status == "failed")
            n_partial = sum(1 for r in mol_records if r.status == "partial")

            times = [r.metadata.get("elapsed_seconds", 0.0) for r in mol_records]
            avg_time = float(np.mean(times)) if times else 0.0
            total_time = float(np.sum(times)) if times else 0.0

            # Estimate n_atoms from first valid xyz
            n_atoms = None
            for r in mol_records:
                if r.xyz_path and os.path.exists(r.xyz_path):
                    try:
                        with open(r.xyz_path) as f:
                            first_line = f.readline().strip()
                            n_atoms = int(first_line)
                        break
                    except Exception:
                        n_atoms = None

            mol_stats.append(
                {
                    "molecule_id": mol_id,
                    "n_atoms": n_atoms,
                    "n_conformers_attempted": len(mol_entries),
                    "n_conformers_output": n_confs,
                    "n_converged": n_success,
                    "n_failed": n_failed,
                    "n_partial": n_partial,
                    "avg_time_per_conf_s": avg_time,
                    "total_time_s": total_time,
                    "engine": engine,
                }
            )

        pd.DataFrame(mol_stats).to_csv(summary_csv, index=False)
        self.log(f"Human-readable optimisation summary written to: {summary_csv}")
