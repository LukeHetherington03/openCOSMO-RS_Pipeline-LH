#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import time
import shutil
import subprocess
import math
from datetime import datetime

import pandas as pd
import numpy as np

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter
from modules.utils.conformers import ConformerRecord, ConformerSet

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


# -------------------------------------------------------------------------
# Normalised optimisation levels
# -------------------------------------------------------------------------
VALID_LEVELS = ("loose", "normal", "tight", "vtight")


class OptimisationStage(BaseStage):
    """
    Geometry optimisation stage.

    Responsibilities:
      - Load conformers from energies.json (ConformerSet)
      - Copy XYZs into local scratch
      - Optimise each conformer with selected engine
      - Parse logs, extract energies, convergence, timings
      - Update ConformerRecord.provenance["optimisation"]
      - Write:
          - outputs/xyz/*.xyz
          - outputs/summary.csv (per conformer)
          - outputs/energies.json (canonical)
          - outputs/optimisation_summary.csv (per molecule)
          - job_state.json
    """

    # -------------------------------------------------------------------------
    # Entry point
    # -------------------------------------------------------------------------
    def execute(self):
        stage_start = time.perf_counter()
        cpu_count = os.cpu_count() or 1

        self.log_header("Starting Optimisation Stage")

        cfg = self.config or {}
        self.strict_mode = bool(cfg.get("optimisation", {}).get("strict", False))

        params = self.parameters
        engine = params.get("engine", "gxtb").lower()
        level = params.get("level", "normal").lower()
        max_iter = params.get("max_iter", 250)
        global_fail_threshold = params.get("global_fail_threshold", 0.8)

        if level not in VALID_LEVELS:
            self.fail(f"Unknown optimisation level: {level}")

        self.log(f"Engine: {engine}, level={level}, max_iter={max_iter}, strict={self.strict_mode}")
        self.log(f"Global fail threshold: {global_fail_threshold}")

        energies_file = params.get("energies_file")
        if energies_file is None:
            energies_file = self._auto_detect_energies()

        if not os.path.isfile(energies_file):
            self.fail(f"OptimisationStage: energies_file does not exist: {energies_file}")

        # Load conformers
        conformer_set = ConformerSet.load(energies_file)
        entries = conformer_set.records

        # Prepare XYZs
        prep, valid_records = self._prepare_inputs(entries)

        # Conformer-level job items
        self.set_items([rec.lookup_id for rec in valid_records])
        self.log(f"Input: {len(valid_records)} conformers with valid XYZ")

        # Group by molecule for summary
        molecules = self._group_by_molecule(valid_records)

        # Optimise
        optimised_records = self._optimise_all(
            valid_records,
            engine,
            level,
            max_iter,
            global_fail_threshold
        )

        # Write outputs
        self._write_outputs(optimised_records)

        # Human-readable molecule summary
        self._write_human_summary_csv(
            optimised_records,
            molecules,
            engine,
            level,
            self.outputs_dir
        )

        # job_state.json
        stage_end = time.perf_counter()
        job_state_path = os.path.join(self.job.job_dir, "job_state.json")

        with AtomicWriter(job_state_path) as f:
            json.dump(
                {
                    "stage": "optimisation",
                    "engine": engine,
                    "level": level,
                    "max_iter": max_iter,
                    "num_input": len(entries),
                    "num_output": len(optimised_records),
                    "missing_xyz": prep["missing"],
                    "cpu_count": cpu_count,
                    "elapsed_seconds": stage_end - stage_start,
                },
                f,
                indent=2,
            )

        self._log_warning_summary()

        self.job.mark_complete()

    # -------------------------------------------------------------------------
    # Input preparation
    # -------------------------------------------------------------------------
    def _auto_detect_energies(self):
        candidates = [
            os.path.join(self.inputs_dir, "energies.json"),
            os.path.join(self.outputs_dir, "energies.json"),
        ]
        for path in candidates:
            if os.path.isfile(path):
                self.log(f"Auto-detected energies.json at: {path}")
                return path
        self.fail(
            "OptimisationStage: energies_file not provided and energies.json "
            "could not be auto-detected."
        )

    def _prepare_inputs(self, records):
        """
        Copies XYZs into inputs/xyz and returns:
            { "inputs_dir": ..., "outputs_dir": ..., "missing": [...] },
            valid_records
        """
        inputs_dir = self.inputs_dir
        outputs_dir = self.outputs_dir

        inputs_xyz_dir = os.path.join(inputs_dir, "xyz")
        outputs_xyz_dir = os.path.join(outputs_dir, "xyz")
        outputs_log_dir = os.path.join(outputs_dir, "log")

        os.makedirs(inputs_xyz_dir, exist_ok=True)
        os.makedirs(outputs_xyz_dir, exist_ok=True)
        os.makedirs(outputs_log_dir, exist_ok=True)

        missing = []
        valid_records = []

        for rec in records:
            src = rec.xyz_path
            if not src or not os.path.exists(src):
                self.log(f"[WARNING] Missing XYZ for {rec.lookup_id}: {src}")
                missing.append(rec.lookup_id)
                continue

            dst = os.path.join(inputs_xyz_dir, os.path.basename(src))
            if os.path.abspath(src) != os.path.abspath(dst):
                shutil.copy(src, dst)

            rec.xyz_path = dst
            valid_records.append(rec)

        if not valid_records:
            self.fail("No valid conformers available for optimisation (all missing XYZ).")

        return {
            "inputs_dir": inputs_dir,
            "outputs_dir": outputs_dir,
            "missing": missing,
        }, valid_records

    # -------------------------------------------------------------------------
    # Group by molecule
    # -------------------------------------------------------------------------
    def _group_by_molecule(self, records):
        groups = {}
        for rec in records:
            groups.setdefault(rec.inchi_key, []).append(rec)
        return groups

    # -------------------------------------------------------------------------
    # Optimise all conformers
    # -------------------------------------------------------------------------
    def _optimise_all(self, records, engine, level, max_iter, global_fail_threshold):
        optimised = []
        total = len(records)
        failures = 0

        for idx, rec in enumerate(records, start=1):
            self.log(f"[{idx}/{total}] Optimising {rec.lookup_id}", indent=2)

            try:
                updated = self._optimise_single(rec, engine, level, max_iter)
                optimised.append(updated)

                status = updated.provenance["optimisation"]["status"]
                if status == "converged":
                    self.update_progress(rec.lookup_id, success=True)
                else:
                    self.update_progress(rec.lookup_id, success=False)
                    failures += 1

            except Exception as e:
                failures += 1
                self.update_progress(rec.lookup_id, success=False)
                self.log(f"[ERROR] Exception during optimisation of {rec.lookup_id}: {e}")

            if total > 0 and failures / total >= global_fail_threshold:
                self.log(f"[ERROR] Global failure threshold reached ({failures}/{total}); aborting optimisation.")
                break

        usable = [r for r in optimised if self._is_record_usable(r)]

        if not usable:
            self.fail("Optimisation produced zero usable conformers.")

        dropped = len(optimised) - len(usable)
        if dropped > 0:
            self.log(f"[WARNING] Dropped {dropped} unusable conformers (missing XYZ or invalid energy).")

        return usable

    # -------------------------------------------------------------------------
    # Usability check
    # -------------------------------------------------------------------------
    def _is_record_usable(self, rec):
        opt = rec.provenance.get("optimisation", {})
        status = opt.get("status")
        xyz_path = opt.get("xyz_path")
        energy = opt.get("energy")

        if status == "converged":
            if not xyz_path or not os.path.exists(xyz_path):
                return False
            if energy is None or math.isnan(energy) or abs(energy) > 1e6:
                return False
            return True

        if not xyz_path or not os.path.exists(xyz_path):
            return False

        return True

    # -------------------------------------------------------------------------
    # Optimise a single conformer
    # -------------------------------------------------------------------------
    def _optimise_single(self, rec, engine, level, max_iter):
        lookup_id = rec.lookup_id
        xyz_path = rec.xyz_path

        start = time.perf_counter()
        backend_result = self._run_engine(xyz_path, engine, level, max_iter)
        elapsed_wall = time.perf_counter() - start

        log_out = backend_result.get("log_out")
        xyz_out = backend_result.get("xyz_out")
        backend_meta = backend_result.get("backend_meta", {})
        engine_command = backend_result.get("engine_command", "unknown")
        backend_version = backend_result.get("backend_version", "unknown")

        parser_cls = PARSER_BACKENDS.get(engine)
        parsed = {
            "energy": None,
            "iterations": None,
            "converged": None,
            "elapsed_seconds": None,
        }

        warnings_list = []

        if parser_cls and log_out and os.path.exists(log_out):
            try:
                parsed = parser_cls.parse(log_out)
            except Exception as e:
                warnings_list.append(f"Parser failed: {e}")
                self.log(f"[WARNING] Parser failed for {lookup_id} ({engine}): {e}")

        energy = parsed.get("energy")
        converged_flag = parsed.get("converged")

        if converged_flag is True:
            status = "converged" if xyz_out and os.path.exists(xyz_out) else "partial"
        elif converged_flag is False:
            status = "failed"
        else:
            status = "partial"

        if energy is not None:
            try:
                val = float(energy)
                if math.isnan(val) or abs(val) > 1e6:
                    warnings_list.append("Energy insane, discarded")
                    self.log(f"[WARNING] Discarding insane energy for {lookup_id}: {energy}")
                    energy = None
                    if status == "converged":
                        status = "partial"
            except Exception:
                warnings_list.append("Energy parse failed")
                self.log(f"[WARNING] Failed to parse energy for {lookup_id}: {energy}")
                energy = None
                if status == "converged":
                    status = "partial"

        # Convergence quality
        if status == "converged":
            if level == "loose":
                quality = "medium"
            elif level == "normal":
                quality = "good"
            elif level == "tight":
                quality = "very_good"
            else:
                quality = "excellent"
        elif status == "partial":
            quality = "poor"
        else:
            quality = "unusable"

        # Success score
        if status == "converged":
            if level == "vtight":
                score = 1.0
            elif level == "tight":
                score = 0.8
            elif level == "normal":
                score = 0.6
            else:
                score = 0.4
        elif status == "partial":
            score = 0.2
        else:
            score = 0.0

        # Geometry version
        prev_version = rec.provenance.get("optimisation", {}).get("geometry_version", 0)
        geometry_version = prev_version + 1

        # Last modified timestamp
        last_modified = datetime.utcnow().isoformat() + "Z"

        # Update conformer record
        rec.energy = energy
        rec.xyz_path = xyz_out

        rec.provenance["optimisation"] = {
            "engine": engine,
            "level": level,
            "status": status,
            "energy": energy,
            "xyz_path": xyz_out,
            "log_path": log_out,
            "elapsed_seconds": elapsed_wall,
            "parser": {
                "iterations": parsed.get("iterations"),
                "converged": parsed.get("converged"),
                "elapsed_seconds": parsed.get("elapsed_seconds"),
            },
            "backend_meta": backend_meta,
            "engine_command": engine_command,
            "backend_version": backend_version,
            "warnings": warnings_list,
            "convergence_quality": quality,
            "success_score": score,
            "geometry_version": geometry_version,
            "last_modified": last_modified,
        }

        return rec

    # -------------------------------------------------------------------------
    # Engine dispatch
    # -------------------------------------------------------------------------
    def _run_engine(self, xyz_path, engine, level, max_iter):
        if engine not in OPTIMISATION_BACKENDS:
            self.fail(f"Unknown optimisation engine: {engine}")

        method_name, legacy_level = OPTIMISATION_BACKENDS[engine]
        runner = getattr(self, method_name)

        if method_name == "_run_orca_opt":
            return runner(xyz_path, level, max_iter)

        return runner(xyz_path, level, max_iter)

    # -------------------------------------------------------------------------
    # Forcefield placeholder
    # -------------------------------------------------------------------------
    def _run_forcefield(self, xyz_path, level, max_iter):
        job = self.job
        lookup = os.path.splitext(os.path.basename(xyz_path))[0]

        log_out = os.path.join(job.outputs_dir, "log", f"{lookup}_ff.log")
        xyz_out = os.path.join(job.outputs_dir, "xyz", f"{lookup}_opt.xyz")

        os.makedirs(os.path.dirname(log_out), exist_ok=True)
        os.makedirs(os.path.dirname(xyz_out), exist_ok=True)

        with open(log_out, "w") as f:
            f.write("Forcefield optimisation placeholder\n")

        shutil.copy(xyz_path, xyz_out)

        return {
            "log_out": log_out,
            "xyz_out": xyz_out,
            "backend_meta": {"mode": "forcefield"},
            "engine_command": "forcefield placeholder",
            "backend_version": "unknown",
        }

    # -------------------------------------------------------------------------
    # ORCA backend
    # -------------------------------------------------------------------------
    def _run_orca_opt(self, xyz_path, level, max_iter):
        job = self.job
        lookup = os.path.splitext(os.path.basename(xyz_path))[0]

        scratch_root = os.path.join(job.job_dir, "orca_opt_scratch")
        scratch_dir = os.path.join(scratch_root, lookup)
        os.makedirs(scratch_dir, exist_ok=True)

        inp_file = os.path.join(scratch_dir, f"{lookup}.inp")
        log_out = os.path.join(job.outputs_dir, "log", f"{lookup}_orca.log")
        xyz_out = os.path.join(job.outputs_dir, "xyz", f"{lookup}_opt.xyz")

        os.makedirs(os.path.dirname(log_out), exist_ok=True)
        os.makedirs(os.path.dirname(xyz_out), exist_ok=True)

        # Level mapping
        if level == "loose":
            opt_keyword = "LooseOpt"
        elif level == "normal":
            opt_keyword = "Opt"
        elif level == "tight":
            opt_keyword = "TightOpt"
        else:
            opt_keyword = "VeryTightOpt"

        method_line = f"! BP86 def2-TZVP {opt_keyword}"

        with open(inp_file, "w") as f:
            f.write(method_line + "\n")
            f.write("* xyz 0 1\n")
            with open(xyz_path) as xyz:
                lines = xyz.readlines()
                f.writelines(lines[2:])
            f.write("*\n")

        backend_meta = {
            "mode": "orca",
            "level": level,
            "scratch_dir": scratch_dir,
            "run_status": "unknown",
        }

        engine_command = f"orca {inp_file}"

        try:
            with open(log_out, "w") as log:
                subprocess.run(
                    ["orca", inp_file],
                    cwd=scratch_dir,
                    stdout=log,
                    stderr=log,
                    check=True,
                )
            backend_meta["run_status"] = "ok"

            final_xyz = os.path.join(scratch_dir, f"{lookup}.xyz")
            if os.path.exists(final_xyz):
                shutil.copy(final_xyz, xyz_out)
            else:
                self.log(f"[WARNING] ORCA did not produce final XYZ for {lookup}")
                xyz_out = None

        except subprocess.CalledProcessError:
            backend_meta["run_status"] = "failed"
            xyz_out = None

        backend_version = "unknown"
        try:
            with open(log_out) as f:
                for line in f:
                    if "Program Version" in line:
                        backend_version = line.strip()
                        break
        except Exception:
            pass

        return {
            "log_out": log_out if os.path.exists(log_out) else None,
            "xyz_out": xyz_out if os.path.exists(xyz_out) else None,
            "backend_meta": backend_meta,
            "engine_command": engine_command,
            "backend_version": backend_version,
        }

    # -------------------------------------------------------------------------
    # gXTB backend
    # -------------------------------------------------------------------------
    def _run_gxtb(self, xyz_path, level, max_iter):
        job = self.job
        lookup = os.path.splitext(os.path.basename(xyz_path))[0]

        xtb_bin = self.config["xtb"]["executable"]
        gxtb_bin = self.config["gxtb"]["executable"]

        tmp_exec = os.path.join(job.job_dir, "gxtb_scratch")
        os.makedirs(tmp_exec, exist_ok=True)

        tmp_xyz = os.path.join(tmp_exec, f"{lookup}.xyz")
        opt_tmp = os.path.join(tmp_exec, "xtbopt.xyz")

        for fpath in (tmp_xyz, opt_tmp):
            if os.path.exists(fpath):
                os.remove(fpath)

        shutil.copy(xyz_path, tmp_xyz)

        log_out = os.path.join(job.outputs_dir, "log", f"{lookup}_gxtb.log")
        xyz_out = os.path.join(job.outputs_dir, "xyz", f"{lookup}_opt.xyz")

        os.makedirs(os.path.dirname(log_out), exist_ok=True)
        os.makedirs(os.path.dirname(xyz_out), exist_ok=True)

        # Level mapping
        if level == "loose":
            opt_flag = "--opt loose"
        elif level == "normal":
            opt_flag = "--opt"
        elif level == "tight":
            opt_flag = "--opt tight"
        else:
            opt_flag = "--opt vtight"

        max_cores = os.cpu_count() or 1
        ncores = max(1, math.floor(max_cores * 0.8))

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(ncores)

        driver_string = f"{gxtb_bin} -grad -c xtbdriver.xyz"

        cmd = [
            xtb_bin,
            tmp_xyz,
            "--driver", driver_string,
            opt_flag,
            "--iterations", str(max_iter),
        ]

        engine_command = " ".join(cmd)

        backend_meta = {
            "mode": "gxtb",
            "level": level,
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
                self.log(f"[WARNING] gXTB did not produce xtbopt.xyz for {lookup}")
                xyz_out = None

        except subprocess.CalledProcessError:
            backend_meta["run_status"] = "failed"
            xyz_out = None

        backend_version = "unknown"
        try:
            version_cmd = [xtb_bin, "--version"]
            out = subprocess.check_output(version_cmd, stderr=subprocess.STDOUT)
            backend_version = out.decode().strip()
        except Exception:
            pass

        if os.path.exists(tmp_xyz):
            os.remove(tmp_xyz)

        return {
            "log_out": log_out if os.path.exists(log_out) else None,
            "xyz_out": xyz_out if os.path.exists(xyz_out) else None,
            "backend_meta": backend_meta,
            "engine_command": engine_command,
            "backend_version": backend_version,
        }

    # -------------------------------------------------------------------------
    # Classic XTB backend
    # -------------------------------------------------------------------------
    def _run_xtb(self, xyz_path, level, max_iter):
        job = self.job
        lookup = os.path.splitext(os.path.basename(xyz_path))[0]

        scratch = os.path.join(job.job_dir, "xtb_scratch")
        os.makedirs(scratch, exist_ok=True)

        tmp_xyz = os.path.join(scratch, f"{lookup}.xyz")
        opt_tmp = os.path.join(scratch, "xtbopt.xyz")

        for fpath in (tmp_xyz, opt_tmp):
            if os.path.exists(fpath):
                os.remove(fpath)

        shutil.copy(xyz_path, tmp_xyz)

        log_out = os.path.join(job.outputs_dir, "log", f"{lookup}_xtb.log")
        xyz_out = os.path.join(job.outputs_dir, "xyz", f"{lookup}_opt.xyz")

        os.makedirs(os.path.dirname(log_out), exist_ok=True)
        os.makedirs(os.path.dirname(xyz_out), exist_ok=True)

        # Level mapping via iterations
        if level == "loose":
            iters = 100
        elif level == "normal":
            iters = 250
        elif level == "tight":
            iters = 500
        else:
            iters = 1000

        cmd = [
            "xtb", tmp_xyz,
            "--opt",
            "--gfn", job.parameters.get("level", "2"),
            "--iterations", str(iters),
        ]

        engine_command = " ".join(cmd)

        backend_meta = {
            "mode": "xtb",
            "level": level,
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
                self.log(f"[WARNING] XTB did not produce xtbopt.xyz for {lookup}")
                xyz_out = None

        except subprocess.CalledProcessError:
            backend_meta["run_status"] = "failed"
            xyz_out = None

        backend_version = "unknown"
        try:
            out = subprocess.check_output(["xtb", "--version"], stderr=subprocess.STDOUT)
            backend_version = out.decode().strip()
        except Exception:
            pass

        if os.path.exists(tmp_xyz):
            os.remove(tmp_xyz)

        return {
            "log_out": log_out if os.path.exists(log_out) else None,
            "xyz_out": xyz_out if os.path.exists(xyz_out) else None,
            "backend_meta": backend_meta,
            "engine_command": engine_command,
            "backend_version": backend_version,
        }

    # -------------------------------------------------------------------------
    # Unified output writer
    # -------------------------------------------------------------------------
    def _write_outputs(self, records):
        outputs_dir = self.outputs_dir
        outputs_xyz_dir = os.path.join(outputs_dir, "xyz")
        os.makedirs(outputs_xyz_dir, exist_ok=True)

        # Copy XYZs into outputs/xyz
        for rec in records:
            opt = rec.provenance["optimisation"]
            src = opt.get("xyz_path")

            if not src or not os.path.exists(src):
                self.log(f"[WARNING] Missing xyz_path for {rec.lookup_id}, skipping copy.")
                rec.xyz_path = None
                continue

            dst = os.path.join(outputs_xyz_dir, f"{rec.lookup_id}_opt.xyz")
            try:
                shutil.copy(src, dst)
                rec.xyz_path = dst
                opt["xyz_path"] = dst
            except Exception as e:
                self.log(f"[ERROR] Failed to copy {src} → {dst}: {e}")
                rec.xyz_path = None
                opt["xyz_path"] = None

        # summary.csv (per conformer)
        summary_csv = os.path.join(outputs_dir, "summary.csv")
        df = pd.DataFrame(
            [
                {
                    "lookup_id": r.lookup_id,
                    "inchi_key": r.inchi_key,
                    "conf_num": r.conf_num,
                    "energy": r.energy,
                    "status": r.provenance["optimisation"]["status"],
                    "quality": r.provenance["optimisation"]["convergence_quality"],
                    "success_score": r.provenance["optimisation"]["success_score"],
                    "xyz_path": r.xyz_path,
                    "log_path": r.provenance["optimisation"]["log_path"],
                    "elapsed_seconds": r.provenance["optimisation"]["elapsed_seconds"],
                }
                for r in records
            ]
        )
        df.to_csv(summary_csv, index=False)

        # energies.json (canonical)
        energies_json = os.path.join(outputs_dir, "energies.json")
        pruned_set = ConformerSet(records)
        pruned_set.save(energies_json)

        self.log(f"Optimisation outputs written to: {outputs_dir}")

    # -------------------------------------------------------------------------
    # Human-readable molecule summary
    # -------------------------------------------------------------------------
    def _write_human_summary_csv(self, records, molecules, engine, level, outputs_dir):
        summary_csv = os.path.join(outputs_dir, "optimisation_summary.csv")
        mol_stats = []

        by_mol = {}
        for r in records:
            by_mol.setdefault(r.inchi_key, []).append(r)

        for mol_id, mol_entries in molecules.items():
            mol_records = by_mol.get(mol_id, [])

            n_confs = len(mol_records)
            n_success = sum(1 for r in mol_records if r.provenance["optimisation"]["status"] == "converged")
            n_failed = sum(1 for r in mol_records if r.provenance["optimisation"]["status"] == "failed")
            n_partial = sum(1 for r in mol_records if r.provenance["optimisation"]["status"] == "partial")

            times = [r.provenance["optimisation"]["elapsed_seconds"] for r in mol_records]
            avg_time = float(np.mean(times)) if times else 0.0
            total_time = float(np.sum(times)) if times else 0.0

            n_atoms = None
            for r in mol_records:
                if r.xyz_path and os.path.exists(r.xyz_path):
                    try:
                        with open(r.xyz_path) as f:
                            n_atoms = int(f.readline().strip())
                        break
                    except Exception:
                        pass

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
                    "level": level,
                }
            )

        pd.DataFrame(mol_stats).to_csv(summary_csv, index=False)
        self.log(f"Human-readable optimisation summary written to: {summary_csv}")

    # -------------------------------------------------------------------------
    # Warning summary
    # -------------------------------------------------------------------------
    def _log_warning_summary(self):
        self.log_header("Optimisation Summary")

        # Count statuses
        converged = 0
        partial = 0
        failed = 0
        missing_xyz = 0

        # We can only count from outputs/summary.csv
        summary_csv = os.path.join(self.outputs_dir, "summary.csv")
        if not os.path.exists(summary_csv):
            self.log("[WARNING] summary.csv missing; cannot compute warning summary.")
            return

        df = pd.read_csv(summary_csv)

        converged = (df["status"] == "converged").sum()
        partial = (df["status"] == "partial").sum()
        failed = (df["status"] == "failed").sum()
        missing_xyz = df["xyz_path"].isna().sum()

        self.log(f"Converged: {converged}")
        self.log(f"Partial: {partial}")
        self.log(f"Failed: {failed}")
        self.log(f"Missing XYZ: {missing_xyz}")

        self.log_header("Optimisation Stage Complete")
