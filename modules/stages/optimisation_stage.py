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
# Backend dispatch table (new API)
# -------------------------------------------------------------------------
OPT_BACKENDS = {
    "xtb": "_backend_xtb",
    "gxtb": "_backend_gxtb",
    "orca": "_backend_orca",
    "orca_fast": "_backend_orca",
    "orca_final": "_backend_orca",
    "forcefield": "_backend_forcefield",
}

PARSERS = {
    "xtb": XTBLogParser,
    "gxtb": GxTBLogParser,
    "orca": ORCALogParser,
    "orca_fast": ORCALogParser,
    "orca_final": ORCALogParser,
}

VALID_LEVELS = ("loose", "normal", "tight", "vtight")


# =========================================================================
#  NEW OPTIMISATION STAGE
# =========================================================================
class OptimisationStage(BaseStage):

    # ---------------------------------------------------------------------
    # Entry point
    # ---------------------------------------------------------------------
    def execute(self):
        self.strict_mode = self.strict("optimisation")
        self.set_stage_output("energies.json")

        # Load input energies.json
        energies_file = self.require_file(self.get_stage_input(), "stage_input energies.json")

        # Copy into inputs/ for reproducibility
        inputs_energies = os.path.join(self.inputs_dir, "energies.json")
        if os.path.abspath(energies_file) != os.path.abspath(inputs_energies):
            shutil.copy(energies_file, inputs_energies)
        self.log(f"Copied energies.json → {inputs_energies}")

        # Parameters
        params = self.parameters or {}
        engine = params.get("engine", "gxtb").lower()
        level = params.get("level", "normal").lower()
        max_iter = params.get("max_iter", 250)
        global_fail_threshold = params.get("global_fail_threshold", 0.8)
        keep_scratch = params.get("keep_scratch", False)
        gfn_level = params.get("gfn", 2)

        if level not in VALID_LEVELS:
            self.fail(f"Unknown optimisation level: {level}")

        self.log(f"Engine: {engine}, level={level}, max_iter={max_iter}")
        self.log(f"Strict mode: {self.strict_mode}")
        self.log(f"Keep scratch: {keep_scratch}")

        # Load conformers
        conformer_set = ConformerSet.load(energies_file)
        entries = conformer_set.records

        # Ensure optimisation_history exists
        for rec in entries:
            if not hasattr(rec, "optimisation_history"):
                rec.optimisation_history = []

        # Prepare XYZs
        prep, valid_records = self._prepare_inputs(entries)
        self.set_items([rec.lookup_id for rec in valid_records])
        self.log(f"Valid conformers: {len(valid_records)}")

        # Group by molecule
        molecules = self._group_by_molecule(valid_records)

        # Optimise
        optimised_records = self._optimise_all(
            valid_records,
            engine,
            level,
            max_iter,
            global_fail_threshold,
            keep_scratch,
            gfn_level,
        )

        # Write outputs
        self._write_outputs(optimised_records)

        # Molecule summary
        self._write_human_summary_csv(
            optimised_records,
            molecules,
            engine,
            level,
            self.outputs_dir,
        )

        # Metadata
        meta_path = os.path.join(self.outputs_dir, "optimisation_metadata.json")
        with AtomicWriter(meta_path) as f:
            json.dump(
                {
                    "stage": "optimisation",
                    "engine": engine,
                    "level": level,
                    "max_iter": max_iter,
                    "num_input": len(entries),
                    "num_output": len(optimised_records),
                    "missing_xyz": prep["missing"],
                    "timestamp": datetime.utcnow().isoformat() + "Z",
                },
                f,
                indent=2,
            )

    # ---------------------------------------------------------------------
    # Input preparation
    # ---------------------------------------------------------------------
    def _prepare_inputs(self, records):
        inputs_xyz = os.path.join(self.inputs_dir, "xyz")
        os.makedirs(inputs_xyz, exist_ok=True)

        missing = []
        valid = []

        for rec in records:
            src = rec.xyz_path
            if not src or not os.path.exists(src):
                self.log(f"[WARNING] Missing XYZ for {rec.lookup_id}: {src}")
                missing.append(rec.lookup_id)
                continue

            dst = os.path.join(inputs_xyz, os.path.basename(src))
            if os.path.abspath(src) != os.path.abspath(dst):
                shutil.copy(src, dst)

            rec.xyz_path = dst
            valid.append(rec)

        if not valid:
            self.fail("No valid conformers available for optimisation.")

        return {"missing": missing}, valid

    # ---------------------------------------------------------------------
    # Group by molecule
    # ---------------------------------------------------------------------
    def _group_by_molecule(self, records):
        groups = {}
        for rec in records:
            groups.setdefault(rec.inchi_key, []).append(rec)
        return groups

    # ---------------------------------------------------------------------
    # Optimise all conformers
    # ---------------------------------------------------------------------
    def _optimise_all(
        self,
        records,
        engine,
        level,
        max_iter,
        global_fail_threshold,
        keep_scratch,
        gfn_level,
    ):
        optimised = []
        total = len(records)
        failures = 0

        for idx, rec in enumerate(records, start=1):
            self.log(f"[{idx}/{total}] Optimising {rec.lookup_id}", indent=2)

            try:
                updated = self._optimise_single(
                    rec, engine, level, max_iter, keep_scratch, gfn_level
                )
                optimised.append(updated)

                last = updated.optimisation_history[-1]
                status = last["status"]
                self.update_progress(rec.lookup_id, success=(status == "converged"))

                if status != "converged":
                    failures += 1

            except Exception as e:
                failures += 1
                self.update_progress(rec.lookup_id, success=False)
                self.log(f"[ERROR] Exception during optimisation of {rec.lookup_id}: {e}")

            if total > 0 and failures / total >= global_fail_threshold:
                self.log(
                    f"[ERROR] Global failure threshold reached "
                    f"({failures}/{total}); aborting optimisation."
                )
                break

        usable = [r for r in optimised if self._is_record_usable(r)]

        if not usable:
            self.fail("Optimisation produced zero usable conformers.")

        dropped = len(optimised) - len(usable)
        if dropped > 0:
            self.log(f"[WARNING] Dropped {dropped} unusable conformers.")

        return usable

    # ---------------------------------------------------------------------
    # Usability check
    # ---------------------------------------------------------------------
    def _is_record_usable(self, rec):
        if not rec.optimisation_history:
            return False

        last = rec.optimisation_history[-1]
        xyz = last["xyz_path"]
        energy = last["energy"]

        if not xyz or not os.path.exists(xyz):
            return False

        if energy is None or math.isnan(energy) or abs(energy) > 1e6:
            return False

        return True

    # ---------------------------------------------------------------------
    # Optimise a single conformer
    # ---------------------------------------------------------------------
    def _optimise_single(self, rec, engine, level, max_iter, keep_scratch, gfn_level):
        lookup = rec.lookup_id
        input_xyz = rec.xyz_path

        # Determine iteration number
        iteration = len(rec.optimisation_history) 

        # Output paths
        output_xyz = os.path.join(self.outputs_dir, "xyz", f"{lookup}_opt{iteration}.xyz")
        output_log = os.path.join(self.outputs_dir, "log", f"{lookup}_opt{iteration}.log")

        os.makedirs(os.path.dirname(output_xyz), exist_ok=True)
        os.makedirs(os.path.dirname(output_log), exist_ok=True)

        # Scratch directory
        scratch_dir = os.path.join(
            self.job.job_dir,
            "scratch",
            engine,
            lookup,
            f"opt{iteration}",
        )
        os.makedirs(scratch_dir, exist_ok=True)

        # Run backend
        start_wall = time.perf_counter()
        start_cpu = time.process_time()

        backend_result = self._run_backend(
            input_xyz=input_xyz,
            output_xyz=output_xyz,
            output_log=output_log,
            scratch_dir=scratch_dir,
            engine=engine,
            level=level,
            max_iter=max_iter,
            gfn_level=gfn_level,
        )

        end_wall = time.perf_counter()
        end_cpu = time.process_time()

        # Parse log
        parsed = self._parse_backend_log(output_log, engine)

        # Determine status
        status = self._determine_status(parsed, output_xyz)

        # Validate energy
        energy = self._validate_energy(parsed.get("energy"))

        # Assign quality + score
        quality = self._assign_quality(status, level)
        score = self._assign_score(status, level)

        # Update history
        self._update_history(
            rec,
            engine,
            level,
            status,
            energy,
            output_xyz,
            output_log,
            backend_result,
            quality,
            score,
            start_wall,
            end_wall,
            start_cpu,
            end_cpu,
        )

        # Update conformer
        rec.energy = energy
        rec.xyz_path = output_xyz if os.path.exists(output_xyz) else None

        # Cleanup scratch
        if not keep_scratch:
            shutil.rmtree(scratch_dir, ignore_errors=True)

        return rec

    # ---------------------------------------------------------------------
    # Backend runner
    # ---------------------------------------------------------------------
    def _run_backend(
        self,
        input_xyz,
        output_xyz,
        output_log,
        scratch_dir,
        engine,
        level,
        max_iter,
        gfn_level,
    ):
        if engine not in OPT_BACKENDS:
            self.fail(f"Unknown optimisation engine: {engine}")

        method = getattr(self, OPT_BACKENDS[engine])

        return method(
            input_xyz=input_xyz,
            output_xyz=output_xyz,
            output_log=output_log,
            scratch_dir=scratch_dir,
            level=level,
            max_iter=max_iter,
            gfn_level=gfn_level,
        )

    # ---------------------------------------------------------------------
    # Backend: Forcefield (placeholder)
    # ---------------------------------------------------------------------
    def _backend_forcefield(self, input_xyz, output_xyz, output_log, scratch_dir, **kwargs):
        with open(output_log, "w") as f:
            f.write("Forcefield optimisation placeholder\n")

        shutil.copy(input_xyz, output_xyz)

        return {
            "run_status": "ok",
            "version": "unknown",  # placeholder
            "command": [],         # placeholder
            "command_str": "",     # placeholder
        }

    # ---------------------------------------------------------------------
    # Backend: ORCA
    # ---------------------------------------------------------------------
    def _backend_orca(self, input_xyz, output_xyz, output_log, scratch_dir, level, max_iter, **kwargs):
        lookup = os.path.splitext(os.path.basename(input_xyz))[0]
        inp_file = os.path.join(scratch_dir, f"{lookup}.inp")

        opt_keyword = {
            "loose": "LooseOpt",
            "normal": "Opt",
            "tight": "TightOpt",
            "vtight": "VeryTightOpt",
        }[level]

        method_line = f"! BP86 def2-TZVP {opt_keyword}"

        # Write ORCA input
        with open(inp_file, "w") as f:
            f.write(method_line + "\n")
            f.write("* xyz 0 1\n")
            with open(input_xyz) as xyz:
                lines = xyz.readlines()
                f.writelines(lines[2:])
            f.write("*\n")

        cmd = ["orca", inp_file]

        try:
            with open(output_log, "w") as log:
                subprocess.run(
                    cmd,
                    cwd=scratch_dir,
                    stdout=log,
                    stderr=log,
                    check=True,
                )

            final_xyz = os.path.join(scratch_dir, f"{lookup}.xyz")
            if os.path.exists(final_xyz):
                shutil.copy(final_xyz, output_xyz)
            else:
                output_xyz = None

            run_status = "ok"

        except subprocess.CalledProcessError:
            run_status = "failed"
            output_xyz = None

        return {
            "run_status": run_status,
            "version": "unknown",  # placeholder
            "command": cmd,
            "command_str": " ".join(cmd),
        }

    # ---------------------------------------------------------------------
    # Backend: gXTB
    # ---------------------------------------------------------------------
    def _backend_gxtb(self, input_xyz, output_xyz, output_log, scratch_dir, level, max_iter, gfn_level, **kwargs):
        xtb_bin = self.config["xtb"]["executable"]
        gxtb_bin = self.config["gxtb"]["executable"]

        tmp_xyz = os.path.join(scratch_dir, "input.xyz")
        opt_tmp = os.path.join(scratch_dir, "xtbopt.xyz")

        if os.path.exists(tmp_xyz):
            os.remove(tmp_xyz)
        if os.path.exists(opt_tmp):
            os.remove(opt_tmp)

        shutil.copy(input_xyz, tmp_xyz)

        opt_flag = {
            "loose": "--opt loose",
            "normal": "--opt",
            "tight": "--opt tight",
            "vtight": "--opt vtight",
        }[level]

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

        try:
            with open(output_log, "w") as log:
                subprocess.run(
                    cmd,
                    cwd=scratch_dir,
                    stdout=log,
                    stderr=log,
                    check=True,
                    env=env,
                )

            if os.path.exists(opt_tmp):
                shutil.copy(opt_tmp, output_xyz)
                run_status = "ok"
            else:
                run_status = "failed"
                output_xyz = None

        except subprocess.CalledProcessError:
            run_status = "failed"
            output_xyz = None

        return {
            "run_status": run_status,
            "version": "unknown",  # placeholder
            "command": cmd,
            "command_str": " ".join(cmd),
        }

    # ---------------------------------------------------------------------
    # Backend: classic XTB
    # ---------------------------------------------------------------------
    def _backend_xtb(self, input_xyz, output_xyz, output_log, scratch_dir, level, max_iter, gfn_level, **kwargs):
        tmp_xyz = os.path.join(scratch_dir, "input.xyz")
        opt_tmp = os.path.join(scratch_dir, "xtbopt.xyz")

        if os.path.exists(tmp_xyz):
            os.remove(tmp_xyz)
        if os.path.exists(opt_tmp):
            os.remove(opt_tmp)

        shutil.copy(input_xyz, tmp_xyz)

        if level == "loose":
            iters = 100
        elif level == "normal":
            iters = 250
        elif level == "tight":
            iters = 500
        else:
            iters = 1000

        cmd = [
            "xtb",
            tmp_xyz,
            "--opt",
            "--gfn", str(gfn_level),
            "--iterations", str(iters),
        ]

        try:
            with open(output_log, "w") as log:
                subprocess.run(
                    cmd,
                    cwd=scratch_dir,
                    stdout=log,
                    stderr=log,
                    check=True,
                )

            if os.path.exists(opt_tmp):
                shutil.copy(opt_tmp, output_xyz)
                run_status = "ok"
            else:
                run_status = "failed"
                output_xyz = None

        except subprocess.CalledProcessError:
            run_status = "failed"
            output_xyz = None

        return {
            "run_status": run_status,
            "version": "unknown",  # placeholder
            "command": cmd,
            "command_str": " ".join(cmd),
        }

    # ---------------------------------------------------------------------
    # Log parsing
    # ---------------------------------------------------------------------
    def _parse_backend_log(self, log_path, engine):
        parser_cls = PARSERS.get(engine)
        if not parser_cls:
            return {"energy": None, "converged": None}

        if not log_path or not os.path.exists(log_path):
            return {"energy": None, "converged": None}

        try:
            return parser_cls.parse(log_path)
        except Exception as e:
            self.log(f"[WARNING] Parser failed for {log_path}: {e}")
            return {"energy": None, "converged": None}

    # ---------------------------------------------------------------------
    # Status determination
    # ---------------------------------------------------------------------
    def _determine_status(self, parsed, xyz_path):
        conv = parsed.get("converged")
        xyz_exists = xyz_path and os.path.exists(xyz_path)

        if conv is True:
            return "converged" if xyz_exists else "partial"
        if conv is False:
            return "failed"
        return "partial"

    # ---------------------------------------------------------------------
    # Energy validation
    # ---------------------------------------------------------------------
    def _validate_energy(self, energy):
        if energy is None:
            return None
        try:
            val = float(energy)
            if math.isnan(val) or abs(val) > 1e6:
                return None
            return val
        except Exception:
            return None

    # ---------------------------------------------------------------------
    # Quality scoring
    # ---------------------------------------------------------------------
    def _assign_quality(self, status, level):
        if status == "converged":
            return {
                "loose": "medium",
                "normal": "good",
                "tight": "very_good",
                "vtight": "excellent",
            }.get(level, "good")

        if status == "partial":
            return "poor"

        return "unusable"

    # ---------------------------------------------------------------------
    # Success score
    # ---------------------------------------------------------------------
    def _assign_score(self, status, level):
        if status == "converged":
            return {
                "loose": 0.4,
                "normal": 0.6,
                "tight": 0.8,
                "vtight": 1.0,
            }.get(level, 0.6)

        if status == "partial":
            return 0.2

        return 0.0

    # ---------------------------------------------------------------------
    # Update optimisation history
    # ---------------------------------------------------------------------
    def _update_history(
        self,
        rec,
        engine,
        level,
        status,
        energy,
        xyz_path,
        log_path,
        backend_result,
        quality,
        score,
        start_wall,
        end_wall,
        start_cpu,
        end_cpu,
    ):
        rec.optimisation_history.append(
            {
                "stage": "optimisation",
                "engine": engine,
                "level": level,
                "status": status,
                "energy": energy,
                "xyz_path": xyz_path,
                "log_path": log_path,
                "elapsed_seconds": end_wall - start_wall,
                "timing": {
                    "start_wall": start_wall,
                    "end_wall": end_wall,
                    "start_cpu": start_cpu,
                    "end_cpu": end_cpu,
                    "wall_seconds": end_wall - start_wall,
                    "cpu_seconds": end_cpu - start_cpu,
                },
                "backend_meta": {
                    "run_status": backend_result.get("run_status"),
                    "version": backend_result.get("version", "unknown"),
                    "command": backend_result.get("command", []),
                    "command_str": backend_result.get("command_str", ""),
                },
                "convergence_quality": quality,
                "success_score": score,
                "timestamp": datetime.utcnow().isoformat() + "Z",
            }
        )

    # ---------------------------------------------------------------------
    # Write outputs (now trivial)
    # ---------------------------------------------------------------------
    def _write_outputs(self, records):
        outputs_dir = self.outputs_dir

        # summary.csv
        summary_csv = os.path.join(outputs_dir, "summary.csv")
        rows = []

        for r in records:
            last = r.optimisation_history[-1]
            rows.append(
                {
                    "lookup_id": r.lookup_id,
                    "inchi_key": r.inchi_key,
                    "conf_num": r.conf_num,
                    "energy": last["energy"],
                    "status": last["status"],
                    "quality": last["convergence_quality"],
                    "success_score": last["success_score"],
                    "xyz_path": last["xyz_path"],
                    "log_path": last["log_path"],
                    "elapsed_seconds": last["elapsed_seconds"],
                    "engine": last["engine"],
                    "level": last["level"],
                }
            )

        pd.DataFrame(rows).to_csv(summary_csv, index=False)

        # energies.json (canonical)
        energies_json = self.get_stage_output()
        ConformerSet(records).save(energies_json)

        self.log(f"Optimisation outputs written to: {outputs_dir}")

    # ---------------------------------------------------------------------
    # Molecule-level summary
    # ---------------------------------------------------------------------
    def _write_human_summary_csv(self, records, molecules, engine, level, outputs_dir):
        summary_csv = os.path.join(outputs_dir, "optimisation_summary.csv")
        mol_stats = []

        by_mol = {}
        for r in records:
            by_mol.setdefault(r.inchi_key, []).append(r)

        for mol_id, mol_entries in molecules.items():
            mol_records = by_mol.get(mol_id, [])

            n_confs = len(mol_records)
            n_success = 0
            n_failed = 0
            n_partial = 0
            times = []
            n_atoms = None

            for r in mol_records:
                last = r.optimisation_history[-1]
                status = last["status"]

                if status == "converged":
                    n_success += 1
                elif status == "failed":
                    n_failed += 1
                else:
                    n_partial += 1

                times.append(last["elapsed_seconds"])

                if last["xyz_path"] and os.path.exists(last["xyz_path"]) and n_atoms is None:
                    try:
                        with open(last["xyz_path"]) as f:
                            n_atoms = int(f.readline().strip())
                    except Exception:
                        pass

            avg_time = float(np.mean(times)) if times else 0.0
            total_time = float(np.sum(times)) if times else 0.0

            mol_stats.append(
                {
                    "molecule_id": mol_id,
                    "n_atoms": n_atoms,
                    "n_conformers": n_confs,
                    "n_converged": n_success,
                    "n_failed": n_failed,
                    "n_partial": n_partial,
                    "avg_time_seconds": avg_time,
                    "total_time_seconds": total_time,
                    "engine": engine,
                    "level": level,
                }
            )

        pd.DataFrame(mol_stats).to_csv(summary_csv, index=False)
