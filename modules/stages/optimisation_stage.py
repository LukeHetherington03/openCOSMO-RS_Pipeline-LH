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


VALID_LEVELS = ("loose", "normal", "tight", "vtight")

PARSERS = {
    "orca": ORCALogParser,
    "xtb": XTBLogParser,
    "gxtb": GxTBLogParser,
}


# =========================================================================
#  OPTIMISATION STAGE (FINAL VERSION)
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

        # Load conformers
        conformer_set = ConformerSet.load(energies_file)
        all_records = conformer_set.records

        # Load optimisation config
        defaults_path = self.config["optimisation"]["defaults"]
        engines_path = self.config["optimisation"]["engines"]

        with open(defaults_path) as f:
            defaults = json.load(f)
        with open(engines_path) as f:
            engine_registry = json.load(f)

        # Merge defaults + user overrides
        params = dict(defaults)
        if self.parameters:
            params.update(self.parameters)

        # Engine selection
        self.engine_name = params.get("engine")
        if self.engine_name not in engine_registry:
            self.fail(f"Unknown optimisation engine: {self.engine_name}")

        self.engine_spec = engine_registry[self.engine_name]

        # Stage-level parameters
        self.level = self.engine_spec.get("level", params.get("level", "normal"))
        self.max_iter = params.get("max_iter", 250)
        self.gfn_level = self.engine_spec.get("gfn", params.get("gfn", 2))
        self.keep_scratch = params.get("keep_scratch", False)
        self.global_fail_threshold = params.get("global_fail_threshold", 0.8)

        # Resume logic
        resume = params.get("resume", True)

        if not resume or not self.job.pending_items:
            self.set_items([rec.lookup_id for rec in all_records])
            self._clear_checkpoints()
            self.log("[INFO] Starting optimisation from scratch")
        else:
            self.log("[INFO] Resuming optimisation")
            self.log(f"[INFO] Pending conformers: {len(self.job.pending_items)}")

        # Filter pending conformers
        pending = [rec for rec in all_records if rec.lookup_id in self.job.pending_items]

        if not pending:
            self.log("No pending conformers — optimisation stage already complete.")
            return

        # Ensure optimisation_history exists
        for rec in pending:
            if not hasattr(rec, "optimisation_history"):
                rec.optimisation_history = []

        # Prepare XYZs
        prep, valid_records = self._prepare_inputs(pending)
        self.log(f"Valid conformers: {len(valid_records)}")

        # Run optimisation
        optimised_records = self._optimise_all(valid_records)

        # Merge checkpoints → final energies.json
        final_set = self._merge_checkpoints()
        final_set.save(self.get_stage_output())

        self.log(f"[INFO] Optimisation complete. Output written to {self.get_stage_output()}")
        self.job.mark_complete()


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
    def _optimise_all(self, records):
        optimised = []
        total = len(records)
        failures = 0

        for idx, rec in enumerate(records, start=1):
            self.log(f"[{idx}/{total}] Optimising {rec.lookup_id}", indent=2)

            try:
                updated = self._optimise_single(rec)
                optimised.append(updated)

                last = updated.optimisation_history[-1]
                status = last["status"]

                self.update_progress(rec.lookup_id, success=(status == "converged"))

                if status != "converged":
                    failures += 1

                # Write checkpoint
                self._write_checkpoint(rec.lookup_id, updated)

            except Exception as e:
                failures += 1
                self.update_progress(rec.lookup_id, success=False)
                self.log(f"[ERROR] Exception during optimisation of {rec.lookup_id}: {e}")

            # Global failure threshold
            if total > 0 and failures / total >= self.global_fail_threshold:
                self.log(f"[ERROR] Global failure threshold reached ({failures}/{total}); aborting optimisation.")
                break

        usable = [r for r in optimised if self._is_record_usable(r)]

        if not usable:
            self.fail("Optimisation produced zero usable conformers.")

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
    def _optimise_single(self, rec):
        lookup = rec.lookup_id
        input_xyz = rec.xyz_path

        iteration = len(rec.optimisation_history)

        output_xyz = os.path.join(self.outputs_dir, "xyz", f"{lookup}_opt{iteration}.xyz")
        output_log = os.path.join(self.outputs_dir, "log", f"{lookup}_opt{iteration}.log")

        os.makedirs(os.path.dirname(output_xyz), exist_ok=True)
        os.makedirs(os.path.dirname(output_log), exist_ok=True)

        scratch_dir = os.path.join(self.job.job_dir, "scratch", lookup, f"opt{iteration}")
        os.makedirs(scratch_dir, exist_ok=True)

        start_wall = time.perf_counter()
        start_cpu = time.process_time()

        backend_result = self._run_backend(
            input_xyz=input_xyz,
            output_xyz=output_xyz,
            output_log=output_log,
            scratch_dir=scratch_dir,
            engine_name=self.engine_name,
            engine_spec=self.engine_spec,
            level=self.level,
            max_iter=self.max_iter,
            gfn_level=self.gfn_level,
        )

        end_wall = time.perf_counter()
        end_cpu = time.process_time()

        parsed = self._parse_backend_log(output_log, self.engine_spec)
        status = self._determine_status(parsed, output_xyz)
        energy = self._validate_energy(parsed.get("energy"))

        self._update_history(
            rec,
            engine=self.engine_name,
            level=self.level,
            status=status,
            energy=energy,
            xyz_path=output_xyz,
            log_path=output_log,
            backend_result=backend_result,
            start_wall=start_wall,
            end_wall=end_wall,
            start_cpu=start_cpu,
            end_cpu=end_cpu,
        )

        rec.energy = energy
        rec.xyz_path = output_xyz if os.path.exists(output_xyz) else None

        if not self.keep_scratch:
            shutil.rmtree(scratch_dir, ignore_errors=True)

        return rec


    
    
    # ---------------------------------------------------------------------
    # Checkpoint Helpers
    # ---------------------------------------------------------------------

    def _checkpoint_path(self, lookup_id):
        cp = os.path.join(self.outputs_dir, "checkpoints")
        os.makedirs(cp, exist_ok=True)
        return os.path.join(cp, f"{lookup_id}.json")

    def _write_checkpoint(self, lookup_id, rec):
        path = self._checkpoint_path(lookup_id)
        ConformerSet([rec]).save(path)
        self.log(f"[INFO] Checkpoint written: {path}")

    def _clear_checkpoints(self):
        cp = os.path.join(self.outputs_dir, "checkpoints")
        if os.path.exists(cp):
            shutil.rmtree(cp)
        os.makedirs(cp, exist_ok=True)

    def _merge_checkpoints(self):
        cp = os.path.join(self.outputs_dir, "checkpoints")
        all_records = []

        for fname in os.listdir(cp):
            if fname.endswith(".json"):
                cs = ConformerSet.load(os.path.join(cp, fname))
                all_records.extend(cs.records)

        return ConformerSet(all_records)



    # ---------------------------------------------------------------------
    # Backend runner
    # ---------------------------------------------------------------------
    def _run_backend(
        self,
        input_xyz,
        output_xyz,
        output_log,
        scratch_dir,
        engine_name,
        engine_spec,
        level,
        max_iter,
        gfn_level,
    ):
        family = engine_spec["family"]

        if family == "orca":
            return self._backend_orca(
                input_xyz=input_xyz,
                output_xyz=output_xyz,
                output_log=output_log,
                scratch_dir=scratch_dir,
                engine_spec=engine_spec,
                level=level,
                max_iter=max_iter,
            )

        if family == "gxtb":
            return self._backend_gxtb(
                input_xyz=input_xyz,
                output_xyz=output_xyz,
                output_log=output_log,
                scratch_dir=scratch_dir,
                level=level,
                max_iter=max_iter,
                gfn_level=gfn_level,
            )
        if family == "xtb":
            return self._backend_xtb(
                input_xyz=input_xyz,
                output_xyz=output_xyz,
                output_log=output_log,
                scratch_dir=scratch_dir,
                level=level,
                max_iter=max_iter,
                gfn_level=gfn_level,
            )

        if family == "forcefield":
            return self._backend_forcefield(
                input_xyz=input_xyz,
                output_xyz=output_xyz,
                output_log=output_log,
                scratch_dir=scratch_dir,
                engine_spec=engine_spec,
            )

        self.fail(f"Unknown engine family: {family}")

    # ---------------------------------------------------------------------
    # Backend: Forcefield
    # ---------------------------------------------------------------------
    def _backend_forcefield(self, input_xyz, output_xyz, output_log, scratch_dir, engine_spec, **kwargs):
        with open(output_log, "w") as f:
            f.write(f"Forcefield optimisation placeholder ({engine_spec.get('forcefield')})\n")

        shutil.copy(input_xyz, output_xyz)

        return {
            "run_status": "ok",
            "version": "unknown",
            "command": [],
            "command_str": "",
        }

    # ---------------------------------------------------------------------
    # Backend: ORCA (final version)
    # ---------------------------------------------------------------------
    def _backend_orca(
        self,
        input_xyz,
        output_xyz,
        output_log,
        scratch_dir,
        engine_spec,
        level,
        max_iter,
        **kwargs,
    ):
        lookup = os.path.splitext(os.path.basename(input_xyz))[0]
        inp_file = os.path.join(scratch_dir, f"{lookup}.inp")

        opt_keywords = {
            "loose": "LooseOpt",
            "normal": "Opt",
            "tight": "TightOpt",
            "vtight": "VeryTightOpt",
        }

        opt = engine_spec.get("opt", True)
        sp = engine_spec.get("sp", False)
        cpcm_mode = engine_spec.get("cpcm", "none")
        alpb = engine_spec.get("alpb", False)
        method = engine_spec.get("method")
        basis = engine_spec.get("basis")
        solvent = engine_spec.get("solvent", "water")

        lines = []
        lines.append("%MaxCore 2000")
        lines.append("")

        # CPCM block
        if cpcm_mode == "default":
            lines.append("%cpcm")
            lines.append("end")
            lines.append("")

        elif cpcm_mode == "custom":
            chem_dir = self.config["constant_files"]["chemistry_dir"]
            cpcm_path = os.path.join(chem_dir, "cpcm_radii.json")
            with open(cpcm_path) as f:
                cpcm_cfg = json.load(f)

            lines.append("%cpcm")
            for Z, r in cpcm_cfg["radii"].items():
                lines.append(f"  radius[{Z}] {r}")
            lines.append(f"  cut_area {cpcm_cfg['cut_area']}")
            lines.append("end")
            lines.append("")

        # Method line
        if alpb:
            opt_kw = opt_keywords[level] if opt and not sp else "SP"
            method_line = f"! XTB2 {opt_kw} ALPB({solvent})"
        else:
            if opt and not sp:
                opt_kw = opt_keywords[level]
            elif sp and not opt:
                opt_kw = "SP"
            else:
                opt_kw = opt_keywords[level]

            if basis:
                method_line = f"! {method} {basis} {opt_kw}"
            else:
                method_line = f"! {method} {opt_kw}"

        lines.append(method_line)
        lines.append("")
        lines.append(f"%base \"{lookup}\"")
        lines.append("")
        lines.append("* xyz 0 1")

        with open(input_xyz) as xyz:
            xyz_lines = xyz.readlines()
            lines.extend(xyz_lines[2:])

        lines.append("*")
        lines.append("")

        with open(inp_file, "w") as f:
            f.write("\n".join(lines))

        cmd = [self.config["orca"]["executable"], inp_file]

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
                run_status = "ok"
            else:
                output_xyz = None
                run_status = "failed"

        except subprocess.CalledProcessError:
            run_status = "failed"
            output_xyz = None

        return {
            "run_status": run_status,
            "version": "unknown",
            "command": cmd,
            "command_str": " ".join(cmd),
        }

    # ---------------------------------------------------------------------
    # Backend: gXTB
    # ---------------------------------------------------------------------
    def _backend_gxtb(
        self,
        input_xyz,
        output_xyz,
        output_log,
        scratch_dir,
        level,
        max_iter,
        gfn_level,
        **kwargs
    ):
        import os
        import shutil
        import subprocess
        import math

        xtb_bin = self.config["xtb"]["executable"]
        gxtb_bin = self.config["gxtb"]["executable"]

        # Prepare scratch
        tmp_xyz = os.path.join(scratch_dir, "molecule.xyz")
        opt_tmp = os.path.join(scratch_dir, "xtbopt.xyz")

        for p in (tmp_xyz, opt_tmp):
            if os.path.exists(p):
                os.remove(p)

        shutil.copy(input_xyz, tmp_xyz)

        # Optimisation level
        opt_flag = {
            "loose": ["--opt", "loose"],
            "normal": ["--opt"],
            "tight": ["--opt", "tight"],
            "vtight": ["--opt", "vtight"],
        }[level]

        # Parallelisation
        max_cores = os.cpu_count() or 1
        ncores = max(1, math.floor(max_cores * 0.8))

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(ncores)

        # Build driver string EXACTLY like your working script
        driver_string = f"{gxtb_bin} -grad -c xtbdriver.xyz"

        # Build command exactly like your working script
        cmd = [
            xtb_bin,
            "molecule.xyz",
            "--driver", driver_string,
            *opt_flag,
            "--iterations", str(max_iter),
        ]

        # Run
        try:
            with open(output_log, "w") as log:
                subprocess.run(
                    cmd,
                    cwd=scratch_dir,
                    stdout=log,
                    stderr=log,
                    text=True,
                    check=False,
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
            "version": "unknown",
            "command": cmd,
            "command_str": " ".join(cmd),
        }



    # ---------------------------------------------------------------------
    # Backend: classic XTB
    # ---------------------------------------------------------------------
    def _backend_xtb(self, input_xyz, output_xyz, output_log, scratch_dir, level, max_iter, gfn_level, **kwargs):
        xtb_bin = self.config["xtb"]["executable"]

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
            xtb_bin,
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
            "version": "unknown",
            "command": cmd,
            "command_str": " ".join(cmd),
        }

    # ---------------------------------------------------------------------
    # Log parsing
    # ---------------------------------------------------------------------
    def _parse_backend_log(self, log_path, engine_spec):
        family = engine_spec["family"]
        parser_cls = PARSERS.get(family)
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
                "timestamp": datetime.utcnow().isoformat() + "Z",
            }
        )

    # ---------------------------------------------------------------------
    # Write outputs
    # ---------------------------------------------------------------------
    def _write_human_summary_csv(self):
        final_set = self._merge_checkpoints()
        records = final_set.records

        summary_csv = os.path.join(self.outputs_dir, "optimisation_summary.csv")
        rows = []

        for r in records:
            last = r.optimisation_history[-1]
            rows.append({
                "lookup_id": r.lookup_id,
                "inchi_key": r.inchi_key,
                "energy": last["energy"],
                "status": last["status"],
                "xyz_path": last["xyz_path"],
                "elapsed_seconds": last["elapsed_seconds"],
            })

        pd.DataFrame(rows).to_csv(summary_csv, index=False)

