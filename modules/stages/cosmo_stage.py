#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import shutil
import subprocess
from dataclasses import dataclass

from modules.parsers.xtb_cosmo_parser import XtbCosmoParser
from modules.parsers.orca_cpcm_parser import OrcaCpcmParser


# ============================================================
#  Result Dataclass
# ============================================================

@dataclass
class CosmoSurfaceResult:
    lookup_id: str
    xyz_input: str
    xyz_output: str
    cosmo_file: str
    cpcm_file: str
    log_file: str
    parsed_json: str
    xtb_cosmo_data: object
    cpcm_data: dict


# ============================================================
#  Exception
# ============================================================

class CosmoSurfaceError(Exception):
    pass


# ============================================================
#  COSMO Surface + CPCM Stage
# ============================================================

class CosmoStage:
    """
    Runs XTB COSMO + ORCA CPCM for each lookup_id in a Job.

    Job model:
      - Job auto-initialises items from input_folder
      - job.pending_items contains lookup_ids to process
      - Stages must NOT initialise items
      - Stages must NOT modify job.items directly
    """

    def __init__(self, job):
        from modules.provenance.job_manager import Job
        if not isinstance(job, Job):
            raise TypeError("CosmoSurfaceGen must be instantiated with a Job object")

        self.job = job

        # Commands (may be overridden in job.parameters)
        self.xtb_command = job.parameters.get("xtb_command", "xtb")
        self.orca_command = job.parameters.get("orca_command", "orca")

    # ------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------
    def run(self):
        job = self.job

        try:
            job.write_pipeline_log("Starting COSMO surface + CPCM stage")
            self._validate_xtb()

            input_folder = job.parameters["input_folder"]
            results = []

            # --------------------------------------------------------
            # Loop over pending lookup_ids
            # --------------------------------------------------------
            for lookup_id in job.pending_items.copy():

                xyz_path = os.path.join(input_folder, f"{lookup_id}.xyz")

                if not os.path.exists(xyz_path):
                    job.write_pipeline_log(
                        f"Skipping {lookup_id}: missing XYZ file ({xyz_path})"
                    )
                    job.update_progress(lookup_id)
                    continue

                job.write_pipeline_log(f"Processing {lookup_id}")

                result = self._run_single(lookup_id, xyz_path)
                results.append(result)

                job.update_progress(lookup_id)

            job.mark_complete()
            job.write_pipeline_log("COSMO surface + CPCM stage completed successfully.")

            return results

        except Exception as e:
            job.fail(str(e))
            raise

    # ------------------------------------------------------------
    # Run XTB COSMO + ORCA CPCM for a single lookup_id
    # ------------------------------------------------------------
    def _run_single(self, lookup_id: str, xyz_path: str) -> CosmoSurfaceResult:
        job = self.job
        root = job.output_dir

        # Directories
        input_dir = os.path.join(root, "input")
        output_raw = os.path.join(root, "output_raw")
        output_parsed = os.path.join(root, "output_parsed")
        tmp_exec = os.path.join(root, "tmp_exec")

        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(output_raw, exist_ok=True)
        os.makedirs(output_parsed, exist_ok=True)
        os.makedirs(tmp_exec, exist_ok=True)

        # ------------------------------------------------------------
        # Copy XYZ into input + tmp_exec
        # ------------------------------------------------------------
        input_xyz_copy = os.path.join(input_dir, f"{lookup_id}.xyz")
        shutil.copy(xyz_path, input_xyz_copy)

        tmp_xyz = os.path.join(tmp_exec, f"{lookup_id}.xyz")
        shutil.copy(xyz_path, tmp_xyz)

        # ------------------------------------------------------------
        # Run XTB COSMO
        # ------------------------------------------------------------
        raw_log = os.path.join(output_raw, f"{lookup_id}_cosmo.log")

        cmd = [
            self.xtb_command,
            tmp_xyz,
            "--cosmo",
            job.parameters.get("solvent", "water"),
            "--norestart"
        ]

        with open(raw_log, "w") as logf:
            process = subprocess.run(
                cmd,
                cwd=tmp_exec,
                stdout=logf,
                stderr=logf,
                text=True
            )

        if process.returncode != 0:
            raise CosmoSurfaceError(
                f"XTB COSMO failed for {lookup_id}. See log: {raw_log}"
            )

        # Collect raw outputs
        xtb_cosmo = os.path.join(tmp_exec, "xtb.cosmo")
        if not os.path.exists(xtb_cosmo):
            raise CosmoSurfaceError(f"xtb.cosmo missing for {lookup_id}")

        final_xyz = os.path.join(output_raw, f"{lookup_id}.xyz")
        final_cosmo = os.path.join(output_raw, f"{lookup_id}.cosmo")

        shutil.copy(tmp_xyz, final_xyz)
        shutil.copy(xtb_cosmo, final_cosmo)

        # ------------------------------------------------------------
        # Parse XTB COSMO
        # ------------------------------------------------------------
        xtb_data = XtbCosmoParser(final_cosmo).get_cosmo_data()
        xtb_json = xtb_data.to_dict()

        # ------------------------------------------------------------
        # Run ORCA CPCM
        # ------------------------------------------------------------
        cpcm_file = self._run_orca_cpcm(lookup_id, tmp_xyz)

        # ------------------------------------------------------------
        # Parse ORCA CPCM
        # ------------------------------------------------------------
        cpcm_data = OrcaCpcmParser(cpcm_file).parse()

        # ------------------------------------------------------------
        # Merge JSON
        # ------------------------------------------------------------
        merged = {**xtb_json, **cpcm_data}
        merged["lookup_id"] = lookup_id
        merged["inchi_key"] = lookup_id.split("_conf")[0]

        parsed_json = os.path.join(output_parsed, f"{lookup_id}_cosmo_full.json")
        with open(parsed_json, "w") as f:
            json.dump(merged, f, indent=2)

        # ------------------------------------------------------------
        # Clean up tmp_exec
        # ------------------------------------------------------------
        shutil.rmtree(tmp_exec, ignore_errors=True)

        return CosmoSurfaceResult(
            lookup_id=lookup_id,
            xyz_input=xyz_path,
            xyz_output=final_xyz,
            cosmo_file=final_cosmo,
            cpcm_file=cpcm_file,
            log_file=raw_log,
            parsed_json=parsed_json,
            xtb_cosmo_data=xtb_data,
            cpcm_data=cpcm_data,
        )

    # ------------------------------------------------------------
    # Run ORCA CPCM
    # ------------------------------------------------------------
    def _run_orca_cpcm(self, lookup_id: str, xyz_path: str) -> str:
        job = self.job
        root = job.output_dir

        cpcm_dir = os.path.join(root, "cpcm_tmp")
        os.makedirs(cpcm_dir, exist_ok=True)

        inp_path = os.path.join(cpcm_dir, f"{lookup_id}.inp")
        log_path = os.path.join(cpcm_dir, "log_output.dat")

        # Write ORCA CPCM input
        with open(inp_path, "w") as f:
            f.write("%MaxCore 2000\n\n")
            f.write("%cpcm\nend\n\n")
            f.write("! CPCM BP86 def2-SVP SP\n\n")
            f.write(f"* xyzfile 0 1 {xyz_path}\n")

        # Run ORCA
        cmd = [self.orca_command, inp_path]
        with open(log_path, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=out, cwd=cpcm_dir, check=True)

        # ORCA writes <lookup_id>.cpcm
        cpcm_file = os.path.join(cpcm_dir, f"{lookup_id}.cpcm")
        if not os.path.exists(cpcm_file):
            raise CosmoSurfaceError(f"CPCM file not generated for {lookup_id}")

        return cpcm_file

    # ------------------------------------------------------------
    # XTB validation
    # ------------------------------------------------------------
    def _validate_xtb(self):
        cmd = [self.xtb_command, "--version"]

        try:
            subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True
            )
        except Exception:
            raise CosmoSurfaceError(
                f"XTB command not found or failed: '{self.xtb_command}'"
            )
