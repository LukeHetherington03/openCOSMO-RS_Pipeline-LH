#!/usr/bin/env python3
import subprocess
from pathlib import Path
import tempfile
import os
import re
import sys
import inspect
import shutil


def run_legacy_cosmors(
    mixture_text: str,
    python_src: str,
    cpp_bindings: str,
    driver_script: str,
    timeout: int = 3600,
):
    """
    Runs the COSMO‑RS Python driver and extracts key results.
    """

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Write mixture_inputs.txt
        mixfile = tmpdir / "mixture_inputs.txt"
        mixfile.write_text(mixture_text)

        # ------------------------------------------------------------
        # Build environment for subprocess
        # ------------------------------------------------------------
        env = os.environ.copy()

        env["PYTHONPATH"] = os.pathsep.join([
            python_src,
            cpp_bindings,
            env.get("PYTHONPATH", "")
        ])

        env["LD_LIBRARY_PATH"] = os.pathsep.join([
            cpp_bindings,
            env.get("LD_LIBRARY_PATH", "")
        ])

        # ------------------------------------------------------------
        # Write import validation file (quiet)
        # ------------------------------------------------------------
        validation_file = tmpdir / "import_validation.txt"

        with open(validation_file, "w") as vf:
            vf.write("=== Import Validation ===\n")
            vf.write(f"python_src: {python_src}\n")
            vf.write(f"cpp_bindings: {cpp_bindings}\n\n")

            try:
                sys.path.insert(0, python_src)
                sys.path.insert(0, cpp_bindings)

                from opencosmorspy import Parameterization, COSMORS
                from opencosmorspy.parameterization import openCOSMORS24a
                import openCOSMORS

                vf.write("Imports: SUCCESS\n")
                vf.write(f"Parameterization: {inspect.getfile(Parameterization)}\n")
                vf.write(f"COSMORS: {inspect.getfile(COSMORS)}\n")
                vf.write(f"openCOSMORS24a: {inspect.getfile(openCOSMORS24a)}\n")
                vf.write(f"openCOSMORS: {openCOSMORS.__file__}\n")

            except Exception as e:
                vf.write("Imports: FAILED\n")
                vf.write(str(e) + "\n")

        # ------------------------------------------------------------
        # Run the Python COSMO‑RS driver
        # ------------------------------------------------------------
        cmd = ["python3", driver_script, "mixture_inputs.txt"]

        # Use Popen + communicate so we have a live process handle and can kill
        # it cleanly on timeout.  subprocess.run() raises TimeoutExpired without
        # a .process attribute, making kill() impossible.
        proc = subprocess.Popen(
            cmd,
            cwd=tmpdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
        )
        try:
            raw_out_bytes, raw_err_bytes = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            raw_out_bytes, raw_err_bytes = proc.communicate()
            return {
                "returncode": None,
                "solubility": None,
                "timed_out": True,
                "saturation_iterations": 0,
                "conformer_iterations": 0,
                "time_taken": None,
                "raw_stdout": raw_out_bytes.decode("utf-8", errors="replace"),
                "raw_stderr": raw_err_bytes.decode("utf-8", errors="replace"),
                "import_validation_file": str(validation_file),
            }

        stdout = raw_out_bytes.decode("utf-8", errors="replace")
        stderr = raw_err_bytes.decode("utf-8", errors="replace")

        # ------------------------------------------------------------
        # Extract key results
        # ------------------------------------------------------------
        solubility = None
        m = re.search(r"Mole fraction solubility:\s+([0-9.Ee+-]+)", stdout)
        if m:
            solubility = float(m.group(1))

        sat_iters = len(re.findall(r"Saturation iteration", stdout))
        conf_iters = len(re.findall(r"Iteration\s+\d+", stdout))

        # Extract per-iteration solubility estimates (the "New mole fraction" line
        # that opens each saturation iteration block).
        saturation_mole_fractions = [
            float(v)
            for v in re.findall(r"New mole fraction:\s+([0-9.Ee+-]+)", stdout)
        ]

        time_taken = None
        m = re.search(r"time taken is\s+([0-9:.\-]+)", stdout)
        if m:
            time_taken = m.group(1)

        return {
            "returncode": proc.returncode,
            "solubility": solubility,
            "saturation_iterations": sat_iters,
            "conformer_iterations": conf_iters,
            "saturation_mole_fractions": saturation_mole_fractions,
            "time_taken": time_taken,
            "raw_stdout": stdout,
            "raw_stderr": stderr,
            "import_validation_file": str(validation_file),
        }
