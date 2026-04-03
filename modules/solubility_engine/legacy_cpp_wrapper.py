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

        try:
            proc = subprocess.run(
                cmd,
                cwd=tmpdir,
                capture_output=True,
                text=True,
                env=env,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired as e:
            e.process.kill()
            e.process.wait()
            raw_out = e.stdout or b""
            raw_err = e.stderr or b""
            if isinstance(raw_out, bytes):
                raw_out = raw_out.decode("utf-8", errors="replace")
            if isinstance(raw_err, bytes):
                raw_err = raw_err.decode("utf-8", errors="replace")
            return {
                "returncode": None,
                "solubility": None,
                "timed_out": True,
                "saturation_iterations": 0,
                "conformer_iterations": 0,
                "time_taken": None,
                "raw_stdout": raw_out,
                "raw_stderr": raw_err,
                "import_validation_file": str(validation_file),
            }

        stdout = proc.stdout
        stderr = proc.stderr

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
