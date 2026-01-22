"""
cosmors_wrapper.py

A clean, pipeline‑friendly wrapper around the legacy runCOSMO_RS_cpp.py driver.

This module is responsible ONLY for:
    - generating mixture_inputs.txt
    - preparing a job‑local mixture directory
    - temporarily chdir'ing into that directory for import
    - importing runCOSMO_RS_cpp (which reads mixture_inputs.txt at import time)
    - restoring the working directory
    - calling runCOSMO_RS_cpp.run() with absolute paths

It does NOT:
    - modify COSMO‑RS physics
    - change iteration behaviour
    - alter numerical settings
    - touch conformer logic
    - interpret results

This keeps the pipeline clean and isolates all legacy behaviour.
"""

from pathlib import Path
import os
import sys
import json
import io
import contextlib


class COSMORSWrapper:
    """
    Wrapper class for running COSMO‑RS solubility calculations using the
    legacy runCOSMO_RS_cpp.py driver.

    Parameters
    ----------
    config_path : str or Path
        Path to a minimal JSON config containing mixture‑file defaults.
        (temperature, calculation_type, saturation settings, etc.)

    opencosmo_paths : dict
        The "opencosmo" block from your global config.json.
        Must contain:
            - python_src
            - cpp_bindings
            - python_home (optional)
    """

    def __init__(self, config_path: Path, opencosmo_paths: dict):
        self.config = self._load_config(config_path)
        self.opencosmo_paths = opencosmo_paths

        # Add Python + C++ bindings to sys.path
        sys.path.insert(0, opencosmo_paths["python_src"])
        sys.path.insert(0, opencosmo_paths["cpp_bindings"])

        # Ensure LD_LIBRARY_PATH contains the C++ bindings
        ld = os.environ.get("LD_LIBRARY_PATH", "")
        if opencosmo_paths["cpp_bindings"] not in ld.split(":"):
            os.environ["LD_LIBRARY_PATH"] = (
                opencosmo_paths["cpp_bindings"] + (":" + ld if ld else "")
            )

    # ------------------------------------------------------------------
    # CONFIG LOADING
    # ------------------------------------------------------------------

    def _load_config(self, path: Path) -> dict:
        """Load the minimal solubility config JSON."""
        with open(path, "r") as f:
            return json.load(f)

    # ------------------------------------------------------------------
    # MIXTURE FILE GENERATION
    # ------------------------------------------------------------------

    def _write_mixture_file(
        self,
        mixture_path: Path,
        solute_name: str,
        solute_smiles: str,
        solute_Tm: float,
        solute_dir: Path,
        solvent_name: str,
        solvent_dir: Path,
        n_solute: int,
        n_solvent: int,
    ):
        """
        Write mixture_inputs.txt using the minimal config and the job's
        solute/solvent directories.
        """

        calc_type = self.config["calculation_type"]
        T = self.config["temperature"]

        sat_cfg = self.config["saturation"]
        saturation_line = (
            f"{solute_name} {solute_smiles}" if sat_cfg["enabled"] else "no"
        )

        lines = [
            f"calculations {calc_type}",
            f"temperature {T}",
            f"saturation {saturation_line}",
            f"meltingtemp {solute_Tm}",
            f"Gfus {sat_cfg['Gfus_mode']}",
            f"Hfus {sat_cfg['Hfus']}",
            f"SORcf {sat_cfg['SORcf']}",
            "# name molfrac path_to_dir nconf multiplicities",
        ]

        # Solute (molfrac = 0.0)
        lines.append(
            f"{solute_name}\t0.0\t{solute_dir}\t{n_solute}\t"
            + "\t".join(["1"] * n_solute)
        )

        # Solvent (molfrac = 1.0)
        lines.append(
            f"{solvent_name}\t1.0\t{solvent_dir}\t{n_solvent}\t"
            + "\t".join(["1"] * n_solvent)
        )

        mixture_path.write_text("\n".join(lines))

    # ------------------------------------------------------------------
    # MAIN ENTRY POINT
    # ------------------------------------------------------------------

    def run(
        self,
        job_dir: Path,
        solute_name: str,
        solute_smiles: str,
        solute_Tm: float,
        solute_dir: Path,
        solvent_name: str,
        solvent_dir: Path,
    ):
        """
        Run a COSMO‑RS solubility calculation.

        Parameters
        ----------
        job_dir : Path
            The job directory for this solubility calculation.
            A mixture directory will be created inside it.

        solute_name : str
            InChIKey or identifier for the solute.

        solute_smiles : str
            SMILES string for the solute.

        solute_Tm : float
            Melting point (K) or "liquid".

        solute_dir : Path
            Directory containing solute conformer .orcacosmo files.

        solvent_name : str
            Name of the solvent.

        solvent_dir : Path
            Directory containing solvent conformer .orcacosmo files.

        Returns
        -------
        result : any
            Whatever runCOSMO_RS_cpp.run() returns.
        """

        # --------------------------------------------------------------
        # 1. Prepare job-local mixture directory
        # --------------------------------------------------------------
        mixture_dir = job_dir / "mixture"
        results_dir = job_dir / "results"

        mixture_dir.mkdir(parents=True, exist_ok=True)
        results_dir.mkdir(parents=True, exist_ok=True)

        mixture_file = mixture_dir / "mixture_inputs.txt"

        # Count conformers
        solute_confs = sorted(solute_dir.glob("*.orcacosmo"))
        solvent_confs = sorted(solvent_dir.glob("*.orcacosmo"))

        n_solute = len(solute_confs)
        n_solvent = len(solvent_confs)

        # --------------------------------------------------------------
        # 2. Write mixture_inputs.txt
        # --------------------------------------------------------------
        self._write_mixture_file(
            mixture_file,
            solute_name,
            solute_smiles,
            solute_Tm,
            solute_dir,
            solvent_name,
            solvent_dir,
            n_solute,
            n_solvent,
        )

        # --------------------------------------------------------------
        # 3. TEMPORARY CHDIR FOR IMPORT
        # --------------------------------------------------------------
        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()

        orig_cwd = os.getcwd()
        os.chdir(mixture_dir)

        try:
            with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
                import modules.solubility_engine.runCOSMO_RS_cpp as runCOSMO_RS_cpp
        finally:
            os.chdir(orig_cwd)

        import_output = stdout_buf.getvalue() + "\n" + stderr_buf.getvalue()


        # --------------------------------------------------------------
        # 4. Capture raw COSMO‑RS output
        # --------------------------------------------------------------
        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()

        with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(stderr_buf):
            result = runCOSMO_RS_cpp.run(
                mixture_file=str(mixture_file),
                results_dir=str(results_dir),
            )

        run_output = stdout_buf.getvalue() + "\n" + stderr_buf.getvalue()
        raw_output = import_output + "\n" + run_output

        # --------------------------------------------------------------
        # 5. Return structured result
        # --------------------------------------------------------------
        return {
            "result": result,
            "raw_output": raw_output,
            "n_solute": n_solute,
            "n_solvent": n_solvent,
        }
