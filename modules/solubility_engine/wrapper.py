import os
import sys
import tempfile
import io
import contextlib

# Import your colleague's script as a module
# (Assumes runCOSMO_RS_cpp.py is in the same folder)
from . import runCOSMO_RS_cpp


def run_solubility(
    solute_name,
    solute_smiles,
    solute_Tm,
    solvent_name,
    solute_dir,
    solvent_dir,
    T=298.15,
    Gfus_mode="MyrdalYalkowsky",
    Hfus="N/A",
    w=1.0,
    calc_type="mixed_only"
):
    """
    Runs your colleague's COSMO‑RS solubility engine using a generated mixture_inputs.txt.

    Parameters
    ----------
    solute_name : str
        Name of the solute (e.g. "tert-butylmethacrylate")
    solute_smiles : str
        SMILES string for Myrdal–Yalkowsky ΔGfus estimation
    solute_Tm : float or "liquid"
        Melting point in Kelvin, or "liquid"
    solvent_name : str
        Name of the solvent
    solute_dir : str
        Directory containing solute conformers (orcacosmo files)
    solvent_dir : str
        Directory containing solvent conformers
    T : float
        Temperature in Kelvin
    Gfus_mode : str
        "MyrdalYalkowsky" or "N/A"
    Hfus : str
        Enthalpy of fusion (unused for MY)
    w : float
        Successive over‑relaxation coefficient
    calc_type : str
        "mixed_only", "pure_only", or "all"

    Returns
    -------
    float
        Solubility (mole fraction)
    """

    # ------------------------------------------------------------
    # 1. Build mixture_inputs.txt content
    # ------------------------------------------------------------
    mixture_text = f"""calculations {calc_type}
Temperature {T}
saturation {solute_name} {solute_smiles}
meltingtemp {solute_Tm}
Gfus {Gfus_mode}
Hfus {Hfus}
SORcf {w}
#name #mol fraction #path_to_dir #nconf
{solute_name} 0.0 {solute_dir} 1
{solvent_name} 1.0 {solvent_dir} 1
"""

    # ------------------------------------------------------------
    # 2. Write to a temporary file
    # ------------------------------------------------------------
    with tempfile.TemporaryDirectory() as tmpdir:
        mixture_path = os.path.join(tmpdir, "mixture_inputs.txt")
        with open(mixture_path, "w") as f:
            f.write(mixture_text)

        # ------------------------------------------------------------
        # 3. Prepare to run your colleague's script
        # ------------------------------------------------------------
        # Their script expects sys.argv[1] to be the mixture_inputs file
        sys.argv = ["runCOSMO_RS_cpp.py", mixture_path]

        # Capture stdout so we can extract the solubility
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            runCOSMO_RS_cpp.main() if hasattr(runCOSMO_RS_cpp, "main") else None

        output = buffer.getvalue()

    # ------------------------------------------------------------
    # 4. Extract solubility from output
    # ------------------------------------------------------------
    # Your colleague prints solubility lines like:
    # "Solubility X = 3.2e-4"
    sol = None
    for line in output.splitlines():
        if "Solubility" in line:
            try:
                sol = float(line.split()[-1])
            except:
                pass

    if sol is None:
        raise RuntimeError(
            "Solubility could not be extracted from runCOSMO_RS_cpp output.\n"
            f"Output was:\n{output}"
        )

    return sol
