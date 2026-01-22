#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

# ============================================================
# 0. Generate mixture_inputs.txt BEFORE importing the engine
# ============================================================

base = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH"

# Updated solute directory (your new folder)
solute_dir  = f"{base}/pipeline_data/input_data/solute_dir"

# Updated solvent directory (your new folder)
solvent_dir = f"{base}/CONSTANT_FILES/solvents/water_dir"

mixture_path = "mixture_inputs.txt"

with open(mixture_path, "w") as f:
    f.write("calculations mixed_only\n")
    f.write("temperature 298.15\n")
    f.write("saturation DMAEA CCN(C)CCOC(=C)C\n")
    f.write("meltingtemp liquid\n")
    f.write("Gfus MyrdalYalkowsky\n")
    f.write("Hfus N/A\n")
    f.write("SORcf 1\n")
    f.write("# name molfrac path_to_dir nconf multiplicities\n")
    f.write(f"DMAEA 0.0 {solute_dir} 2 1 1\n")
    f.write(f"water 1.0 {solvent_dir} 1 1\n")

# ============================================================
# 1. Ensure openCOSMORS C++ extension is importable
# ============================================================

cpp_bindings = "/home/lunet/cglh4/resources/openCOSMO-RS_cpp/bindings"
if cpp_bindings not in sys.path:
    sys.path.insert(0, cpp_bindings)

# ============================================================
# 2. Import environment + wrapper (safe now that file exists)
# ============================================================

from modules.setup.environment_manager import EnvironmentManager
from modules.solubility_engine.wrapper import run_solubility


def main():
    print("\n==============================================")
    print(" COSMO‑RS Solubility Test (Pipeline Wrapper)  ")
    print("==============================================\n")

    # ------------------------------------------------------------
    # Load environment paths (Python + C++ bindings)
    # ------------------------------------------------------------
    env = EnvironmentManager()
    env.add_to_python_path()

    # ------------------------------------------------------------
    # Run solubility at 298.15 K
    # ------------------------------------------------------------
    solubility = run_solubility(
        solute_name="DMAEA",
        solute_smiles="CCN(C)CCOC(=C)C",
        solute_Tm="liquid",
        solvent_name="water",
        solute_dir=solute_dir,
        solvent_dir=solvent_dir,
        T=298.15,
        Gfus_mode="MyrdalYalkowsky",
        Hfus="N/A",
        w=1.0,
        calc_type="mixed_only"
    )

    # ------------------------------------------------------------
    # Print result
    # ------------------------------------------------------------
    print("\n==============================================")
    print(" Solubility result")
    print("==============================================")
    print(f"Solubility (mole fraction): {solubility}")
    print("==============================================\n")


if __name__ == "__main__":
    main()
