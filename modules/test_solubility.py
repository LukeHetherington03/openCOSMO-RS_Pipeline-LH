#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from modules.setup.environment_manager import EnvironmentManager
from modules.solubility_engine.wrapper import run_solubility


def main():
    print("\n==============================================")
    print(" COSMO‑RS Solubility Test (Pipeline Wrapper)  ")
    print("==============================================\n")

    # ------------------------------------------------------------
    # 1. Ensure environment paths are loaded (openCOSMORS + opencosmorspy)
    # ------------------------------------------------------------
    env = EnvironmentManager()
    env.add_to_python_path()

    # ------------------------------------------------------------
    # 2. Define your test files
    # ------------------------------------------------------------
    base = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH"

    solute_c000 = f"{base}/pipeline_data/input_data/2-(dimethylamino)ethylacrylate_c000.orcacosmo"
    solute_c001 = f"{base}/pipeline_data/input_data/2-(dimethylamino)ethylacrylate_c001.orcacosmo"
    water_file = f"{base}/CONSTANT_FILES/solvents/water.orcacosmo"

    # For your colleague’s engine, each molecule must be a directory of conformers.
    # So we create temporary “directories” by using the parent folder.
    solute_dir = os.path.dirname(solute_c000)
    solvent_dir = os.path.dirname(water_file)

    # ------------------------------------------------------------
    # 3. Run solubility at 298.15 K
    # ------------------------------------------------------------
    solubility = run_solubility(
        solute_name="DMAEA",                     # arbitrary label
        solute_smiles="CCN(C)CCOC(=C)C",         # approximate SMILES (not critical unless using Myrdal–Yalkowsky)
        solute_Tm="liquid",                      # or numeric melting point if known
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
    # 4. Print result
    # ------------------------------------------------------------
    print("\n==============================================")
    print(" Solubility result")
    print("==============================================")
    print(f"Solubility (mole fraction): {solubility}")
    print("==============================================\n")


if __name__ == "__main__":
    main()
