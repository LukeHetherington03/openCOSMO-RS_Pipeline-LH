#!/usr/bin/env python3
# cd ./openCOSMO-RS_Pipeline-LH/
# python3 -m modules.roundrobin

import os
import pprint
from .workflow_manager import WorkflowManager

def main():
    project_root = os.path.dirname(os.path.dirname(__file__))
    base_dir = os.path.join(project_root, "pipeline_data")

    dataset = "acrylates"
    conf_prod = "rdkit"

    # Paths to pruned lookup and conformer xyz directory
    lookup_csv = os.path.join(
        base_dir,
        "4_pruned_conformers",
        dataset,
        conf_prod,
        "topN_step1",
        f"lookup_{dataset}_{conf_prod}_topN_step1.csv"
    )
    input_xyz_dir = os.path.join(base_dir, "3_conformer_xyz", dataset)

    wf = WorkflowManager(base_dir=base_dir)

    # Run round-robin optimisation
    results = wf.run_roundrobin_optimisation(
        lookup_csv=lookup_csv,
        input_xyz_dir=input_xyz_dir,
        opt_engine="dft_fast",   # or "dft", "gxtb", etc.
        top_n=3,                 # number of conformers per molecule
        functional="BP86",       # engine-specific parameters
        basis="def2-SVP"
    )

    print("\nRound-robin optimisation finished. Outputs:\n")
    pprint.pprint(results, indent=2, width=80)

if __name__ == "__main__":
    main()
