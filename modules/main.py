#!/usr/bin/env python3
# cd ./openCOSMO-RS_Pipeline-LH/
# python3 -m modules.main

import os
import pprint
from .workflow_manager import WorkflowManager

def main():
    project_root = os.path.dirname(os.path.dirname(__file__))
    base_dir = os.path.join(project_root, "pipeline_data")

    dataset = "acrylates"
    conf_prod = "rdkit"

    wf = WorkflowManager(base_dir=base_dir)

    # Define argument dictionaries for each stage
    cleaning_args = {
        "dataset": dataset,
        "cleaned_filename": f"{dataset}/{dataset}_clean.csv",
    }

    generation_args = {
        "clean_csv": f"{dataset}/{dataset}_clean.csv",
        "dataset": dataset,
        "engine": conf_prod,
        "n_conformers": 100,
        "seed": 42,
        "append": False,
    }

    pruning_args = {
        "dataset": "acrylates",
        "engine": "rdkit",
        "method": "topN",
        "top_n": 20,
    }

    optimisation_args = {
        "dataset": dataset,
        "engine": "dft",
        "xyz_dir": os.path.join(base_dir, "3_conformer_xyz", dataset, conf_prod),
        "lookup_csv": os.path.join(base_dir,
                                "4_pruned_conformers",
                                dataset,
                                conf_prod,
                                "topN_20_step1",
                                f"lookup_{dataset}_{conf_prod}_topN_20_step1.csv")
    }

    opt_args = {
        "dataset": "acrylates",
        "engine": "dft",                 # optimisation engine
        "conf_prod": "rdkit",             # conformer generation engine
        "lookup_csv": (
            "pipeline_data/4_pruned_conformers/"
            "acrylates/rdkit/topN_20_step1/"
            "lookup_acrylates_rdkit_topN_20_step2.csv"
        ),
        "xyz_dir": (
            "pipeline_data/3_conformer_xyz/"
            "acrylates/rdkit"
        ),
        "max_iter": 250
    }

    cosmo_args = {
        "dataset": dataset,
        "engine": "orca",
    }

    # Build the pipeline sequence explicitly in the order you want
    pipeline_sequence = [
        #{"stage": "cleaning", "args": cleaning_args},
        #{"stage": "generation", "args": generation_args},
        #{"stage": "pruning", "args": pruning_args},
        {"stage": "optimisation", "args": opt_args},
        #{"stage": "pruning", "args": pruning_args2},
        #{"stage": "optimisation", "args": optimisation_args2},
        #{"stage": "cosmo", "args": cosmo_args},
    ]

    # Run pipeline with the ordered sequence
    results = wf.run_pipeline(pipeline_sequence)

    print("\nPipeline finished. Outputs:\n")
    for stage, out in results.items():
        print(f"--- {stage.upper()} ---")
        pprint.pprint(out, indent=2, width=80)  # pretty-print dictionary
        print()  # blank line between stages

if __name__ == "__main__":
    main()
