#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### cd openCOSMO-RS_Pipeline-LH/
### python3 -m modules.main

import os
from modules.provenance.request_manager import Request
from modules.execution.pipeline_runner import PipelineRunner


def main():
    base_dir = os.path.abspath("pipeline_data")

    input_csv = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/pipeline_data/input_data/acr_t2.csv"

    # Declarative pipeline specification
    pipeline_spec = [

        {
            "stage": "cleaning",
            "args": {
                "input_csv": input_csv,
                "num_confs": 25,
                "seed": 42,
                "engine": "rdkit",
            },
        },
        # --------------------------------------------------------
        # 1. RDKit conformer generation (25)
        # --------------------------------------------------------
        {
            "stage": "generation",
            "args": {
                "num_confs": 25,
                "seed": 42,
                "engine": "rdkit",
            },
        },

        # --------------------------------------------------------
        # 2. Prune to top 3
        # --------------------------------------------------------
        {
            "stage": "pruning",
            "args": {
                "strategy": "top_n",
                "strategy_params": {"n": 3},
            },
        },

        # --------------------------------------------------------
        # 3. gXTB optimisation
        # --------------------------------------------------------
        {
            "stage": "optimisation",
            "args": {
                "engine": "gxtb",
                "max_iter": 250,
            },
        },

        # --------------------------------------------------------
        # 4. Prune to top 1
        # --------------------------------------------------------
        {
            "stage": "pruning",
            "args": {
                "strategy": "top_n",
                "strategy_params": {"n": 1},
            },
        },

        # --------------------------------------------------------
        # 5. DFT optimisation (ORCA)
        # --------------------------------------------------------
#        {
#            "stage": "optimisation",
#             "args": {
#                 "engine": "orca_final",
#                 "max_iter": 250,
#             },
#         },

        # --------------------------------------------------------
        # 6. ORCA COSMO calculation
        # --------------------------------------------------------
        {
            "stage": "orcacosmo",
            "args": {
                "orca_command": "orca",
                "do_optimization": False,  # Already optimized with gXTB
                "calculate_polarizabilities": True,  # Can set False to skip if needed
            },
        },
    ]

    # Global request parameters
    parameters = {
        "title": "acr_t5_full_pipeline",
        "detached": False,  # run in background
        "resources": {
            "cpus": 20,
            "memory_gb": 64,
        },
    }

    # Create the Request
    req = Request.create_new(
        base_dir=base_dir,
        dataset="acr_t5",
        pipeline_spec=pipeline_spec,
        parameters=parameters,
    )

    print(f"Created Request: {req.request_id}")
    print("Launching pipeline...")

    # Run the pipeline (detached mode will background it)
    PipelineRunner.run_request(req.request_id, base_dir)


if __name__ == "__main__":
    main()
