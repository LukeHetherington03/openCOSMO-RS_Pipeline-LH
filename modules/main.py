#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### cd openCOSMO-RS_Pipeline-LH/
### python3 -m modules.main

import os
from modules.provenance.request_manager import Request
from modules.execution.pipeline_runner import PipelineRunner


def main():
    base_dir = os.path.abspath("pipeline_data")
    input_csv = os.path.join(base_dir, "input_data", "test_3_mols.csv")


    # Declarative pipeline specification
    pipeline_spec = [
        {
            "stage": "generation",
            "args": {
                "input_csv": input_csv,
                "num_confs": 30,
                "seed": 42,
            },
        },
        {
            "stage": "pruning",
            "args": {
                "strategy": "top_n",
                "strategy_params": {"n": 3},
            },
        },
        {
            "stage": "optimisation",
            "args": {
                "engine": "gxtb",
                "max_iter": 250,
            },
        },
    ]

    # Global request parameters
    parameters = {
        "title": "test_3_mols_detached",
        "detached": True,  # run in background
        "resources": {
            "cpus": 19,
            "memory_gb": 32,
        },
    }

    # Create the Request
    req = Request.create_new(
        base_dir=base_dir,
        dataset="test_3_mols",
        pipeline_spec=pipeline_spec,
        parameters=parameters,
    )

    print(f"Created Request: {req.request_id}")
    print("Launching pipeline...")

    # Run the pipeline (detached mode will background it)
    PipelineRunner.run_request(req.request_id, base_dir)


if __name__ == "__main__":
    main()
