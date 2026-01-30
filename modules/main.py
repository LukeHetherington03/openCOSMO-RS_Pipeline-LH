#!/usr/bin/env python3#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### cd openCOSMO-RS_Pipeline-LH/
### python3 -m modules.main

import os
from modules.provenance.request_manager import Request
from modules.execution.pipeline_runner import PipelineRunner
import json

def main():
    config_path = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/config/paths.json"
    with open(config_path) as f:
        config = json.load(f)


    base_dir = os.path.abspath("pipeline_data")

    input_csv = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/pipeline_data/input_data/NTU_AZ_medicinesol.csv"

    # Declarative pipeline specification
    pipeline_spec = [

        {
            "stage": "cleaning",
            "args": {
                "input_csv": input_csv,
                "seed": 42,
                "engine": "rdkit",
            },
        },

        {
            "stage": "generation",
            "args": {
                "num_confs": 200,
                "seed": 42,
                "engine": "rdkit",
            },
        },

        {
            "stage": "pruning",
            "args": {
                "strategy": "top_n",
                "strategy_params": {"n": 10},
            },
        },

        {
            "stage": "optimisation",
            "args": {
                "engine": "gxtb",
                "max_iter": 250,
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
                "engine": "dft_final",
                "max_iter": 250,
            },
        },

        {
            "stage": "orcacosmo",
            "args": {
                "orca_command": "orca",
                "do_optimization": False,
                "calculate_polarizabilities": True,
            },
        },

        {
            "stage": "solubility",
            "args": {
                "opencosmo_binary": "openCOSMORS",
                "solvent_name": "water",
                "solvent_smiles": "O",
                "temperature": 298.15,
                "calculations": "all",
                "SORcf": 1.0,
            },
        },
    ]
                    

    # Global request parameters
    parameters = {
        "title": "NTU_PHARMA",
        "resources": {
            "cpus": 20,
            "memory_gb": 64,
        },
        "config": config
    }

    # Create the Request
    req = Request.create_new(
        base_dir=base_dir,
        dataset="NTU_PHARMA",
        pipeline_spec=pipeline_spec,
        parameters=parameters,
    )

    print(f"Created Request: {req.request_id}")
    print("Launching pipeline...")

    # Run the pipeline (detached mode will background it)
    PipelineRunner.run_request(req.request_id, base_dir)


if __name__ == "__main__":
    main()