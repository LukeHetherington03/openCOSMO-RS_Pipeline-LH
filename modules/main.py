#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json

from modules.build.request_manager import Request
from modules.execution.runner import PipelineRunner


def main():
    # Load config
    config_path = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/config/paths.json"
    with open(config_path) as f:
        config = json.load(f)

    base_dir = os.path.abspath("pipeline_data")

    input_csv= ("/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/pipeline_data/input_data/acrylate_initiator_set_1.csv")

    # Declarative pipeline specification
    pipeline_spec = [
        {"stage": "cleaning", "args": {"input_csv":input_csv,"overwrite_metadata":True}},
        #{"stage": "generation", "args": {"engine": "rdkit", "num_confs": 100}},
        #{"stage": "pruning", "args": {"n": 2}},
        #{"stage": "optimisation", "args": {"engine": "gxtb"}},
        #{"stage": "orcacosmo", "args": {}},
        #{"stage": "solubility", "args": {}},
    ]

    # Global request parameters
    parameters = {
        "title": "overwrite_metadata",
        "resources": {"cpus": 20, "memory_gb": 64},
        "config": config,
    }

    # Create the Request
    req = Request.create_new(
        base_dir=base_dir,
        dataset="overwrite_metadata",
        pipeline_spec=pipeline_spec,
        parameters=parameters,
    )

    print(f"Created Request: {req.request_id}")
    print("Running pipeline directly (no queue, no worker)...\n")

    # ------------------------------------------------------------
    # RUN THE PIPELINE DIRECTLY
    # ------------------------------------------------------------
    PipelineRunner.run_request(req.request_id, base_dir)

    print("\nPipeline finished.")


if __name__ == "__main__":
    main()