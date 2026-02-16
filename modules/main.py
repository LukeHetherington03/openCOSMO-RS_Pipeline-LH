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

    input_csv= ("/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/pipeline_data/input_data/acr_t1.csv")

    # Declarative pipeline specification
    pipeline_spec = [
        {"stage": "cleaning", "args": {"input_csv":input_csv,"overwrite_metadata":False}},
        {"stage": "generation", "args": {"engine": "openbabel", "num_confs": 10}},
        {"stage": "pruning", "args": {"n": 2}},
        {"stage": "optimisation", "args": {"engine": "xtb_opt_normal"}},
        {"stage": "orcacosmo", "args": {}},
        {"stage": "solubility", "args": {}},
    ]

    # Global request parameters
    parameters = {
        "title": "new_run",
        "resources": {"cpus": 20, "memory_gb": 64},
        "config": config,
    }

    # Create the Request
    req = Request.create_new(
        base_dir=base_dir,
        dataset="new_run",
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