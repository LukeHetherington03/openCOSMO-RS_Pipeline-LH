#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json

from modules.build.request_manager import Request
from modules.execution.runner import PipelineRunner
from modules.execution.queue import QueueManager

from modules.utils.git_version import get_git_version


# ------------------------------------------------------------
# Choose execution mode
# ------------------------------------------------------------
USE_QUEUE = True   # Set to True to enqueue instead of running directly


def main():
    # ------------------------------------------------------------
    # Load global paths config
    # ------------------------------------------------------------
    config_path = "config/paths.json"
    with open(config_path) as f:
        paths = json.load(f)

    base_dir = paths["base_dir"]

    # ------------------------------------------------------------
    # Use the ENTIRE paths.json as the pipeline config
    # ------------------------------------------------------------
    config = paths
    config["pipeline_version"] = get_git_version()

    input_csv = os.path.join(base_dir, "input_data", "ml_50_1.csv")

    title = "ml_50_1"

    pipeline_spec = [
        {"stage": "cleaning", "args": {"input_csv": input_csv, "overwrite_metadata": False}},
        {"stage": "generation", "args": {"engine": "rdkit", "num_confs": 200}},
        {"stage": "pruning", "args": {"n": 1}},
        #{"stage": "optimisation", "args": {"engine": "xtb_opt_normal"}},
        {"stage": "optimisation", "args": {"engine": "gxtb_opt_normal"}},
        #{"stage": "optimisation", "args": {"engine": "orca_opt_fast"}},
        #{"stage": "optimisation", "args": {"engine": "orca_opt_cpcm_fast"}},
        #{"stage": "optimisation", "args": {"engine": "orca_opt_final"}},
        #{"stage": "optimisation", "args": {"engine": "orca_opt_cpcm_final"}},
        {"stage": "orcacosmo", "args": {}},
        {"stage": "solubility", "args": {}},
    ]


    parameters = {
        "title": title,
        "resources": {"cpus": 20, "memory_gb": 64},
        "config": config,
    }

    # ------------------------------------------------------------
    # Create the Request
    # ------------------------------------------------------------
    req = Request.create_new(
        base_dir=base_dir,
        dataset=title,
        pipeline_spec=pipeline_spec,
        parameters=parameters,
    )

    print(f"Created Request: {req.request_id}")

    # ------------------------------------------------------------
    # EXECUTION MODE
    # ------------------------------------------------------------
    if USE_QUEUE:
        print("Adding to queue...\n")
        QueueManager.init_base_dir(base_dir)
        QueueManager.enqueue(req.request_id)
        print("Done. Worker will process this request.")
    else:
        print("Running pipeline directly...\n")
        PipelineRunner.run_request(req.request_id, base_dir)
        print("\nPipeline finished.")


if __name__ == "__main__":
    main()
