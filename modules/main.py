#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json

from modules.build.request_manager import Request
from modules.execution.runner import PipelineRunner
from modules.execution.queue import QueueManager

from modules.utils.git_version import get_git_version


# ------------------------------------------------------------
# Execution options
# ------------------------------------------------------------
USE_QUEUE = True    # False → run in-process (blocks until complete)
VERBOSE   = False   # True  → print job progress to stdout
                    # (queue mode: start worker with -v instead)


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

    input_csv = os.path.join(base_dir, "input_data", "acr_t2.csv")

    title = "ff_acr_t2"

    pipeline_spec = [
        {"stage": "cleaning", "args": {"input_csv": input_csv, "overwrite_metadata": True}},
        {"stage": "generation", "args": {"engine": "rdkit", "n": 200}},
        {"stage": "pruning",      "args": {"rmsd_threshold": 1.0}},
        {"stage": "pruning", "args": {"n": 10}},
        {"stage": "optimisation", "args": {"engine": "forcefield_uff"}},
        {"stage": "optimisation", "args": {"engine": "forcefield_mmff"}},
        {"stage": "pruning",      "args": {"energy_window": 6}},
        {"stage": "optimisation", "args": {"engine": "orca_opt_cpcm_fast"}},
        #{"stage": "optimisation", "args": {"engine": "orca_opt_final"}},
        {"stage": "optimisation", "args": {"engine": "orca_opt_cpcm_final"}},
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
        import modules.execution.runner as _runner
        _runner._VERBOSE = VERBOSE
        PipelineRunner.run_request(req.request_id, base_dir)
        print("\nPipeline finished.")


if __name__ == "__main__":
    main()
