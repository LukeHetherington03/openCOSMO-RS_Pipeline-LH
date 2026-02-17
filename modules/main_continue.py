#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json

from modules.build.request_manager import Request
from modules.execution.runner import PipelineRunner
from modules.execution.queue import QueueManager


# ------------------------------------------------------------
# USER-EDITABLE CONFIGURATION
# ------------------------------------------------------------

# Parent request + job to continue from
PARENT_REQUEST_ID = "R-17Feb26-1204-53419-orca_fast_t2"     # <-- EDIT THIS
PARENT_JOB_ID     = "J-17Feb26-1213-10503-optimisation"     # <-- EDIT THIS

# Title for the new continuation request
NEW_TITLE = "continued_run"

# Whether to enqueue instead of running directly
USE_QUEUE = False

# Define the new pipeline spec you want to run
PIPELINE_SPEC = [
    {"stage": "orcacosmo", "args": {}},
    {"stage": "solubility", "args": {}},
]


# ------------------------------------------------------------
# MAIN LOGIC
# ------------------------------------------------------------
def main():
    # ------------------------------------------------------------
    # Load global paths config
    # ------------------------------------------------------------
    config_path = "config/paths.json"
    with open(config_path) as f:
        paths = json.load(f)

    base_dir = paths["base_dir"]

    # ------------------------------------------------------------
    # Parameters for the new request
    # ------------------------------------------------------------
    parameters = {
        "title": NEW_TITLE,
        "config": paths,
        "resources": {"cpus": 20, "memory_gb": 64},
    }

    # ------------------------------------------------------------
    # Create continuation request
    # ------------------------------------------------------------
    req = Request.continue_from(
        base_dir=base_dir,
        parent_request_id=PARENT_REQUEST_ID,
        parent_job_id=PARENT_JOB_ID,
        dataset=NEW_TITLE,
        pipeline_spec=PIPELINE_SPEC,
        parameters=parameters,
    )


    print(f"Created continuation Request: {req.request_id}")
    print(f"Continuing from: {PARENT_REQUEST_ID}/{PARENT_JOB_ID}")

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
