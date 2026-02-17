#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
from modules.build.request_manager import Request
from modules.build.job_manager import Job
from modules.execution.runner import PipelineRunner


# ------------------------------------------------------------
# USER-EDITABLE CONFIGURATION
# ------------------------------------------------------------

# The request you want to resume
REQUEST_ID = "R-17Feb26-1349-21475-test_abort"   # <-- EDIT THIS

# Path to global config
CONFIG_PATH = "config/paths.json"

RETRY_FAILED = True

# ------------------------------------------------------------
# MAIN LOGIC
# ------------------------------------------------------------
def main():
    # Load global paths config
    with open(CONFIG_PATH) as f:
        paths = json.load(f)

    base_dir = paths["base_dir"]

    print(f"Resuming request: {REQUEST_ID}")
    print(f"Base directory: {base_dir}")

    # Load existing request (does NOT create a new one)
    req = Request.load(base_dir, REQUEST_ID)

    if RETRY_FAILED:
        job = Job.load(req, req.current_job_id)
        job.retry_failed_items()
        print("Moved failed items back into pending_items.")


    # Hand off to the runner — it will resume automatically
    PipelineRunner.run_request(REQUEST_ID, base_dir)

    print("Resume complete.")


if __name__ == "__main__":
    main()
