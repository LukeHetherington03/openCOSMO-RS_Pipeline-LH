#!/usr/bin/env python3

import json
from modules.provenance.request_manager import Request
from modules.execution.pipeline_runner import PipelineRunner

# ------------------------------------------------------------
# USER INPUTS
# ------------------------------------------------------------

BASE_DIR = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/pipeline_data"
GLOBAL_CONFIG_PATH = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/config/paths.json"

with open(GLOBAL_CONFIG_PATH) as f:
    global_config = json.load(f)

# Summary file from the *unoptimised* stage
SUMMARY_FILE = (
    "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/pipeline_data/requests/"
    "R-22Jan26-1604-45135-acr_t2_test/jobs/"
    "J-22Jan26-1625-27784-pruning/outputs/energies.json"
)

pipeline_spec = [
    {
        "stage": "orcacosmo",
        "args": {
            "summary_file": SUMMARY_FILE
        }
    }
]

# ------------------------------------------------------------
# CREATE CONTINUATION REQUEST
# ------------------------------------------------------------

request = Request.create_new(
    base_dir=BASE_DIR,
    dataset=None,
    pipeline_spec=pipeline_spec,
    parameters={
        "title": "continue_orcacosmo",
        "summary_file": SUMMARY_FILE,
        "config": global_config,
    }
)

print(f"Created continuation request: {request.request_id}")
print("Launching ORCA-COSMO stage under nohup...")

# ------------------------------------------------------------
# RUN THROUGH PIPELINE RUNNER (DETACHED)
# ------------------------------------------------------------

PipelineRunner.run_request(request.request_id, BASE_DIR)
