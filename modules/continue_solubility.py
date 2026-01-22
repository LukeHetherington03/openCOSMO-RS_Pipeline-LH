#!/usr/bin/env python3

import os
from modules.provenance.request_manager import Request
from modules.provenance.job_manager import Job
import json

# ------------------------------------------------------------
# USER INPUTS
# ------------------------------------------------------------

BASE_DIR = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/pipeline_data"
GLOBAL_CONFIG_PATH = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/config/paths.json"
with open(GLOBAL_CONFIG_PATH) as f:
    config = json.load(f)


SUMMARY_FILE = (
    "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/pipeline_data/requests/"
    "R-22Jan26-1604-45135-acr_t2_test/jobs/"
    "J-22Jan26-1625-27784-orcacosmo/outputs/orcacosmo_summary.json"
)

# The stage you want to run
TARGET_STAGE = "solubility"

# ------------------------------------------------------------
# PIPELINE SPEC FOR CONTINUATION
# ------------------------------------------------------------
# A pipeline spec is a list of dicts: [{"stage": "...", "args": {...}}, ...]

pipeline_spec = [
    {
        "stage": TARGET_STAGE,
        "args": {
            "summary_file": SUMMARY_FILE
        }
    }
]

# ------------------------------------------------------------
# CREATE A NEW REQUEST THAT STARTS AT SOLUBILITY
# ------------------------------------------------------------

with open(GLOBAL_CONFIG_PATH) as f:
    global_config = json.load(f)

request = Request.create_new(
    base_dir=BASE_DIR,
    dataset=None,
    pipeline_spec=pipeline_spec,
    parameters={
        "title": "continue_solubility",
        "summary_file": SUMMARY_FILE,
        "config": global_config,   # <-- REQUIRED
        "detached": True,
    }
)


print(f"Created continuation request: {request.request_id}")
print("Now running the solubility stage...")

# ------------------------------------------------------------
# RUN THE FIRST (AND ONLY) JOB
# ------------------------------------------------------------

job_id = request.jobs[0]
job = Job.load(request=request, job_id=job_id)

job.run()

print(f"Solubility job complete: {job_id}")
