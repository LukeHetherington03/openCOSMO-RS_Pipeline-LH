#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Runs gXTB optimisation → orcacosmo → solubility starting from a custom
# pre-built energies.json (pre-DFT MMFF94 geometries, 1 conformer/molecule).
#
# Used to compare: DFT_fast-optimised geometry + COSMO
#              vs: MMFF94 geometry + gXTB opt + COSMO
# for the same conformers (min or max by DFT energy ranking).
#
# Set STAGE_INPUT and TITLE below, then run twice (once for min, once for max).

import os
import json

from modules.build.request_manager import Request
from modules.execution.runner import PipelineRunner
from modules.execution.queue import QueueManager


# ------------------------------------------------------------
# USER-EDITABLE CONFIGURATION
# ------------------------------------------------------------

# Path to the custom pre-DFT energies.json (1 conformer per molecule)
# Options:
#   pipeline_data/custom_inputs/acr19_preDFT_min.json   <- DFT-min conformers, MMFF94 geom
#   pipeline_data/custom_inputs/acr19_preDFT_max.json   <- DFT-max conformers, MMFF94 geom
STAGE_INPUT = "pipeline_data/custom_inputs/acr19_preDFT_max.json"

TITLE = "gxtb_preDFT_max"   # change to "gxtb_preDFT_max" for the max run

USE_QUEUE = True

PIPELINE_SPEC = [
    {"stage": "optimisation", "args": {"engine": "gxtb_opt_normal",
                                       "stage_input": STAGE_INPUT}},
    {"stage": "orcacosmo",    "args": {}},
    {"stage": "solubility",   "args": {}},
]


# ------------------------------------------------------------
# MAIN LOGIC
# ------------------------------------------------------------
def main():
    config_path = "config/paths.json"
    with open(config_path) as f:
        paths = json.load(f)

    base_dir = paths["base_dir"]

    # Resolve STAGE_INPUT to absolute path
    abs_stage_input = os.path.join(base_dir, "..", STAGE_INPUT)
    abs_stage_input = os.path.normpath(abs_stage_input)
    if not os.path.exists(abs_stage_input):
        # Try as already-absolute
        abs_stage_input = os.path.abspath(STAGE_INPUT)
    if not os.path.exists(abs_stage_input):
        raise FileNotFoundError(f"STAGE_INPUT not found: {STAGE_INPUT}")

    # Patch the stage_input into the first pipeline step
    PIPELINE_SPEC[0]["args"]["stage_input"] = abs_stage_input

    parameters = {
        "title": TITLE,
        "resources": {"cpus": 20, "memory_gb": 64},
        "config": paths,
    }

    req = Request.create_new(
        base_dir=base_dir,
        dataset=TITLE,
        pipeline_spec=PIPELINE_SPEC,
        parameters=parameters,
    )

    print(f"Created Request: {req.request_id}")
    print(f"Stage input:     {abs_stage_input}")
    print(f"Title:           {TITLE}")

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
