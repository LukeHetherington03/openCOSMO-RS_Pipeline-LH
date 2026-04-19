#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Investigation of the 6 RMSD-surviving conformers of PFPUZMSQZJFLBK-UHFFFAOYSA-N.

Continues from the RMSD-pruning job of request R-10Apr26-1139-43214-pruning_high_dft_fast,
isolates this one molecule, and forks into 6 separate pipelines — one per conformer
(indices 0–5 in energy-sorted order post-pruning).

Each job runs:
  pruning   → select_conf (isolate molecule) + conf_select [i] (pick one conformer)
  optimisation → gxtb_opt_normal
  orcacosmo
  solubility

Set USE_QUEUE = True (after restarting the worker) to enqueue all 6.
"""

import os
import json

from modules.build.request_manager import Request
from modules.execution.queue import QueueManager

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

PARENT_REQUEST_ID = "R-10Apr26-1139-43214-pruning_high_dft_fast"
PARENT_JOB_ID     = "J-10Apr26-1140-13760-pruning"

TARGET_INCHI_KEY  = "PFPUZMSQZJFLBK-UHFFFAOYSA-N"
N_CONFORMERS      = 6

USE_QUEUE = True   # set True after restarting the worker


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    config_path = "config/paths.json"
    with open(config_path) as f:
        paths = json.load(f)

    base_dir = paths["base_dir"]

    if USE_QUEUE:
        QueueManager.init_base_dir(base_dir)

    created = []

    for i in range(N_CONFORMERS):
        title = f"pfpuz_conf{i}"

        pipeline_spec = [
            {"stage": "pruning", "args": {
                "select_conf": TARGET_INCHI_KEY,
                "conf_select": [i],
            }},
            {"stage": "optimisation", "args": {"engine": "gxtb_opt_normal"}},
            {"stage": "orcacosmo",    "args": {}},
            {"stage": "solubility",   "args": {}},
        ]

        parameters = {
            "title":     title,
            "config":    paths,
            "resources": {"cpus": 20, "memory_gb": 64},
        }

        req = Request.continue_from(
            base_dir=base_dir,
            parent_request_id=PARENT_REQUEST_ID,
            parent_job_id=PARENT_JOB_ID,
            dataset=title,
            pipeline_spec=pipeline_spec,
            parameters=parameters,
        )

        print(f"[conf {i}] Created: {req.request_id}")
        created.append(req.request_id)

        if USE_QUEUE:
            QueueManager.enqueue(req.request_id)
            print(f"[conf {i}] Enqueued.")

    print(f"\n{N_CONFORMERS} requests created for {TARGET_INCHI_KEY}:")
    for rid in created:
        print(f"  {rid}")

    if not USE_QUEUE:
        print("\nUSE_QUEUE=False — set to True and re-run to enqueue.")


if __name__ == "__main__":
    main()
