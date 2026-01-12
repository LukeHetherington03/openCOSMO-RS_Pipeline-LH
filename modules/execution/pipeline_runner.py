#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import time

from modules.provenance.request_manager import Request
from modules.provenance.job_manager import Job
from modules.execution.scheduler import Scheduler
from modules.execution.resource_manager import ResourceManager
from modules.execution.pipeline_control import PipelineControl


class PipelineRunner:
    """
    Orchestrates execution of a Request:
      - supports detached/background mode
      - ensures only one pipeline runs at a time (Scheduler)
      - applies resource allocation (ResourceManager)
      - runs each Job in sequence
      - obeys pause/resume/stop signals (PipelineControl)
      - logs progress
    """

    # ------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------
    @classmethod
    def run_request(cls, request_id, base_dir):
        """
        Run or resume a Request end-to-end.
        Handles detached mode and scheduling.
        """
        req = Request.load(base_dir, request_id)

        # Detached mode: re-launch in background
        if req.parameters.get("detached", False):
            cls._run_detached(request_id, base_dir)
            return

        # Foreground mode
        cls._run_attached(req)

    # ------------------------------------------------------------
    # Detached execution
    # ------------------------------------------------------------
    @classmethod
    def _run_detached(cls, request_id, base_dir):
        """
        Relaunch the pipeline in the background using nohup-like behaviour.
        """
        cmd = [
            sys.executable,
            "-m",
            "modules.execution.pipeline_runner",
            request_id,
            base_dir
        ]

        subprocess.Popen(
            cmd,
            stdout=open(os.devnull, "w"),
            stderr=open(os.devnull, "w"),
            preexec_fn=os.setpgrp
        )

        print(f"Pipeline for {request_id} is now running in the background.")

    # ------------------------------------------------------------
    # Attached execution
    # ------------------------------------------------------------
    @classmethod
    def _run_attached(cls, req):
        """
        Run the pipeline in the foreground.
        Ensures only one pipeline runs at a time.
        """
        base_dir = req.base_dir

        print(f"Starting pipeline for Request {req.request_id}")
        req.log_header("Pipeline started")

        # Acquire scheduler lock (blocks until free)
        Scheduler.acquire_lock(base_dir)

        try:
            cls._execute_pipeline(req)
        finally:
            Scheduler.release_lock(base_dir)

    # ------------------------------------------------------------
    # Core pipeline execution
    # ------------------------------------------------------------
    @classmethod
    def _execute_pipeline(cls, req):
        """
        Run all jobs in the Request sequentially.
        Obeys pause/resume/stop signals.
        """
        control = PipelineControl(req.base_dir)

        # Load first job
        if not req.jobs:
            raise RuntimeError("Request has no jobs registered.")

        current_job_id = req.jobs[0]
        current_job = Job.load(req, current_job_id)

        # Apply resource manager
        req.resources = ResourceManager.get_resources(req.parameters)
        req.log(f"Using resources: {req.resources}")

        # Sequential job execution
        while True:

            # --------------------------------------------------------
            # STOP CHECK
            # --------------------------------------------------------
            if control.is_stopped(req.request_id):
                req.log("Pipeline STOPPED by user.")
                print("Pipeline stopped by user.")
                return

            # --------------------------------------------------------
            # PAUSE CHECK
            # --------------------------------------------------------
            if control.is_paused(req.request_id):
                req.log("Pipeline PAUSED by user.")
                print("Pipeline paused. Waiting for resume...")
                while control.is_paused(req.request_id):
                    time.sleep(1)
                req.log("Pipeline RESUMED by user.")
                print("Pipeline resumed.")

            # --------------------------------------------------------
            # RUN JOB
            # --------------------------------------------------------
            req.log(f"Running job: {current_job.job_id}")
            print(f"Running job: {current_job.job_id}")
            current_job.run()

            # --------------------------------------------------------
            # NEXT STAGE
            # --------------------------------------------------------
            next_stage = req.get_next_stage(current_job.stage)
            if next_stage is None:
                req.log_header("Pipeline complete")
                print("Pipeline complete.")
                return

            # Create next job
            next_job = req.create_next_job(current_job)
            current_job = next_job

    # ------------------------------------------------------------
    # CLI entry
    # ------------------------------------------------------------
    @classmethod
    def cli(cls):
        """
        CLI wrapper:
            python3 -m modules.execution.pipeline_runner <request_id> <base_dir>
        """
        if len(sys.argv) != 3:
            print("Usage: python3 -m modules.execution.pipeline_runner <request_id> <base_dir>")
            sys.exit(1)

        request_id = sys.argv[1]
        base_dir = sys.argv[2]

        cls.run_request(request_id, base_dir)


if __name__ == "__main__":
    PipelineRunner.cli()
