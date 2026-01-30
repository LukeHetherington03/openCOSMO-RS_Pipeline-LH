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
      - automatically detaches from SSH using nohup
      - ensures only one pipeline runs at a time (Scheduler)
      - applies resource allocation (ResourceManager)
      - runs each Job in sequence (index-based)
      - obeys pause/resume/stop signals (PipelineControl)
      - logs progress
    """

    # ------------------------------------------------------------
    # Auto-detach using nohup
    # ------------------------------------------------------------
    @classmethod
    def _ensure_nohup(cls, request_id, base_dir):
        print("[DEBUG] Entered _ensure_nohup()")

        # If stdout is a TTY, we are running interactively (SSH)
        if os.isatty(sys.stdout.fileno()):
            print("[DEBUG] Detected interactive TTY — performing nohup relaunch")

            cmd = [
                "nohup",
                sys.executable,
                "-m",
                "modules.execution.pipeline_runner",
                request_id,
                base_dir
            ]

            print(f"[DEBUG] Relaunch command: {' '.join(cmd)}")

            subprocess.Popen(
                cmd,
                stdout=open(os.devnull, "w"),
                stderr=open(os.devnull, "w")
            )

            print(f"[INFO] Pipeline re-launched under nohup for request {request_id}. Safe to disconnect.")
            sys.exit(0)

        print("[DEBUG] No TTY detected — already running under nohup")

    # ------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------
    @classmethod
    def run_request(cls, request_id, base_dir):
        print("[DEBUG] Entered run_request()")
        print(f"[DEBUG] request_id={request_id}, base_dir={base_dir}")

        # Ensure nohup detachment
        cls._ensure_nohup(request_id, base_dir)

        print("[DEBUG] Loading Request object...")
        req = Request.load(base_dir, request_id)
        print("[DEBUG] Request loaded successfully")

        cls._run_attached(req)

    # ------------------------------------------------------------
    # Attached execution (nohup-safe)
    # ------------------------------------------------------------
    @classmethod
    def _run_attached(cls, req):
        print("[DEBUG] Entered _run_attached()")
        base_dir = req.base_dir

        print(f"Starting pipeline for Request {req.request_id}")
        req.log_header("Pipeline started")

        print("[DEBUG] Attempting to acquire scheduler lock...")
        Scheduler.acquire_lock(base_dir)
        print("[DEBUG] Scheduler lock acquired")

        try:
            cls._execute_pipeline(req)
        finally:
            print("[DEBUG] Releasing scheduler lock")
            Scheduler.release_lock(base_dir)

    # ------------------------------------------------------------
    # Core pipeline execution (INDEX-BASED)
    # ------------------------------------------------------------
    @classmethod
    def _execute_pipeline(cls, req):
        print("[DEBUG] Entered _execute_pipeline()")

        control = PipelineControl(req.base_dir)

        if not req.jobs:
            raise RuntimeError("Request has no jobs registered.")

        print(f"[DEBUG] Jobs registered: {req.jobs}")

        # Start at the first job (index 0)
        current_index = 0
        current_job_id = req.jobs[current_index]
        print(f"[DEBUG] Loading first job: {current_job_id}")

        current_job = Job.load(req, current_job_id)

        # Apply resource manager
        print("[DEBUG] Applying resource manager...")
        req.resources = ResourceManager.get_resources(req.parameters)
        req.log(f"Using resources: {req.resources}")
        print(f"[DEBUG] Resources applied: {req.resources}")

        # Sequential job execution
        while True:
            print(f"[DEBUG] Loop start — current_index={current_index}, job_id={current_job.job_id}")

            # STOP CHECK
            if control.is_stopped(req.request_id):
                print("[DEBUG] STOP detected")
                req.log("Pipeline STOPPED by user.")
                print("Pipeline stopped by user.")
                return

            # PAUSE CHECK
            if control.is_paused(req.request_id):
                print("[DEBUG] PAUSE detected")
                req.log("Pipeline PAUSED by user.")
                print("Pipeline paused. Waiting for resume...")
                while control.is_paused(req.request_id):
                    time.sleep(1)
                req.log("Pipeline RESUMED by user.")
                print("Pipeline resumed.")

            # RUN JOB
            print(f"[DEBUG] Running job: {current_job.job_id}")
            req.log(f"Running job: {current_job.job_id}")
            current_job.run()
            print(f"[DEBUG] Job completed: {current_job.job_id}")

            # NEXT STAGE (index-based)
            print("[DEBUG] Creating next job...")
            next_job = req.create_next_job_by_index(current_index)

            if next_job is None:
                print("[DEBUG] No next job — pipeline complete")
                req.log_header("Pipeline complete")
                print("Pipeline complete.")
                return

            # Advance
            current_index += 1
            current_job = next_job
            print(f"[DEBUG] Advanced to next job: {current_job.job_id}")

    # ------------------------------------------------------------
    # CLI entry
    # ------------------------------------------------------------
    @classmethod
    def cli(cls):
        print("[DEBUG] Entered CLI mode")
        print(f"[DEBUG] sys.argv={sys.argv}")

        if len(sys.argv) != 3:
            print("Usage: python3 -m modules.execution.pipeline_runner <request_id> <base_dir>")
            sys.exit(1)

        request_id = sys.argv[1]
        base_dir = sys.argv[2]

        cls.run_request(request_id, base_dir)


if __name__ == "__main__":
    PipelineRunner.cli()
