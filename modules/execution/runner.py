#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from modules.build.request_manager import Request
from modules.build.job_manager import Job
import os
import json

class PipelineRunner:
    """
    Minimal, pure executor for a Request.

    Responsibilities:
      - load the Request
      - run each Job sequentially
      - write logs to the request directory
      - update pipeline_state.json via Request/Job APIs

    No:
      - nohup
      - scheduler lock
      - pause/resume/stop
      - control.json
      - resource manager
      - CLI entrypoint
    """

    @classmethod
    def run_request(cls, request_id: str, base_dir: str):
        """
        Entry point used exclusively by QueueWorker.
        """
        req = Request.load(base_dir, request_id)
        req.log_header("Pipeline started")
        cls._execute_pipeline(req)

    @classmethod
    def _execute_pipeline(cls, req: Request):
        """
        Sequential job execution.
        """
        if not req.jobs:
            raise RuntimeError(f"Request {req.request_id} has no jobs registered.")

        # Determine starting point
        start_index = cls._find_resume_index(req)

        current_index = start_index
        current_job_id = req.jobs[current_index]
        current_job = Job.load(req, current_job_id)


        while True:
            req.log(f"Running job: {current_job.job_id}")
            current_job.run()
            req.log(f"Completed job: {current_job.job_id}")

            # Create next job based on index
            next_job = req.create_next_job_by_index(current_index)
            if next_job is None:
                req.log_header("Pipeline complete")
                return

            current_index += 1
            current_job = next_job

    @classmethod
    def _find_resume_index(cls, req: Request):
        """
        Determine which job index to start from:
        - If no jobs completed → start at 0
        - If some jobs completed → resume from next unfinished job
        """
        for i, job_id in enumerate(req.jobs):
            job = Job.load(req, job_id)
            state_path = job.job_state_path

            if not os.path.exists(state_path):
                return i  # job never started

            with open(state_path) as f:
                state = json.load(f)

            if state.get("status") != "completed":
                return i  # resume here

        # All jobs completed
        return len(req.jobs) - 1
