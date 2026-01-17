#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import re
import socket
import getpass
from datetime import datetime

from modules.provenance.log_helper import LogHelper
from modules.provenance.job_manager import Job


class RequestError(Exception):
    pass


class Request:
    """
    A Request represents a top-level pipeline invocation.
    It stores:
      - immutable user intent (dataset, pipeline_sequence, stage_args, parameters)
      - provenance (user, hostname, timestamp, pipeline version)
      - job lineage (list of job_ids)
      - parent_request (for branching)
      - continuation info (continued_from_request, continued_from_job)
    """

    # ------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------

    @classmethod
    def create_new(cls, base_dir, dataset, pipeline_spec, parameters=None):
        parameters = parameters or {}

        # Preserve order AND duplicates
        pipeline_sequence = [step["stage"] for step in pipeline_spec]
        stage_args = [step.get("args", {}) for step in pipeline_spec]

        # 1. Generate ID
        request_id = cls._generate_request_id(parameters)

        # 2. Construct the Request object
        request = cls(
            base_dir=base_dir,
            request_id=request_id,
            dataset=dataset,
            pipeline_sequence=pipeline_sequence,
            stage_args=stage_args,
            parameters=parameters,
            parent_request=None,
            continued_from_request=None,
            continued_from_job=None,
            load_existing=False,
        )

        # 3. Logging
        request.log_header(f"Request {request.request_id} created")
        request.log(f"User: {getpass.getuser()}")
        request.log(f"Host: {socket.gethostname()}")
        request.log(f"Pipeline sequence: {request.pipeline_sequence}")

        # 4. Create first job (index 0) — unified logic
        first_stage = pipeline_sequence[0]
        first_args = {**parameters, **stage_args[0]}

        job = Job.create_new(
            request=request,
            stage=first_stage,
            parameters=first_args,
        )

        # 5. Register job ONCE
        request.register_job(job.job_id)
        request.log(f"Created first job: {job.job_id}")

        return request

    @classmethod
    def continue_from(cls, base_dir, parent_request_id, parent_job_id,
                      dataset, pipeline_spec, parameters=None):

        parameters = parameters or {}

        pipeline_sequence = [step["stage"] for step in pipeline_spec]
        stage_args = [step.get("args", {}) for step in pipeline_spec]

        # Generate new request ID
        request_id = cls._generate_request_id(parameters)

        # Construct the Request object
        request = cls(
            base_dir=base_dir,
            request_id=request_id,
            dataset=dataset,
            pipeline_sequence=pipeline_sequence,
            stage_args=stage_args,
            parameters=parameters,
            parent_request=parent_request_id,
            continued_from_request=parent_request_id,
            continued_from_job=parent_job_id,
            load_existing=False,
        )

        # Log creation header
        request.log_header(f"Request {request.request_id} created (continuation)")
        request.log(f"User: {getpass.getuser()}")
        request.log(f"Host: {socket.gethostname()}")
        request.log(f"Pipeline sequence: {request.pipeline_sequence}")

        # Log continuation info
        request.log_section(f"Continuing from {parent_request_id}/{parent_job_id}")

        # First job uses previous job's energies.json
        first_stage = pipeline_sequence[0]
        summary_file = request.get_job_output(parent_job_id, "energies.json")

        first_args = {
            **parameters,
            **stage_args[0],
            "summary_file": summary_file,
        }

        job = Job.create_new(
            request=request,
            stage=first_stage,
            parameters=first_args,
        )

        request.register_job(job.job_id)
        request.log(f"Created first job: {job.job_id}")

        return request

    @classmethod
    def load(cls, base_dir, request_id):
        """Load an existing Request for resuming."""
        return cls(
            base_dir=base_dir,
            request_id=request_id,
            dataset=None,
            pipeline_sequence=None,
            stage_args=None,
            parameters=None,
            parent_request=None,
            continued_from_request=None,
            continued_from_job=None,
            load_existing=True,
        )

    # ------------------------------------------------------------
    # Internal constructor
    # ------------------------------------------------------------
    def __init__(self, base_dir, request_id, dataset, pipeline_sequence,
                 stage_args, parameters, parent_request, continued_from_request,
                 continued_from_job, load_existing):

        self.base_dir = base_dir
        self.requests_root = os.path.join(base_dir, "requests")
        os.makedirs(self.requests_root, exist_ok=True)

        self.request_id = request_id
        self.request_dir = os.path.join(self.requests_root, request_id)
        os.makedirs(self.request_dir, exist_ok=True)

        self.request_log_path = os.path.join(self.request_dir, "request.log")

        self.jobs_dir = os.path.join(self.request_dir, "jobs")
        os.makedirs(self.jobs_dir, exist_ok=True)

        self.request_json_path = os.path.join(self.request_dir, "request.json")

        if load_existing:
            self._load_request_json()
        else:
            self.dataset = dataset
            self.pipeline_sequence = pipeline_sequence
            self.stage_args = stage_args or []
            self.parameters = parameters or {}

            # Resource + detached handling
            from modules.execution.resource_manager import ResourceManager
            self.resources = ResourceManager.get_resources(self.parameters)
            self.detached = self.parameters.get("detached", False)

            self.parent_request = parent_request
            self.continued_from_request = continued_from_request
            self.continued_from_job = continued_from_job

            self.jobs = []
            self._write_request_json()

    # ------------------------------------------------------------
    # ID generation
    # ------------------------------------------------------------
    @classmethod
    def _generate_request_id(cls, parameters):
        timestamp = datetime.now().strftime("%d%b%y-%H%M-%S%f")[:-3]
        title = parameters.get("title") if parameters else None

        if title:
            safe_title = re.sub(r"[^A-Za-z0-9_]+", "_", str(title).strip())
            return f"R-{timestamp}-{safe_title}"

        return f"R-{timestamp}"

    # ------------------------------------------------------------
    # JSON handling
    # ------------------------------------------------------------
    def _write_request_json(self):
        data = {
            "request_id": self.request_id,
            "created_at": datetime.now().isoformat(),
            "user": getpass.getuser(),
            "hostname": socket.gethostname(),
            "pipeline_version": "1.0",

            "dataset": self.dataset,
            "pipeline_sequence": self.pipeline_sequence,
            "stage_args": self.stage_args,
            "parameters": self.parameters,

            "parent_request": self.parent_request,
            "continued_from_request": self.continued_from_request,
            "continued_from_job": self.continued_from_job,

            "jobs": getattr(self, "jobs", []),
        }

        with open(self.request_json_path, "w") as f:
            json.dump(data, f, indent=2)

    def _load_request_json(self):
        with open(self.request_json_path) as f:
            data = json.load(f)

        LogHelper.write(self.request_log_path, "Loaded existing request")

        self.dataset = data["dataset"]
        self.pipeline_sequence = data["pipeline_sequence"]

        stage_args = data.get("stage_args", [])
        if isinstance(stage_args, dict):
            # legacy dict format
            self.stage_args = [stage_args.get(stage, {}) for stage in self.pipeline_sequence]
        else:
            self.stage_args = stage_args

        self.parameters = data.get("parameters", {})

        self.parent_request = data.get("parent_request")
        self.continued_from_request = data.get("continued_from_request")
        self.continued_from_job = data.get("continued_from_job")

        self.jobs = data.get("jobs", [])

        # Recompute resources + detached
        from modules.execution.resource_manager import ResourceManager
        self.resources = ResourceManager.get_resources(self.parameters)
        self.detached = self.parameters.get("detached", False)

    # ------------------------------------------------------------
    # Job registration
    # ------------------------------------------------------------
    def register_job(self, job_id):
        self.jobs.append(job_id)
        self._write_request_json()
        LogHelper.write(self.request_log_path, f"Registered job: {job_id}")

    # ------------------------------------------------------------
    # Legacy name-based sequencing (kept for compatibility)
    # ------------------------------------------------------------
    def get_next_stage(self, current_stage):
        if self.pipeline_sequence is None:
            raise RequestError("pipeline_sequence is not defined.")

        try:
            idx = self.pipeline_sequence.index(current_stage)
        except ValueError:
            raise RequestError(f"Stage {current_stage} not in pipeline_sequence.")

        return self.pipeline_sequence[idx + 1] if idx + 1 < len(self.pipeline_sequence) else None

    # ------------------------------------------------------------
    # Preferred index-based sequencing
    # ------------------------------------------------------------
    def create_next_job_by_index(self, current_index, extra_parameters=None):
        next_index = current_index + 1
        if next_index >= len(self.pipeline_sequence):
            return None

        next_stage = self.pipeline_sequence[next_index]
        stage_args = self.stage_args[next_index] if next_index < len(self.stage_args) else {}
        extra_parameters = extra_parameters or {}

        # Get previous job's summary file using Job.STAGE_OUTPUTS
        summary_file = None
        if 0 <= current_index < len(self.jobs):
            prev_job_id = self.jobs[current_index]
            prev_stage = self.pipeline_sequence[current_index]
            
            # Use Job's stage output mapping
            if prev_stage in Job.STAGE_OUTPUTS:
                summary_filename = Job.STAGE_OUTPUTS[prev_stage]
                prev_outputs = os.path.join(self.jobs_dir, prev_job_id, "outputs")
                summary_file = os.path.join(prev_outputs, summary_filename)

        params = {
            **self.parameters,
            **stage_args,
            **extra_parameters,
        }

        if summary_file is not None:
            params["summary_file"] = summary_file

        job = Job.create_new(
            request=self,
            stage=next_stage,
            parameters=params,
        )

        self.register_job(job.job_id)
        self.log(f"Created next job: {job.job_id} (stage: {next_stage}, index: {next_index})")

        return job


    # ------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------
    def summary_path(self):
        return os.path.join(self.request_dir, f"_pipeline_summary_{self.request_id}.log")

    def log(self, message, indent=0, echo=False):
        LogHelper.write(self.request_log_path, message, indent, echo)

    def log_header(self, title, echo=False):
        LogHelper.header(self.request_log_path, title, echo)

    def log_section(self, title, echo=False):
        LogHelper.section(self.request_log_path, title, echo)

    def get_job_output(self, job_id, filename):
        return os.path.join(self.jobs_dir, job_id, "outputs", filename)
