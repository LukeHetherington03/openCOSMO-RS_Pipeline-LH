#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import re
from datetime import datetime
import shutil
from modules.build.log_helper import LogHelper


class JobError(Exception):
    pass


class Job:
    """
    A Job represents a single stage execution within a Request.

    New architecture:
      - Each stage produces exactly one canonical output file.
      - The next stage receives exactly one input: stage_input.
      - No auto-detection, no summary logic, no guessing.
      - Job handles input→output chaining cleanly and deterministically.
    """

    # Canonical output filenames for each stage
    STAGE_OUTPUTS = {
        "cleaning": "cleaned.csv",
        "generation": "energies.json",
        "pruning": "energies.json",
        "optimisation": "energies.json",
        "orcacosmo": "orcacosmo_summary.json",
        "solubility": "solubility_results.json",
    }

    # ------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------
    @classmethod
    def create_new(cls, request, stage, parameters=None, previous_job=None):
        """
        Create a new job, optionally chaining from a previous job.

        Automatically injects:
            stage_input = previous_job.output_file
        """
        job_id = cls._generate_job_id(stage)
        job_dir = os.path.join(request.jobs_dir, job_id)
        os.makedirs(job_dir, exist_ok=True)

        final_parameters = parameters or {}

        # Inject stage_input from previous job
        if previous_job is not None:
            output_file = cls._get_stage_output(previous_job)
            if output_file:
                # Copy previous stage output into this job's inputs
                dst = os.path.join(job_dir, "inputs", os.path.basename(output_file))
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy(output_file, dst)

                # Use the copied file as stage_input
                final_parameters["stage_input"] = dst
                final_parameters["previous_job_id"] = previous_job.job_id
                final_parameters["previous_stage"] = previous_job.stage

        job = cls(
            request=request,
            job_id=job_id,
            stage=stage,
            job_dir=job_dir,
            parameters=final_parameters,
            load_existing=False,
        )

        job._initialise_directories()
        job._write_parameters_json()

        # Initial job state
        job._write_job_state(
            status="initialised",
            created_at=datetime.now().isoformat(),
            completed_at=None,
            error_message=None,
            output_file=None,
        )

        return job

    @classmethod
    def load(cls, request, job_id):
        job_dir = os.path.join(request.jobs_dir, job_id)
        if not os.path.exists(job_dir):
            raise JobError(f"Job directory not found: {job_dir}")

        return cls(
            request=request,
            job_id=job_id,
            stage=None,
            job_dir=job_dir,
            parameters=None,
            load_existing=True,
        )

    # ------------------------------------------------------------
    # Output file lookup
    # ------------------------------------------------------------
    @classmethod
    def _get_stage_output(cls, job):
        """
        Return the canonical output file for a completed job.
        """
        if job.stage not in cls.STAGE_OUTPUTS:
            return None

        output_filename = cls.STAGE_OUTPUTS[job.stage]
        output_path = os.path.join(job.outputs_dir, output_filename)

        return output_path if os.path.exists(output_path) else None

    # ------------------------------------------------------------
    # Internal constructor
    # ------------------------------------------------------------
    def __init__(self, request, job_id, stage, job_dir, parameters, load_existing):
        self.request = request
        self.request_id = request.request_id
        self.job_id = job_id
        self.job_dir = job_dir
        self.stage = stage
        self.parameters = parameters or {}
        self.config = request.parameters.get("config", {})

        self.job_state_path = os.path.join(job_dir, "job_state.json")
        self.stage_log_path = os.path.join(job_dir, "stage.log")

        self.inputs_dir = os.path.join(job_dir, "inputs")
        self.outputs_dir = os.path.join(job_dir, "outputs")

        # Item tracking
        self.items = []
        self.pending_items = []
        self.completed_items = []
        self.failed_items = []

        if load_existing:
            self._load_job_state()

    # ------------------------------------------------------------
    # Parameter handling
    # ------------------------------------------------------------
    def _write_parameters_json(self):
        params_path = os.path.join(self.inputs_dir, "parameters.json")
        with open(params_path, "w") as f:
            json.dump(self.parameters, f, indent=2)

    # ------------------------------------------------------------
    # ID generation
    # ------------------------------------------------------------
    @classmethod
    def _generate_job_id(cls, stage):
        timestamp = datetime.now().strftime("%d%b%y-%H%M-%S%f")[:-3]
        safe_stage = re.sub(r"[^A-Za-z0-9_]+", "_", stage.strip().lower())
        return f"J-{timestamp}-{safe_stage}"

    # ------------------------------------------------------------
    # Directory initialisation
    # ------------------------------------------------------------
    def _initialise_directories(self):
        os.makedirs(self.inputs_dir, exist_ok=True)
        os.makedirs(self.outputs_dir, exist_ok=True)

    # ------------------------------------------------------------
    # Job state handling
    # ------------------------------------------------------------
    def _write_job_state(
        self,
        status=None,
        created_at=None,
        completed_at=None,
        error_message=None,
        output_file=None,
    ):
        state = {
            "job_id": self.job_id,
            "stage": self.stage,
            "parameters": self.parameters,
            "status": status,
            "created_at": created_at,
            "completed_at": completed_at,
            "error_message": error_message,
            "output_file": output_file,
            "items": self.items,
            "pending_items": self.pending_items,
            "completed_items": self.completed_items,
            "failed_items": self.failed_items,
        }

        with open(self.job_state_path, "w") as f:
            json.dump(state, f, indent=2)

    def _load_job_state(self):
        if not os.path.exists(self.job_state_path):
            raise JobError(f"Missing job_state.json for {self.job_id}")

        with open(self.job_state_path) as f:
            state = json.load(f)

        self.stage = state["stage"]
        self.parameters = state["parameters"]

        # Restore item tracking
        self.items = state.get("items", [])
        self.pending_items = state.get("pending_items", [])
        self.completed_items = state.get("completed_items", [])
        self.failed_items = state.get("failed_items", [])

    def mark_running(self):
        self._write_job_state(
            status="running",
            created_at=datetime.now().isoformat(),
            completed_at=None,
            error_message=None,
            output_file=None,
        )

    def mark_failed(self, error_message):
        self._write_job_state(
            status="failed",
            created_at=None,
            completed_at=datetime.now().isoformat(),
            error_message=str(error_message),
            output_file=None,
        )

    def mark_complete(self):
        output_file = self.get_output_file()
        self._write_job_state(
            status="completed",
            created_at=None,
            completed_at=datetime.now().isoformat(),
            error_message=None,
            output_file=output_file,
        )

    # ------------------------------------------------------------
    # Item tracking API
    # ------------------------------------------------------------
    def set_items(self, items):
        self.items = list(items)
        self.pending_items = list(items)
        self.completed_items = []
        self.failed_items = []

        self._write_job_state(
            status="running",
            created_at=datetime.now().isoformat(),
            completed_at=None,
            error_message=None,
            output_file=None,
        )

    def update_progress(self, item, success=True):
        if item in self.pending_items:
            self.pending_items.remove(item)

        if success:
            self.completed_items.append(item)
        else:
            self.failed_items.append(item)

        self._write_job_state(
            status="running",
            created_at=None,
            completed_at=None,
            error_message=None,
            output_file=None,
        )

    # ------------------------------------------------------------
    # Logging
    # ------------------------------------------------------------
    def log(self, message, indent=0, echo=False):
        LogHelper.write(self.stage_log_path, message, indent, echo)

    def log_header(self, title, echo=False):
        LogHelper.header(self.stage_log_path, title, echo)

    def log_section(self, title, echo=False):
        LogHelper.section(self.stage_log_path, title, echo)

    # ------------------------------------------------------------
    # Stage execution
    # ------------------------------------------------------------
    def run(self):
        module_name = f"modules.stages.{self.stage}_stage"
        class_name = f"{self.stage.capitalize()}Stage"

        try:
            module = __import__(module_name, fromlist=[class_name])
            StageClass = getattr(module, class_name)
        except Exception as e:
            self.mark_failed(f"Failed to load stage module: {e}")
            raise

        stage_instance = StageClass(self)

        self.mark_running()

        try:
            stage_instance.run()
            self.mark_complete()
        except Exception as e:
            self.mark_failed(e)
            raise

    # ------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------
    def input_path(self, *parts):
        return os.path.join(self.inputs_dir, *parts)

    def output_path(self, *parts):
        return os.path.join(self.outputs_dir, *parts)

    def get_output_file(self):
        """Return the canonical output file for this job."""
        if self.stage in self.STAGE_OUTPUTS:
            return os.path.join(self.outputs_dir, self.STAGE_OUTPUTS[self.stage])
        return None

    def get_request_id(self):
        return self.request.request_id

    def get_job_id(self):
        return self.job_id
