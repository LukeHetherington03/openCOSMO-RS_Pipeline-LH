#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
from datetime import datetime
import re
from modules.provenance.log_helper import LogHelper


class JobError(Exception):
    pass


class Job:
    """
    A Job represents a single stage execution within a Request.

    Improvements:
      - job_state.json is written immediately
      - status transitions: initialised → running → complete/failed
      - crash-safe: failures are recorded with error_message
      - item tracking preserved across crashes
      - stage output chaining for pipeline flow
    """

    # Stage output mapping - defines which summary file each stage produces
    STAGE_OUTPUTS = {
        "cleaning": "clean_dataset.csv",
        "generation": "energies.json",
        "pruning": "energies.json",
        "optimisation": "energies.json",
        "orcacosmo": "orcacosmo_summary.json",
        "solubility": "solubility_summary.json",
    }

    # ------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------
    @classmethod
    def create_new(cls, request, stage, parameters=None, previous_job=None):
        """
        Create a new job, optionally chaining from a previous job.
        
        Args:
            request: Parent Request object
            stage: Stage name (e.g., "orcacosmo", "solubility")
            parameters: Dict of stage parameters
            previous_job: Previous Job object to chain from (optional)
        """
        job_id = cls._generate_job_id(stage)
        job_dir = os.path.join(request.jobs_dir, job_id)
        os.makedirs(job_dir, exist_ok=True)

        # Merge parameters with previous stage outputs
        final_parameters = parameters or {}
        
        if previous_job is not None:
            # Auto-detect and pass the previous stage's summary file
            summary_file = cls._get_stage_summary(previous_job)
            if summary_file:
                final_parameters["summary_file"] = summary_file
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

        # Write initial job state immediately
        job._write_job_state(
            status="initialised",
            created_at=datetime.now().isoformat(),
            completed_at=None,
            error_message=None,
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

    @classmethod
    def _get_stage_summary(cls, job):
        """
        Get the summary output file from a completed job.
        
        Args:
            job: Job object
            
        Returns:
            Path to summary file if it exists, None otherwise
        """
        if job.stage not in cls.STAGE_OUTPUTS:
            return None
            
        summary_filename = cls.STAGE_OUTPUTS[job.stage]
        summary_path = os.path.join(job.outputs_dir, summary_filename)
        
        if os.path.exists(summary_path):
            return summary_path
        else:
            return None

    # ------------------------------------------------------------
    # Internal constructor
    # ------------------------------------------------------------
    def __init__(self, request, job_id, stage, job_dir, parameters, load_existing):
        self.request = request
        self.job_id = job_id
        self.job_dir = job_dir
        self.stage = stage
        self.parameters = parameters or {}

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
    ):
        state = {
            "job_id": self.job_id,
            "stage": self.stage,
            "parameters": self.parameters,
            "status": status,
            "created_at": created_at,
            "completed_at": completed_at,
            "error_message": error_message,
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
        )

    def mark_failed(self, error_message):
        self._write_job_state(
            status="failed",
            created_at=None,
            completed_at=datetime.now().isoformat(),
            error_message=str(error_message),
        )

    def mark_complete(self):
        self._write_job_state(
            status="completed",
            created_at=None,
            completed_at=datetime.now().isoformat(),
            error_message=None,
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
    # Stage execution (crash-safe)
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

        # Mark job as running before execution
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

    def get_summary_path(self):
        """Get the expected summary output path for this job's stage."""
        if self.stage in self.STAGE_OUTPUTS:
            return os.path.join(self.outputs_dir, self.STAGE_OUTPUTS[self.stage])
        return None