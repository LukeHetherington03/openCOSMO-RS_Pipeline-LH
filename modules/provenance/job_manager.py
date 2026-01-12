#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import shutil
from datetime import datetime
import re
from modules.provenance.log_helper import LogHelper


class JobError(Exception):
    pass


class Job:
    """
    A Job represents a single stage execution within a Request.

    Directory structure:

        job/
            job_state.json
            stage.log
            inputs/
                parameters.json
                ... (stage-specific inputs)
            outputs/
                energies.json
                summary.csv
                xyz/
                ... (stage-specific outputs)

    The Job object:
      - creates the job directory
      - manages inputs/ and outputs/ folders
      - writes job_state.json
      - writes stage.log
      - runs the stage module
    """

    # ------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Create a new job
    # ------------------------------------------------------------
    @classmethod
    def create_new(cls, request, stage, parameters=None):
        """
        Create a new Job belonging to a Request.
        """

        # Generate job ID with milliseconds
        job_id = cls._generate_job_id(stage)

        # Job directory (FIXED: must use request.jobs_dir)
        job_dir = os.path.join(request.jobs_dir, job_id)
        os.makedirs(job_dir, exist_ok=True)

        # Instantiate job
        job = cls(
            request=request,
            job_id=job_id,
            stage=stage,
            job_dir=job_dir,
            parameters=parameters or {},
            load_existing=False,
        )

        # Initialise directories
        job._initialise_directories()

        # Write parameters.json
        job._write_parameters_json()

        # Write initial job_state.json
        job._write_job_state(initial=True)

        return job


    @classmethod
    def load(cls, request, job_id):
        """
        Load an existing job from disk.
        """
        job_dir = os.path.join(request.jobs_dir, job_id)
        if not os.path.exists(job_dir):
            raise JobError(f"Job directory not found: {job_dir}")

        return cls(
            request=request,
            job_id=job_id,
            stage=None,
            job_dir=job_dir,
            parameters=None,
            load_existing=True
        )

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

        if load_existing:
            self._load_job_state()




    # ------------------------------------------------------------
    # ID generation
    # ------------------------------------------------------------


    @classmethod
    def _generate_job_id(cls, stage):
        """
        Generate a job ID of the form:
            J-10Jan26-1447-02014-generation
        """

        # Timestamp with milliseconds
        timestamp = datetime.now().strftime("%d%b%y-%H%M-%S%f")[:-3]

        # Clean stage name
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
    def _write_job_state(self, initial=False):
        state = {
            "job_id": self.job_id,
            "stage": self.stage,
            "created_at": datetime.now().isoformat(),
            "parameters": self.parameters,
            "status": "initialised" if initial else "completed",
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

    def mark_complete(self):
        self._write_job_state(initial=False)

    # ------------------------------------------------------------
    # Parameter handling
    # ------------------------------------------------------------
    def _write_parameters_json(self):
        params_path = os.path.join(self.inputs_dir, "parameters.json")
        with open(params_path, "w") as f:
            json.dump(self.parameters, f, indent=2)

    # ------------------------------------------------------------
    # Logging
    # ------------------------------------------------------------
    def write_pipeline_log(self, message):
        with open(self.stage_log_path, "a") as f:
            f.write(f"[{datetime.now().isoformat()}] {message}\n")

    # ------------------------------------------------------------
    # Stage execution
    # ------------------------------------------------------------
    def run(self):
        """
        Dynamically import and run the stage module.
        """
        module_name = f"modules.stages.{self.stage}_stage"
        class_name = f"{self.stage.capitalize()}Stage"

        try:
            module = __import__(module_name, fromlist=[class_name])
            StageClass = getattr(module, class_name)
        except Exception as e:
            raise JobError(f"Failed to load stage module {module_name}: {e}")

        stage_instance = StageClass()
        stage_instance.run(self)

    # ------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------
    def input_path(self, *parts):
        return os.path.join(self.inputs_dir, *parts)

    def output_path(self, *parts):
        return os.path.join(self.outputs_dir, *parts)

    def log(self, message, indent=0, echo=False):
        LogHelper.write(self.stage_log_path, message, indent, echo)

    def log_header(self, title, echo=False):
        LogHelper.header(self.stage_log_path, title, echo)

    def log_section(self, title, echo=False):
        LogHelper.section(self.stage_log_path, title, echo)
