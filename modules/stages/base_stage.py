from datetime import datetime
import os
import time
import json

class BaseStage:
    """
    Base class for all pipeline stages.
    Provides:
      - unified lifecycle (running → complete/failed)
      - logging helpers
      - item tracking
      - path helpers
      - strict mode helpers
      - stage_input / stage_output helpers
    """

    def __init__(self, job):
        self.job = job
        self.request = job.request
        self.request_id = job.request_id
        self.job_id = job.job_id

        self.inputs_dir = job.inputs_dir
        self.outputs_dir = job.outputs_dir
        self.parameters = job.parameters
        self.config = job.config
        self._stage_output = None

    # ------------------------------------------------------------
    # Unified lifecycle wrapper
    # ------------------------------------------------------------
    def run(self):
        # Debug context
        self.log_header(f"Starting {self.__class__.__name__}")
        self.log("[DEBUG] Stage context:")
        self.log(f"  request_id: {self.request_id}")
        self.log(f"  job_id: {self.job_id}")
        self.log(f"  inputs_dir: {self.inputs_dir}")
        self.log(f"  outputs_dir: {self.outputs_dir}")
        self.log(f"  parameters: {json.dumps(self.parameters, indent=2)}")

        self.job.mark_running()

        start = time.perf_counter()


        try:
            self.execute()
            elapsed = time.perf_counter() - start
            self.log(f"[INFO] Stage completed in {elapsed:.2f} seconds")
            self.job.mark_complete()

            # SUCCESS: update pipeline state
            self.request.update_pipeline_state(
                last_completed_stage=self.job.stage,
                last_completed_job_id=self.job.job_id,
            )

        except Exception as e:
            self.log_header("FATAL ERROR — STAGE HALTED")
            self.log(f"[ERROR] {e}")

            # Mark job failed
            self.job.mark_failed(str(e))

            # FAILURE: update pipeline state
            self.request.update_pipeline_state(
                state="failed",
                halt_reason=str(e),
                current_stage=self.job.stage,
                current_job_id=self.job.job_id,
            )

            raise




    # ------------------------------------------------------------
    # Child classes must override this
    # ------------------------------------------------------------
    def execute(self):
        raise NotImplementedError("Stage must implement execute()")

    # ------------------------------------------------------------
    # Stage input / output helpers
    # ------------------------------------------------------------
    def get_stage_input(self):
        stage_input = self.parameters.get("stage_input")
        if not stage_input:
            self.fail("Stage requires 'stage_input' but none was provided.")
        return stage_input

    def set_stage_output(self, filename):
        """Declare the canonical output file for this stage."""
        self._stage_output = self.output_path(filename)
        return self._stage_output

    def get_stage_output(self):
        return self._stage_output

    # ------------------------------------------------------------
    # Strict mode helper
    # ------------------------------------------------------------
    def strict(self, section):
        return bool(self.config.get(section, {}).get("strict", False))

    # ------------------------------------------------------------
    # Logging helpers
    # ------------------------------------------------------------
    def log(self, msg, indent=0):
        self.job.log(msg, indent)

    def log_header(self, title):
        self.job.log_header(title)

    def log_section(self, title):
        self.job.log_section(title)

    # ------------------------------------------------------------
    # Item tracking
    # ------------------------------------------------------------
    def set_items(self, items):
        self.job.set_items(items)

    def update_progress(self, item, success=True):
        self.job.update_progress(item, success)

    # ------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------
    def input_path(self, *parts):
        return self.job.input_path(*parts)

    def output_path(self, *parts):
        return self.job.output_path(*parts)

    # ------------------------------------------------------------
    # File existence helper
    # ------------------------------------------------------------
    def require_file(self, path, description="required file"):
        if not os.path.exists(path):
            self.fail(f"{description} not found: {path}")
        return path

    # ------------------------------------------------------------
    # Error helper
    # ------------------------------------------------------------
    def fail(self, message):
        self.log(f"[ERROR] {message}")
        raise RuntimeError(message)
