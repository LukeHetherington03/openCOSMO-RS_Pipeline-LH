# modules/stages/base_stage.py

from datetime import datetime

class BaseStage:
    """
    Base class for all pipeline stages.
    Provides:
      - job lifecycle integration (running → complete/failed)
      - logging helpers
      - item tracking helpers
      - path helpers
      - fail() helper
    """

    def __init__(self, job):
        self.job = job
        self.inputs_dir = job.inputs_dir
        self.outputs_dir = job.outputs_dir
        self.parameters = job.parameters

    # ------------------------------------------------------------
    # Unified lifecycle wrapper
    # ------------------------------------------------------------
    def run(self):
        """
        DO NOT override this in child classes.
        Child stages must implement execute().
        """
        # Mark job as running
        self.job.mark_running()

        try:
            # Stage-specific logic
            self.execute()

            # Mark success
            self.job.mark_complete()

        except Exception as e:
            # Mark failure
            self.job.mark_failed(str(e))
            raise

    # ------------------------------------------------------------
    # Child classes must override this
    # ------------------------------------------------------------
    def execute(self):
        raise NotImplementedError("Stage must implement execute()")

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
    # Error helper
    # ------------------------------------------------------------
    def fail(self, message):
        self.log(f"[ERROR] {message}")
        raise RuntimeError(message)
