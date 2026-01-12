import os
import json
import time

class PipelineControl:
    """
    Control-plane operations for the pipeline:
      - pause(request_id)
      - resume(request_id)
      - stop(request_id)
      - is_paused(request_id)
      - is_stopped(request_id)
    """

    def __init__(self, base_dir):
        self.base_dir = base_dir

    # ------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------
    def _control_path(self, request_id):
        return os.path.join(self.base_dir, "requests", request_id, "control.json")

    def _load_control(self, request_id):
        path = self._control_path(request_id)
        if not os.path.exists(path):
            return {"pause": False, "stop": False}
        with open(path) as f:
            return json.load(f)

    def _write_control(self, request_id, data):
        path = self._control_path(request_id)
        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    # ------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------
    def pause(self, request_id):
        ctrl = self._load_control(request_id)
        ctrl["pause"] = True
        self._write_control(request_id, ctrl)

    def resume(self, request_id):
        ctrl = self._load_control(request_id)
        ctrl["pause"] = False
        self._write_control(request_id, ctrl)

    def stop(self, request_id):
        ctrl = self._load_control(request_id)
        ctrl["stop"] = True
        self._write_control(request_id, ctrl)

    def is_paused(self, request_id):
        return self._load_control(request_id).get("pause", False)

    def is_stopped(self, request_id):
        return self._load_control(request_id).get("stop", False)
