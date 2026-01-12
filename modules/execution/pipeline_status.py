import os
import json
import time
from datetime import datetime

class PipelineStatus:
    """
    Status-plane operations:
      - list_active()
      - get_status(request_id)
      - summary()
    """

    def __init__(self, base_dir):
        self.base_dir = base_dir

    # ------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------
    def _requests_root(self):
        return os.path.join(self.base_dir, "requests")

    def _job_state_path(self, req_id, job_id):
        return os.path.join(self._requests_root(), req_id, "jobs", job_id, "job_state.json")

    def _stage_log_path(self, req_id, job_id):
        return os.path.join(self._requests_root(), req_id, "jobs", job_id, "stage.log")

    # ------------------------------------------------------------
    # Active requests
    # ------------------------------------------------------------
    def list_active(self):
        root = self._requests_root()
        if not os.path.exists(root):
            return []

        active = []
        for req_id in os.listdir(root):
            req_dir = os.path.join(root, req_id)
            jobs_dir = os.path.join(req_dir, "jobs")
            if not os.path.exists(jobs_dir):
                continue

            # Look for any job without completion marker
            for job_id in os.listdir(jobs_dir):
                job_state_path = self._job_state_path(req_id, job_id)
                if not os.path.exists(job_state_path):
                    continue
                state = json.load(open(job_state_path))
                if not state.get("complete", False):
                    active.append(req_id)
                    break

        return active

    # ------------------------------------------------------------
    # Parse stage.log for progress
    # ------------------------------------------------------------
    def _parse_progress(self, log_path):
        if not os.path.exists(log_path):
            return ("unknown", None)

        lines = open(log_path).read().splitlines()
        for line in reversed(lines):
            if "[" in line and "/" in line:
                try:
                    part = line.split("]")[0].strip("[")
                    done, total = map(int, part.split("/"))
                    return (f"{done}/{total}", (done, total))
                except:
                    pass

        return ("unknown", None)

    # ------------------------------------------------------------
    # ETA estimation
    # ------------------------------------------------------------
    def _estimate_eta(self, req_id, job_id, progress_tuple):
        if progress_tuple is None:
            return "unknown"

        done, total = progress_tuple
        remaining = total - done

        # Look for timing in stage.log
        log_path = self._stage_log_path(req_id, job_id)
        lines = open(log_path).read().splitlines()

        times = []
        for line in lines:
            if "seconds" in line.lower():
                try:
                    t = float(line.split()[-2])
                    times.append(t)
                except:
                    pass

        if not times:
            return "unknown"

        avg = sum(times) / len(times)
        eta_seconds = remaining * avg

        return f"~{int(eta_seconds)}s"

    # ------------------------------------------------------------
    # Status for a single request
    # ------------------------------------------------------------
    def get_status(self, req_id):
        req_dir = os.path.join(self._requests_root(), req_id)
        jobs_dir = os.path.join(req_dir, "jobs")

        job_ids = sorted(os.listdir(jobs_dir))
        current_job = job_ids[-1]

        job_state = json.load(open(self._job_state_path(req_id, current_job)))

        progress_str, progress_tuple = self._parse_progress(
            self._stage_log_path(req_id, current_job)
        )

        eta = self._estimate_eta(req_id, current_job, progress_tuple)

        return {
            "request_id": req_id,
            "stage": job_state.get("stage"),
            "job_id": current_job,
            "progress": progress_str,
            "eta": eta,
            "resources": job_state.get("resources", {}),
        }

    # ------------------------------------------------------------
    # Human-readable summary
    # ------------------------------------------------------------
    def summary(self):
        active = self.list_active()
        if not active:
            return "No active pipeline runs."

        out = []
        for req_id in active:
            st = self.get_status(req_id)
            out.append(
                f"Request: {st['request_id']}\n"
                f"  Stage: {st['stage']}\n"
                f"  Job: {st['job_id']}\n"
                f"  Progress: {st['progress']}\n"
                f"  ETA: {st['eta']}\n"
            )

        return "\n".join(out)
