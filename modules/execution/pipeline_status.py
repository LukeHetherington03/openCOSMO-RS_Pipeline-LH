#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
from datetime import datetime


class PipelineStatus:
    """
    Reports status of pipeline requests.
    """

    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.requests_root = os.path.join(base_dir, "requests")

    # ------------------------------------------------------------
    # Request ID resolver (NEW)
    # ------------------------------------------------------------
    def _resolve_request_id(self, partial_id):
        """
        Resolve a partial request ID to the full folder name.
        If multiple matches exist, return None and let caller handle it.
        """
        if not os.path.exists(self.requests_root):
            return None

        matches = [f for f in os.listdir(self.requests_root) if f.startswith(partial_id)]

        if len(matches) == 1:
            return matches[0]

        if len(matches) == 0:
            return None

        # Multiple matches — caller will handle messaging
        return matches


    # ------------------------------------------------------------
    # List active requests
    # ------------------------------------------------------------
    def list_active(self):
        if not os.path.exists(self.requests_root):
            return []

        active = []
        for folder in os.listdir(self.requests_root):
            control_path = os.path.join(self.requests_root, folder, "control.json")
            if not os.path.exists(control_path):
                active.append(folder)
                continue

            with open(control_path) as f:
                ctrl = json.load(f)

            if not ctrl.get("stop", False):
                active.append(folder)

        return active

    # ------------------------------------------------------------
    # Get status for a specific request
    # ------------------------------------------------------------
    def get_status(self, request_id):
        folder = self._resolve_request_id(request_id)

        if folder is None:
            return f"No request matches prefix '{request_id}'"

        if isinstance(folder, list):
            msg = [f"Multiple requests match prefix '{request_id}':"]
            for f in folder:
                msg.append(f"  - {f}")
            msg.append("Please specify further.")
            return "\n".join(msg)

        if folder is None:
            return f"No request matches prefix '{request_id}'"

        req_dir = os.path.join(self.requests_root, folder)
        jobs_dir = os.path.join(req_dir, "jobs")

        if not os.path.exists(jobs_dir):
            return f"Request {folder}: no jobs directory (request incomplete or corrupted)"

        job_ids = sorted(os.listdir(jobs_dir))

        lines = [f"Status for request {folder}:"]
        lines.append(f"  Jobs: {len(job_ids)}")

        for job_id in job_ids:
            job_dir = os.path.join(jobs_dir, job_id)
            state_file = os.path.join(job_dir, "job_state.json")

            if os.path.exists(state_file):
                with open(state_file) as f:
                    state = json.load(f)
                stage = state.get("stage", "unknown")
                num_in = state.get("num_input", "?")
                num_out = state.get("num_output", "?")
                lines.append(f"    - {job_id}: stage={stage}, {num_in}→{num_out}")
            else:
                lines.append(f"    - {job_id}: (no job_state.json)")

        return "\n".join(lines)

    # ------------------------------------------------------------
    # Summary of all active requests
    # ------------------------------------------------------------
    def summary(self):
        active = self.list_active()
        if not active:
            return "No active requests."

        lines = ["Active requests:"]
        for r in active:
            lines.append(f"  - {r}")
        return "\n".join(lines)
