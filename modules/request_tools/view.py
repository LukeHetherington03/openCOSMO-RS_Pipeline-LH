#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
from datetime import datetime, timezone


def now_iso():
    return datetime.now(timezone.utc).isoformat()


class RequestView:
    """
    Read-only view of a pipeline request.
    Loads request.json, pipeline_state.json, and job logs.
    """

    def __init__(self, base_dir, request_id, data, pipeline_state):
        self.base_dir = base_dir
        self.request_id = request_id
        self.data = data
        self.pipeline_state = self._normalise_pipeline_state(pipeline_state)

    # ------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------
    @classmethod
    def load(cls, base_dir, request_id):
        req_dir = os.path.join(base_dir, "requests", request_id)
        meta_path = os.path.join(req_dir, "request.json")
        state_path = os.path.join(req_dir, "pipeline_state.json")

        if not os.path.exists(meta_path):
            raise FileNotFoundError(f"No request.json for {request_id}")

        with open(meta_path) as f:
            data = json.load(f)

        if os.path.exists(state_path):
            with open(state_path) as f:
                pipeline_state = json.load(f)
        else:
            pipeline_state = {}

        return cls(base_dir, request_id, data, pipeline_state)

    # ------------------------------------------------------------
    # Normalisation
    # ------------------------------------------------------------
    def _normalise_pipeline_state(self, ps):
        if isinstance(ps, str):
            return {
                "state": ps,
                "stage": "n/a",
                "updated_at": "n/a",
            }

        if isinstance(ps, dict):
            ps.setdefault("state", "unknown")
            ps.setdefault("stage", "n/a")
            ps.setdefault("updated_at", "n/a")
            return ps

        return {
            "state": "unknown",
            "stage": "n/a",
            "updated_at": "n/a",
        }

    # ------------------------------------------------------------
    # CLI helpers
    # ------------------------------------------------------------
    def pretty_status(self):
        ps = self.pipeline_state
        return (
            f"Request:      {self.request_id}\n"
            f"State:        {ps['state']}\n"
            f"Stage:        {ps['stage']}\n"
            f"Updated:      {ps['updated_at']}\n"
            f"User status:  {self.data.get('state', {}).get('user_status', 'n/a')}\n"
            f"Notes:        {self.data.get('state', {}).get('notes', '')}\n"
        )

    def pretty_info(self):
        return (
            f"Request ID:     {self.request_id}\n"
            f"Parameters:     {self.data.get('parameters')}\n"
            f"Jobs:           {self.data.get('jobs')}\n"
            f"User status:    {self.data.get('state', {}).get('user_status')}\n"
            f"Pinned:         {self.data.get('state', {}).get('pinned', False)}\n"
            f"Notes:          {self.data.get('state', {}).get('notes', '')}\n"
            f"Created at:     {self.data.get('created_at')}\n"
            f"Updated at:     {self.data.get('updated_at')}\n"
        )

    def latest_logs(self):
        jobs_dir = os.path.join(self.base_dir, "requests", self.request_id, "jobs")
        if not os.path.exists(jobs_dir):
            return "No jobs directory."

        job_ids = sorted(os.listdir(jobs_dir))
        if not job_ids:
            return "No jobs found."

        latest = job_ids[-1]
        log_path = os.path.join(jobs_dir, latest, "stage.log")

        if not os.path.exists(log_path):
            return f"No stage.log for job {latest}"

        with open(log_path) as f:
            return f.read()

    # ------------------------------------------------------------
    # Listing
    # ------------------------------------------------------------
    @staticmethod
    def list_all(base_dir):
        req_root = os.path.join(base_dir, "requests")
        if not os.path.exists(req_root):
            return []

        return sorted(
            d for d in os.listdir(req_root)
            if os.path.isdir(os.path.join(req_root, d))
        )
