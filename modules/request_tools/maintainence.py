#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
from datetime import datetime, timezone


def now_iso():
    return datetime.now(timezone.utc).isoformat()


class RequestAnnotation:
    """
    Mutates user-facing metadata in request.json:
      - notes
      - tags
      - pinned
      - user_status (published, archived, trash)
    """

    def __init__(self, base_dir):
        self.base_dir = base_dir

    # ------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------
    def _load(self, request_id):
        path = os.path.join(self.base_dir, "requests", request_id, "request.json")
        if not os.path.exists(path):
            raise FileNotFoundError(f"No request.json for {request_id}")

        with open(path) as f:
            return path, json.load(f)

    def _save(self, path, data):
        data.setdefault("state", {})
        data["state"]["updated_at"] = now_iso()

        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    # ------------------------------------------------------------
    # Annotation operations
    # ------------------------------------------------------------
    def note(self, request_id, text):
        path, data = self._load(request_id)
        data.setdefault("state", {})
        data["state"]["notes"] = text
        self._save(path, data)

    def tag(self, request_id, tags):
        path, data = self._load(request_id)
        data.setdefault("state", {})
        existing = data["state"].get("tags", [])
        data["state"]["tags"] = sorted(set(existing + tags))
        self._save(path, data)

    def pin(self, request_id):
        path, data = self._load(request_id)
        data.setdefault("state", {})
        data["state"]["pinned"] = True
        self._save(path, data)

    def publish(self, request_id):
        path, data = self._load(request_id)
        data.setdefault("state", {})
        data["state"]["user_status"] = "published"
        self._save(path, data)

    def archive(self, request_id):
        path, data = self._load(request_id)
        data.setdefault("state", {})
        data["state"]["user_status"] = "archived"
        self._save(path, data)

    def trash(self, request_id):
        path, data = self._load(request_id)
        data.setdefault("state", {})
        data["state"]["user_status"] = "trash"
        self._save(path, data)
