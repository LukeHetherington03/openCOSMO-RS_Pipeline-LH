#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
from datetime import datetime, timezone

class QueueManager:
    """
    File-based queue manager for pipeline requests.
    BASE_DIR is loaded once from paths.json and stored internally.
    """

    BASE_DIR = None
    QUEUE_ROOT = None
    DIRS = None
    PID_PATH = None
    HEARTBEAT_PATH = None
    WORKER_LOG = None

    # ------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------
    @staticmethod
    def init_base_dir(base_dir: str):
        """Initialise all queue paths from a single base directory."""
        QueueManager.BASE_DIR = os.path.abspath(base_dir)
        QueueManager.QUEUE_ROOT = os.path.join(QueueManager.BASE_DIR, "queue")

        QueueManager.DIRS = {
            "pending":   os.path.join(QueueManager.QUEUE_ROOT, "pending"),
            "running":   os.path.join(QueueManager.QUEUE_ROOT, "running"),
            "completed": os.path.join(QueueManager.QUEUE_ROOT, "completed"),
            "failed":    os.path.join(QueueManager.QUEUE_ROOT, "failed"),
            "cancelled": os.path.join(QueueManager.QUEUE_ROOT, "cancelled"),
        }

        QueueManager.PID_PATH = os.path.join(QueueManager.QUEUE_ROOT, "worker.pid")
        QueueManager.HEARTBEAT_PATH = os.path.join(QueueManager.QUEUE_ROOT, "worker_heartbeat")
        QueueManager.WORKER_LOG = os.path.join(QueueManager.QUEUE_ROOT, "worker.log")

        os.makedirs(QueueManager.QUEUE_ROOT, exist_ok=True)
        for d in QueueManager.DIRS.values():
            os.makedirs(d, exist_ok=True)

    # ------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------
    @staticmethod
    def now_iso():
        return datetime.now(timezone.utc).isoformat()

    # ------------------------------------------------------------
    # Enqueue
    # ------------------------------------------------------------
    @staticmethod
    def enqueue(request_id: str, priority: int = 0):
        entry = {
            "request_id": request_id,
            "priority": int(priority),
            "created_at": QueueManager.now_iso(),
        }

        path = os.path.join(QueueManager.DIRS["pending"], f"{request_id}.json")
        with open(path, "w") as f:
            json.dump(entry, f, indent=2)

        return path

    # ------------------------------------------------------------
    # Fetch next request
    # ------------------------------------------------------------
    @staticmethod
    def next_request():
        pending_files = [
            os.path.join(QueueManager.DIRS["pending"], f)
            for f in os.listdir(QueueManager.DIRS["pending"])
            if f.endswith(".json")
        ]

        if not pending_files:
            return None

        entries = []
        for path in pending_files:
            try:
                with open(path) as f:
                    data = json.load(f)
                entries.append((path, data))
            except Exception:
                bad = os.path.basename(path)
                os.rename(path, os.path.join(QueueManager.DIRS["failed"], bad))

        if not entries:
            return None

        entries.sort(key=lambda x: (x[1].get("priority", 0), x[1].get("created_at", "")))

        path, data = entries[0]
        request_id = data["request_id"]

        running_path = os.path.join(QueueManager.DIRS["running"], f"{request_id}.json")
        os.rename(path, running_path)

        return request_id

    # ------------------------------------------------------------
    # State transitions
    # ------------------------------------------------------------
    @staticmethod
    def _move(request_id: str, src: str, dst: str):
        src_path = os.path.join(QueueManager.DIRS[src], f"{request_id}.json")
        dst_path = os.path.join(QueueManager.DIRS[dst], f"{request_id}.json")
        if os.path.exists(src_path):
            os.rename(src_path, dst_path)

    @staticmethod
    def mark_completed(request_id: str):
        QueueManager._move(request_id, "running", "completed")

    @staticmethod
    def mark_failed(request_id: str):
        QueueManager._move(request_id, "running", "failed")

    @staticmethod
    def cancel(request_id: str):
        QueueManager._move(request_id, "pending", "cancelled")

    @staticmethod
    def reprioritise(request_id: str, new_priority: int):
        path = os.path.join(QueueManager.DIRS["pending"], f"{request_id}.json")
        if not os.path.exists(path):
            return False

        with open(path) as f:
            data = json.load(f)

        data["priority"] = int(new_priority)

        with open(path, "w") as f:
            json.dump(data, f, indent=2)

        return True

    # ------------------------------------------------------------
    # PID + heartbeat
    # ------------------------------------------------------------
    @staticmethod
    def write_pid(pid: int):
        with open(QueueManager.PID_PATH, "w") as f:
            f.write(str(pid))

    @staticmethod
    def read_pid():
        try:
            with open(QueueManager.PID_PATH) as f:
                pid_str = f.read().strip()
                return int(pid_str) if pid_str else None
        except:
            return None

    @staticmethod
    def clear_pid():
        if os.path.exists(QueueManager.PID_PATH):
            os.remove(QueueManager.PID_PATH)

    @staticmethod
    def update_heartbeat():
        with open(QueueManager.HEARTBEAT_PATH, "w") as f:
            f.write(QueueManager.now_iso())

    @staticmethod
    def read_heartbeat():
        if not os.path.exists(QueueManager.HEARTBEAT_PATH):
            return None
        return open(QueueManager.HEARTBEAT_PATH).read().strip()

    @staticmethod
    def is_worker_running():
        pid = QueueManager.read_pid()
        if not pid:
            return False

        try:
            os.kill(pid, 0)
            return True
        except ProcessLookupError:
            return False
        except PermissionError:
            return True

    # ------------------------------------------------------------
    # Logging
    # ------------------------------------------------------------
    @staticmethod
    def log_worker(message: str):
        line = f"{QueueManager.now_iso()} {message}\n"
        with open(QueueManager.WORKER_LOG, "a") as f:
            f.write(line)
        print(f"[worker] {message}")
