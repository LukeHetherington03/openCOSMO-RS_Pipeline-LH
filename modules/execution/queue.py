#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
from datetime import datetime, timezone

BASE_DIR = os.path.abspath("pipeline_data")
QUEUE_ROOT = os.path.join(BASE_DIR, "queue")

DIRS = {
    "pending":   os.path.join(QUEUE_ROOT, "pending"),
    "running":   os.path.join(QUEUE_ROOT, "running"),
    "completed": os.path.join(QUEUE_ROOT, "completed"),
    "failed":    os.path.join(QUEUE_ROOT, "failed"),
    "cancelled": os.path.join(QUEUE_ROOT, "cancelled"),
}

PID_PATH = os.path.join(QUEUE_ROOT, "worker.pid")
HEARTBEAT_PATH = os.path.join(QUEUE_ROOT, "worker_heartbeat")
WORKER_LOG = os.path.join(QUEUE_ROOT, "worker.log")


def now_iso():
    return datetime.now(timezone.utc).isoformat()


class QueueManager:
    """
    File-based queue manager for pipeline requests.
    """

    # ------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------
    @staticmethod
    def init():
        os.makedirs(QUEUE_ROOT, exist_ok=True)
        for d in DIRS.values():
            os.makedirs(d, exist_ok=True)

    # ------------------------------------------------------------
    # Enqueue
    # ------------------------------------------------------------
    @staticmethod
    def enqueue(request_id: str, priority: int = 0):
        QueueManager.init()

        entry = {
            "request_id": request_id,
            "priority": int(priority),
            "created_at": now_iso(),
        }

        path = os.path.join(DIRS["pending"], f"{request_id}.json")
        with open(path, "w") as f:
            json.dump(entry, f, indent=2)

        return path

    # ------------------------------------------------------------
    # Fetch next request (atomic pending → running)
    # ------------------------------------------------------------
    @staticmethod
    def next_request():
        QueueManager.init()

        pending_files = [
            os.path.join(DIRS["pending"], f)
            for f in os.listdir(DIRS["pending"])
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
                # Corrupt entry → move to failed
                bad = os.path.basename(path)
                os.rename(path, os.path.join(DIRS["failed"], bad))

        if not entries:
            return None

        # Sort by priority then created_at
        entries.sort(key=lambda x: (x[1].get("priority", 0), x[1].get("created_at", "")))

        path, data = entries[0]
        request_id = data["request_id"]

        running_path = os.path.join(DIRS["running"], f"{request_id}.json")
        os.rename(path, running_path)

        return request_id

    # ------------------------------------------------------------
    # State transitions
    # ------------------------------------------------------------
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
        path = os.path.join(DIRS["pending"], f"{request_id}.json")
        if not os.path.exists(path):
            return False

        with open(path) as f:
            data = json.load(f)

        data["priority"] = int(new_priority)

        with open(path, "w") as f:
            json.dump(data, f, indent=2)

        return True

    @staticmethod
    def _move(request_id: str, src: str, dst: str):
        src_path = os.path.join(DIRS[src], f"{request_id}.json")
        dst_path = os.path.join(DIRS[dst], f"{request_id}.json")
        if os.path.exists(src_path):
            os.rename(src_path, dst_path)

    # ------------------------------------------------------------
    # Listing / stats
    # ------------------------------------------------------------
    @staticmethod
    def list_pending():
        return sorted(f[:-5] for f in os.listdir(DIRS["pending"]) if f.endswith(".json"))

    @staticmethod
    def list_running():
        return sorted(f[:-5] for f in os.listdir(DIRS["running"]) if f.endswith(".json"))

    @staticmethod
    def list_completed():
        return sorted(f[:-5] for f in os.listdir(DIRS["completed"]) if f.endswith(".json"))

    @staticmethod
    def list_failed():
        return sorted(f[:-5] for f in os.listdir(DIRS["failed"]) if f.endswith(".json"))

    @staticmethod
    def stats():
        return {
            "pending":   len(QueueManager.list_pending()),
            "running":   len(QueueManager.list_running()),
            "completed": len(QueueManager.list_completed()),
            "failed":    len(QueueManager.list_failed()),
        }

    # ------------------------------------------------------------
    # Worker PID + heartbeat
    # ------------------------------------------------------------
    @staticmethod
    def write_pid(pid: int):
        with open(PID_PATH, "w") as f:
            f.write(str(pid))

    @staticmethod
    def read_pid():
        try:
            with open(QueueManager.pid_file) as f:
                pid_str = f.read().strip()
                return int(pid_str) if pid_str else None
        except:
            return None


    @staticmethod
    def clear_pid():
        if os.path.exists(PID_PATH):
            os.remove(PID_PATH)

    @staticmethod
    def update_heartbeat():
        with open(HEARTBEAT_PATH, "w") as f:
            f.write(now_iso())

    @staticmethod
    def read_heartbeat():
        if not os.path.exists(HEARTBEAT_PATH):
            return None
        return open(HEARTBEAT_PATH).read().strip()

    @staticmethod
    def is_worker_running():
        pid = QueueManager.read_pid()
        if not pid:
            return False

        try:
            os.kill(pid, 0)  # check if process exists
            return True
        except ProcessLookupError:
            return False
        except PermissionError:
            return True

    # ------------------------------------------------------------
    # Worker logging
    # ------------------------------------------------------------
    @staticmethod
    def log_worker(message: str):
        line = f"{now_iso()} {message}\n"
        with open(WORKER_LOG, "a") as f:
            f.write(line)
        print(f"[worker] {message}")
