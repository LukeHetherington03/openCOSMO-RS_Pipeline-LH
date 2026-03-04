#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modules/execution/queue.py

Queue manager for the pipeline worker system.

Manages the filesystem-based queue used to schedule, track, and control
pipeline requests.  The queue is a set of directories under the base
queue path:

    queue/
        pending/    requests waiting to run  (<request_id>.json)
        running/    request currently executing  (<request_id>.json)
        completed/  finished requests
        failed/     requests that raised an unhandled exception
        worker.pid  PID of the running worker process
        worker.pgid process group ID of the running worker  [new 03/03/2026]
        worker.log  append-only log of worker events
        heartbeat   timestamp updated every idle loop by the worker

One request runs at a time.  The worker picks the highest-priority
pending request, moves it to running/, and blocks until it finishes.

Priority is an integer stored in the queue entry JSON.  Lower values
run first.  Default priority is 100.
"""

import json
import os
import time


class QueueManager:
    """
    Static interface to the filesystem queue.

    All methods are class methods.  Call init_base_dir() once at
    process startup before using any other method.

    Directory layout
    ----------------
    DIRS["pending"]    queue/pending/
    DIRS["running"]    queue/running/
    DIRS["completed"]  queue/completed/
    DIRS["failed"]     queue/failed/

    Control files
    -------------
    PID_FILE           queue/worker.pid     — worker process ID
    PGID_FILE          queue/worker.pgid    — worker process group ID  [new 03/03/2026]
    LOG_FILE           queue/worker.log     — append-only worker events
    HEARTBEAT_FILE     queue/heartbeat      — last-alive timestamp

    Queue entry format  (each .json file)
    --------------------------------------
    {
        "request_id": "abc123",
        "priority":   100,
        "submitted":  "2026-01-01T00:00:00Z"
    }
    """

    DIRS: dict         = {}
    PID_FILE: str      = None
    PGID_FILE: str     = None   # [new 03/03/2026]
    LOG_FILE: str      = None
    HEARTBEAT_FILE:str = None

    # -------------------------------------------------------------------------
    # Initialisation
    # -------------------------------------------------------------------------

    @classmethod
    def init_base_dir(cls, base_dir: str):
        """
        Set the base queue directory and create all subdirectories.

        Must be called once at worker startup and once in each CLI command
        before any other QueueManager method is used.

        Parameters
        ----------
        base_dir : str
            Root directory for all pipeline data.  The queue lives at
            <base_dir>/queue/.
        """
        queue_dir = os.path.join(base_dir, "queue")

        cls.DIRS = {
            "pending":   os.path.join(queue_dir, "pending"),
            "running":   os.path.join(queue_dir, "running"),
            "completed": os.path.join(queue_dir, "completed"),
            "failed":    os.path.join(queue_dir, "failed"),
        }

        cls.PID_FILE       = os.path.join(queue_dir, "worker.pid")
        cls.PGID_FILE      = os.path.join(queue_dir, "worker.pgid")   # [new 03/03/2026]
        cls.LOG_FILE       = os.path.join(queue_dir, "worker.log")
        cls.HEARTBEAT_FILE = os.path.join(queue_dir, "heartbeat")

        for d in cls.DIRS.values():
            os.makedirs(d, exist_ok=True)

    # -------------------------------------------------------------------------
    # Queue entry operations
    # -------------------------------------------------------------------------

    @classmethod
    def enqueue(cls, request_id: str, priority: int = 100):
        """
        Add a request to the pending queue.

        Parameters
        ----------
        request_id : str
            Unique identifier for the request.
        priority : int
            Scheduling priority.  Lower values run first.  Default 100.
        """
        entry = {
            "request_id": request_id,
            "priority":   priority,
            "submitted":  _now(),
        }
        path = os.path.join(cls.DIRS["pending"], f"{request_id}.json")
        with open(path, "w") as f:
            json.dump(entry, f, indent=2)

    @classmethod
    def next_request(cls) -> str | None:
        """
        Return the request_id of the highest-priority pending request,
        or None if the queue is empty.

        Moves the entry from pending/ to running/ atomically.
        """
        pending = cls.DIRS["pending"]
        entries = []

        for fname in os.listdir(pending):
            if not fname.endswith(".json"):
                continue
            path = os.path.join(pending, fname)
            try:
                with open(path) as f:
                    entry = json.load(f)
                entries.append((entry.get("priority", 100), fname, entry))
            except Exception:
                continue

        if not entries:
            return None

        entries.sort(key=lambda x: x[0])
        _, fname, entry = entries[0]

        src = os.path.join(pending, fname)
        dst = os.path.join(cls.DIRS["running"], fname)
        os.rename(src, dst)

        return entry["request_id"]

    @classmethod
    def mark_completed(cls, request_id: str):
        """Move a running request to the completed directory."""
        cls._move_request(request_id, "running", "completed")

    @classmethod
    def mark_failed(cls, request_id: str):
        """Move a running request to the failed directory."""
        cls._move_request(request_id, "running", "failed")

    @classmethod
    def cancel(cls, request_id: str):
        """
        Remove a pending request from the queue.

        Does nothing if the request is already running, completed, or
        not found — cancelling a running request requires pl q stop.
        """
        path = os.path.join(cls.DIRS["pending"], f"{request_id}.json")
        if os.path.exists(path):
            os.remove(path)

    @classmethod
    def reprioritise(cls, request_id: str, priority: int):
        """
        Update the priority of a pending request.

        Has no effect if the request is not in the pending queue.
        """
        path = os.path.join(cls.DIRS["pending"], f"{request_id}.json")
        if not os.path.exists(path):
            return
        with open(path) as f:
            entry = json.load(f)
        entry["priority"] = priority
        with open(path, "w") as f:
            json.dump(entry, f, indent=2)

    @classmethod
    def _move_request(cls, request_id: str, src_dir: str, dst_dir: str):
        src = os.path.join(cls.DIRS[src_dir], f"{request_id}.json")
        dst = os.path.join(cls.DIRS[dst_dir], f"{request_id}.json")
        if os.path.exists(src):
            os.rename(src, dst)

    # -------------------------------------------------------------------------
    # Worker PID
    # -------------------------------------------------------------------------

    @classmethod
    def write_pid(cls, pid: int):
        """Write the worker process ID to queue/worker.pid."""
        with open(cls.PID_FILE, "w") as f:
            f.write(str(pid))

    @classmethod
    def read_pid(cls) -> int | None:
        """
        Read the worker PID from queue/worker.pid.

        Returns None if the file is absent or unreadable.
        """
        if cls.PID_FILE and os.path.exists(cls.PID_FILE):
            try:
                return int(open(cls.PID_FILE).read().strip())
            except (ValueError, OSError):
                return None
        return None

    @classmethod
    def clear_pid(cls):
        """Remove queue/worker.pid."""
        if cls.PID_FILE and os.path.exists(cls.PID_FILE):
            try:
                os.remove(cls.PID_FILE)
            except OSError:
                pass

    @classmethod
    def is_worker_running(cls) -> bool:
        """
        Return True if a worker process is currently live.

        Reads the PID file and checks whether that process exists.
        Returns False if no PID file exists or the process is gone.
        """
        pid = cls.read_pid()
        if pid is None:
            return False
        try:
            os.kill(pid, 0)   # signal 0 — existence check only
            return True
        except ProcessLookupError:
            return False
        except PermissionError:
            return True       # process exists but we can't signal it

    # -------------------------------------------------------------------------
    # Worker PGID  [new 03/03/2026]
    # -------------------------------------------------------------------------
    #
    # The process group ID is used by pl q stop and pl q kill to send signals
    # to the entire worker subtree (pool workers, ORCA, XTB) in one call.
    #
    # The worker writes its PGID on startup after calling os.setpgrp() to
    # place itself in a new process group.  All child processes it spawns
    # automatically inherit the same PGID.
    #
    # pl q stop  — os.killpg(pgid, SIGTERM)  — graceful, propagates to tree
    # pl q kill  — os.killpg(pgid, SIGKILL)  — immediate, propagates to tree

    @classmethod
    def write_pgid(cls, pgid: int):
        """Write the worker process group ID to queue/worker.pgid."""
        with open(cls.PGID_FILE, "w") as f:
            f.write(str(pgid))

    @classmethod
    def read_pgid(cls) -> int | None:
        """
        Read the worker PGID from queue/worker.pgid.

        Returns None if the file is absent or unreadable.
        """
        if cls.PGID_FILE and os.path.exists(cls.PGID_FILE):
            try:
                return int(open(cls.PGID_FILE).read().strip())
            except (ValueError, OSError):
                return None
        return None

    @classmethod
    def clear_pgid(cls):
        """Remove queue/worker.pgid."""
        if cls.PGID_FILE and os.path.exists(cls.PGID_FILE):
            try:
                os.remove(cls.PGID_FILE)
            except OSError:
                pass

    # -------------------------------------------------------------------------
    # Heartbeat
    # -------------------------------------------------------------------------

    @classmethod
    def update_heartbeat(cls):
        """
        Write the current timestamp to queue/heartbeat.

        Called by the worker on every idle loop iteration.  Can be used
        externally to detect a stalled worker (heartbeat older than N seconds
        while the worker is supposedly running).
        """
        with open(cls.HEARTBEAT_FILE, "w") as f:
            f.write(_now())

    @classmethod
    def read_heartbeat(cls) -> str | None:
        """Return the last heartbeat timestamp string, or None if absent."""
        if cls.HEARTBEAT_FILE and os.path.exists(cls.HEARTBEAT_FILE):
            try:
                return open(cls.HEARTBEAT_FILE).read().strip()
            except OSError:
                return None
        return None

    # -------------------------------------------------------------------------
    # Worker log
    # -------------------------------------------------------------------------

    @classmethod
    def log_worker(cls, message: str):
        """
        Append a timestamped message to queue/worker.log.

        This is the worker's own event log — separate from per-request
        stage logs.  Written by the worker and by signal handlers.
        """
        line = f"[{_now()}] {message}\n"
        with open(cls.LOG_FILE, "a") as f:
            f.write(line)


# -----------------------------------------------------------------------------
# Module helpers
# -----------------------------------------------------------------------------

def _now() -> str:
    """Return current UTC time as an ISO 8601 string."""
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).isoformat(timespec="seconds")