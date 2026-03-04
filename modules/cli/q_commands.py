#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modules/cli/q_commands.py

Queue CLI commands.

Changes from previous version
------------------------------
stop()
  - Sends SIGTERM to the PROCESS GROUP (os.killpg) not just the worker PID.
    This propagates to pool workers and ORCA/XTB subprocesses.
  - No longer calls clear_pid() immediately — the worker clears its own
    PID/PGID files in run_forever() after it actually exits.
  - Waits briefly and reports whether the worker is still running.

kill()  (new)
  - Sends SIGKILL to the process group.  Instant.  Use when stop() is
    taking too long or you don't care about clean shutdown.
  - Cleans up PID/PGID files itself since the worker can't (SIGKILL is
    uncatchable).
"""

import json
import os
import signal
import subprocess
import sys
import time

from modules.execution.queue import QueueManager

CONFIG_PATH = os.path.abspath("config/paths.json")


def load_base_dir() -> str:
    with open(CONFIG_PATH) as f:
        cfg = json.load(f)
    return cfg["base_dir"]


class QueueCommands:

    @staticmethod
    def dispatch(args):
        if not args:
            print("Usage: pl q <start|stop|kill|status|list|cancel|reprio>")
            return

        sub = args[0]

        if sub == "start":  return QueueCommands.start()
        if sub == "stop":   return QueueCommands.stop()
        if sub == "kill":   return QueueCommands.kill()    # new
        if sub == "status": return QueueCommands.status()
        if sub == "list":   return QueueCommands.list_queue()
        if sub == "cancel": return QueueCommands.cancel(args[1:])
        if sub == "reprio": return QueueCommands.reprio(args[1:])

        print(f"Unknown queue command: {sub}")

    # ─────────────────────────────────────────────────────────────
    # Worker control
    # ─────────────────────────────────────────────────────────────

    @staticmethod
    def start():
        base_dir = load_base_dir()
        QueueManager.init_base_dir(base_dir)

        if QueueManager.is_worker_running():
            print(f"Worker already running (PID {QueueManager.read_pid()})")
            return

        subprocess.Popen(
            [sys.executable, "-m", "modules.execution.worker", base_dir],
            stdout=open(os.devnull, "w"),
            stderr=open(os.devnull, "w"),
        )
        print("Queue worker started.")

    @staticmethod
    def stop():
        """
        Graceful stop.

        Sends SIGTERM to the worker's process group.  This reaches:
          - the worker process itself
          - all ProcessPoolExecutor worker processes
          - all ORCA / XTB / gXTB child processes

        The worker's signal handler cancels pending futures and sets
        _stop_now = True.  Running futures will receive SIGTERM, terminate,
        and their items will re-run on next resume (safe — checkpoints are
        only written on success).

        The worker clears its own PID/PGID files before exiting.
        This command waits up to 15 s and reports the outcome.
        """
        base_dir = load_base_dir()
        QueueManager.init_base_dir(base_dir)

        pgid = QueueManager.read_pgid()
        pid  = QueueManager.read_pid()

        if not pgid and not pid:
            print("No worker is currently running.")
            return

        target = pgid or pid

        try:
            if pgid:
                os.killpg(pgid, signal.SIGTERM)
                print(f"SIGTERM sent to process group {pgid} "
                      f"(worker PID {pid}).")
            else:
                # Fallback: no PGID file — send to PID only
                os.kill(pid, signal.SIGTERM)
                print(f"SIGTERM sent to worker PID {pid} "
                      f"(no PGID file — child processes may persist).")
        except ProcessLookupError:
            print("Worker process not found — already stopped.")
            QueueManager.clear_pid()
            QueueManager.clear_pgid()
            return

        # Wait for the worker to exit (it clears the PID file itself)
        print("Waiting for worker to finish current item and exit...")
        for _ in range(15):
            time.sleep(1)
            if not QueueManager.is_worker_running():
                print("Worker stopped cleanly.")
                return

        print(
            "Worker has not exited after 15 s.  It may still be finishing "
            "the current item.  Use  pl q kill  to force-stop immediately."
        )

    @staticmethod
    def kill():
        """
        Immediate hard stop.

        Sends SIGKILL to the worker's process group.  Everything in the
        tree dies instantly — no cleanup code runs.  In-flight items will
        NOT have written their checkpoints, so they will re-run on resume.
        This is safe because checkpoints are only written on item success.

        Use this when  pl q stop  is taking too long or after a crash where
        zombie pool workers are still consuming cores.
        """
        base_dir = load_base_dir()
        QueueManager.init_base_dir(base_dir)

        pgid = QueueManager.read_pgid()
        pid  = QueueManager.read_pid()

        if not pgid and not pid:
            print("No worker is currently running.")
            return

        try:
            if pgid:
                os.killpg(pgid, signal.SIGKILL)
                print(f"SIGKILL sent to process group {pgid} — all processes terminated.")
            else:
                os.kill(pid, signal.SIGKILL)
                print(f"SIGKILL sent to PID {pid} — "
                      f"child processes may still be running (no PGID file).")
        except ProcessLookupError:
            print("Worker process not found — already stopped.")

        # Worker can't clean up after SIGKILL — do it here
        QueueManager.clear_pid()
        QueueManager.clear_pgid()
        print("PID/PGID files cleared.")

    # ─────────────────────────────────────────────────────────────
    # Queue inspection
    # ─────────────────────────────────────────────────────────────

    @staticmethod
    def status():
        base_dir = load_base_dir()
        QueueManager.init_base_dir(base_dir)

        pending   = len(os.listdir(QueueManager.DIRS["pending"]))
        running   = len(os.listdir(QueueManager.DIRS["running"]))
        completed = len(os.listdir(QueueManager.DIRS["completed"]))
        failed    = len(os.listdir(QueueManager.DIRS["failed"]))

        pid  = QueueManager.read_pid()
        pgid = QueueManager.read_pgid()

        print(f"Worker running : {QueueManager.is_worker_running()}"
              + (f"  (PID {pid}, PGID {pgid})" if pid else ""))
        print(f"Pending        : {pending}")
        print(f"Running        : {running}")
        print(f"Completed      : {completed}")
        print(f"Failed         : {failed}")

    @staticmethod
    def list_queue():
        base_dir = load_base_dir()
        QueueManager.init_base_dir(base_dir)

        print("Pending:")
        for f in sorted(os.listdir(QueueManager.DIRS["pending"])):
            print("  -", f[:-5])

        print("\nRunning:")
        for f in sorted(os.listdir(QueueManager.DIRS["running"])):
            print("  -", f[:-5])

    # ─────────────────────────────────────────────────────────────
    # Queue operations
    # ─────────────────────────────────────────────────────────────

    @staticmethod
    def cancel(args):
        base_dir = load_base_dir()
        QueueManager.init_base_dir(base_dir)
        rid = args[0]
        QueueManager.cancel(rid)
        print(f"Cancelled {rid}")

    @staticmethod
    def reprio(args):
        base_dir = load_base_dir()
        QueueManager.init_base_dir(base_dir)
        rid, prio = args[0], int(args[1])
        QueueManager.reprioritise(rid, prio)
        print(f"Updated priority for {rid} → {prio}")