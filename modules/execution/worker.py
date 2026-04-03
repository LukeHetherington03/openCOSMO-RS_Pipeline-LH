#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modules/execution/worker.py

Queue worker — runs one request at a time, sequentially.

=====================================================================
SHUTDOWN DESIGN
=====================================================================

Process group
-------------
On startup the worker calls os.setpgrp() to place itself in a new
process group (PGID = its own PID).  All child processes it spawns
(ProcessPoolExecutor workers, ORCA subprocesses) automatically
inherit this group.

This means a single os.killpg(pgid, signal) reaches the entire
tree: worker, pool workers, ORCA, XTB — everything.

PGID file
---------
The PGID is written alongside the PID file so q_commands.py can
target the group without having to read /proc or use psutil.

  queue/worker.pid   — worker PID  (existing)
  queue/worker.pgid  — process group ID  (new)

pl q stop  — graceful
---------
Sends SIGTERM to the process group.

  1. Worker's _handle_sigterm fires.
  2. Sets self._stop_now = True.
  3. Calls _cancel_active_executor() from base_stage — cancels
     pending futures, running futures finish and write checkpoints.
  4. PipelineRunner detects the interrupted state (pool returns
     cleanly with interrupted=True) and marks the job as
     "interrupted" rather than "failed", so resume works.
  5. Worker exits the run_forever loop at next opportunity.

Child ORCA/XTB processes that are already running will receive
SIGTERM (because they're in the same process group) and will
terminate.  Their items did NOT write a checkpoint so they will
re-run on next resume — this is safe by design.

pl q kill  — immediate
---------
Sends SIGKILL to the process group.  Instant death.  Same resume
safety applies — incomplete checkpoints are simply absent.

No cleanup code runs on SIGKILL — that's intentional.  The job
state file will still show "running"; PipelineRunner._find_resume_index
handles this by treating any non-"completed" job as resumable.

=====================================================================
"""

import os
import sys
import time
import signal

from modules.execution.queue import QueueManager
from modules.execution.runner import PipelineRunner
import modules.stages.base_stage as _base_stage_module


class QueueWorker:
    """
    Single-worker queue processor.

    Runs in its own process group.  SIGTERM cancels any active
    ProcessPoolExecutor and lets the current stage finish cleanly
    before the worker exits.
    """

    def __init__(self, base_dir: str, verbose: bool = False):
        self.base_dir  = os.path.abspath(base_dir)
        self._stop_now = False
        self._verbose  = verbose

    # ─────────────────────────────────────────────────────────────
    # Signal handler
    # ─────────────────────────────────────────────────────────────

    def _handle_sigterm(self, signum, frame):
        """
        SIGTERM handler.

        1. Sets _stop_now so the run_forever loop exits after the
           current request finishes (or is interrupted).
        2. Calls _cancel_active_executor() — cancels pending futures
           in any live ProcessPoolExecutor so no new ORCA/XTB jobs
           are submitted.  Running futures receive SIGTERM via the
           process group and will terminate; their items will re-run
           on resume.
        """
        QueueManager.log_worker(
            "SIGTERM received — cancelling active pool and stopping."
        )
        self._stop_now = True
        cancelled = _base_stage_module._cancel_active_executor()
        if cancelled:
            QueueManager.log_worker(
                "Active executor cancelled — running items will terminate "
                "via process group signal.  Checkpointed items are safe."
            )

    # ─────────────────────────────────────────────────────────────
    # Main loop
    # ─────────────────────────────────────────────────────────────

    def run_forever(self, idle_sleep: float = 1.0):
        import modules.execution.runner as _runner
        _runner._VERBOSE = self._verbose

        # Register signal handler
        signal.signal(signal.SIGTERM, self._handle_sigterm)

        QueueManager.log_worker("Queue worker started.")

        while True:
            QueueManager.update_heartbeat()

            if self._stop_now:
                QueueManager.log_worker("Worker stopping — exiting loop.")
                QueueManager.clear_pid()
                QueueManager.clear_pgid()
                return

            req_id = QueueManager.next_request()

            if req_id is None:
                time.sleep(idle_sleep)
                continue

            QueueManager.log_worker(f"Starting request {req_id}")

            try:
                PipelineRunner.run_request(req_id, self.base_dir)
                QueueManager.mark_completed(req_id)
                QueueManager.log_worker(f"Request {req_id} completed.")
            except InterruptedError:
                # Clean shutdown mid-request — leave in "running" state
                # so it resumes correctly next time.
                QueueManager.log_worker(
                    f"Request {req_id} interrupted — will resume on next start."
                )
                # Do NOT mark failed.  _stop_now will be True so the loop
                # exits cleanly on the next iteration.
            except Exception as e:
                QueueManager.mark_failed(req_id)
                QueueManager.log_worker(f"Request {req_id} failed: {e}")


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def start_worker(base_dir: str, verbose: bool = False):
    # ── New process group ─────────────────────────────────────────────────
    # Place this process in its own process group so SIGTERM/SIGKILL to
    # the group kills the entire subtree (pool workers, ORCA, XTB).
    os.setpgrp()

    QueueManager.init_base_dir(base_dir)

    pid  = os.getpid()
    pgid = os.getpgrp()

    QueueManager.write_pid(pid)
    QueueManager.write_pgid(pgid)   # new — needed by pl q stop/kill

    worker = QueueWorker(base_dir, verbose=verbose)
    worker.run_forever()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Pipeline queue worker")
    parser.add_argument("base_dir")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print job progress to stdout")
    args = parser.parse_args()
    start_worker(args.base_dir, verbose=args.verbose)