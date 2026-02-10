#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import time
import signal

from modules.execution.queue import QueueManager, PID_PATH, HEARTBEAT_PATH
from modules.execution.runner import PipelineRunner

BASE_DIR = os.path.abspath("pipeline_data")
DRAIN_FLAG = os.path.join(os.path.abspath("pipeline_data"), "queue", "drain")


class QueueWorker:
    """
    Single-worker queue processor.
    """

    def __init__(self, base_dir: str):
        self.base_dir = os.path.abspath(base_dir)
        self._stop_now = False

    # ------------------------------------------------------------
    # Signal handling
    # ------------------------------------------------------------
    def _handle_sigterm(self, signum, frame):
        QueueManager.log_worker("SIGTERM received — stopping immediately.")
        self._stop_now = True

    # ------------------------------------------------------------
    # Main loop
    # ------------------------------------------------------------
    def run_forever(self, idle_sleep: float = 1.0):
        signal.signal(signal.SIGTERM, self._handle_sigterm)

        QueueManager.log_worker("Queue worker started.")
        while True:
            QueueManager.update_heartbeat()

            # Immediate stop
            if self._stop_now:
                QueueManager.log_worker("Worker stopping now.")
                QueueManager.clear_pid()
                return

            # Drain mode: finish current request, then exit
            if os.path.exists(DRAIN_FLAG):
                QueueManager.log_worker("Drain mode active — finishing current request only.")

            req_id = QueueManager.next_request()

            if req_id is None:
                # If draining and no request running → exit
                if os.path.exists(DRAIN_FLAG):
                    QueueManager.log_worker("Drain complete — exiting worker.")
                    os.remove(DRAIN_FLAG)
                    QueueManager.clear_pid()
                    return

                time.sleep(idle_sleep)
                continue

            QueueManager.log_worker(f"Starting request {req_id}")

            try:
                PipelineRunner.run_request(req_id, self.base_dir)
                QueueManager.mark_completed(req_id)
                QueueManager.log_worker(f"Request {req_id} completed.")
            except Exception as e:
                QueueManager.mark_failed(req_id)
                QueueManager.log_worker(f"Request {req_id} failed: {e}")

            # If draining → exit after finishing this request
            if os.path.exists(DRAIN_FLAG):
                QueueManager.log_worker("Drain mode: exiting after current request.")
                os.remove(DRAIN_FLAG)
                QueueManager.clear_pid()
                return


@staticmethod
def start_worker():
    pid = QueueManager.read_pid()

    # If PID file exists but process is gone → clear it
    if pid and not QueueManager.is_worker_running():
        QueueManager.clear_pid()

    if QueueManager.is_worker_running():
        print(f"Worker already running (PID {pid})")
        return


    pid = os.getpid()
    QueueManager.write_pid(pid)

    worker = QueueWorker(BASE_DIR)
    worker.run_forever()


if __name__ == "__main__":
    start_worker()
