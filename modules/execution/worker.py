#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import time
import signal

from modules.execution.queue import QueueManager
from modules.execution.runner import PipelineRunner

class QueueWorker:
    """
    Single-worker queue processor.
    """

    def __init__(self, base_dir: str):
        self.base_dir = os.path.abspath(base_dir)
        self._stop_now = False

    def _handle_sigterm(self, signum, frame):
        QueueManager.log_worker("SIGTERM received — stopping immediately.")
        self._stop_now = True

    def run_forever(self, idle_sleep: float = 1.0):
        signal.signal(signal.SIGTERM, self._handle_sigterm)

        QueueManager.log_worker("Queue worker started.")
        while True:
            QueueManager.update_heartbeat()

            if self._stop_now:
                QueueManager.log_worker("Worker stopping now.")
                QueueManager.clear_pid()
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
            except Exception as e:
                QueueManager.mark_failed(req_id)
                QueueManager.log_worker(f"Request {req_id} failed: {e}")


def start_worker(base_dir: str):
    QueueManager.init_base_dir(base_dir)

    pid = os.getpid()
    QueueManager.write_pid(pid)

    worker = QueueWorker(base_dir)
    worker.run_forever()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python -m modules.execution.worker <base_dir>")
        sys.exit(1)

    start_worker(sys.argv[1])
