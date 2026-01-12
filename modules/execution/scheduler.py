"""
Scheduler
---------
Ensures that only one pipeline runs at a time.
Prevents resource collisions on a single shared machine.

Implements:
  - Lock file scheduling
  - Optional process-based detection
"""

import os
import psutil
import time


class Scheduler:
    LOCK_FILENAME = ".pipeline_lock"

    @staticmethod
    def lock_path(base_dir):
        return os.path.join(base_dir, Scheduler.LOCK_FILENAME)

    # ------------------------------------------------------------
    # Lock file scheduling
    # ------------------------------------------------------------
    @staticmethod
    def acquire_lock(base_dir, wait=True, poll_interval=5):
        """
        Acquire the pipeline lock.
        If wait=True, block until the lock becomes available.
        """
        lock = Scheduler.lock_path(base_dir)

        while True:
            if not os.path.exists(lock):
                # Create lock file
                with open(lock, "w") as f:
                    f.write(str(os.getpid()))
                return True

            if not wait:
                return False

            time.sleep(poll_interval)

    @staticmethod
    def release_lock(base_dir):
        """Remove the lock file."""
        lock = Scheduler.lock_path(base_dir)
        if os.path.exists(lock):
            os.remove(lock)

    # ------------------------------------------------------------
    # Optional: detect running pipeline processes
    # ------------------------------------------------------------
    @staticmethod
    def count_running_pipelines():
        """
        Count active pipeline_runner processes.
        Useful for diagnostics or future multi-pipeline scheduling.
        """
        count = 0
        for p in psutil.process_iter(['cmdline']):
            try:
                cmd = p.info['cmdline']
                if cmd and "modules.pipeline_runner" in " ".join(cmd):
                    count += 1
            except psutil.NoSuchProcess:
                pass
        return count
