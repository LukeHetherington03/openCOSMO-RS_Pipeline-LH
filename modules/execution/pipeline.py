#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os

from modules.execution.pipeline_runner import PipelineRunner
from modules.execution.pipeline_control import PipelineControl
from modules.execution.pipeline_status import PipelineStatus

DEFAULT_BASE = os.path.abspath("pipeline_data")


class PipelineCLI:
    """
    Commands:
        run <request_id>
        pause [request_id]
        resume [request_id]
        stop [request_id]
        status [request_id]
        list
        logs <request_id>
        clean <request_id>
    """

    def __init__(self, argv):
        self.argv = argv
        self.cmd = argv[1] if len(argv) > 1 else None
        self.args = argv[2:]

    # ------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------
    def _base_dir(self, args):
        return os.path.abspath(args[0]) if args else DEFAULT_BASE

    def _resolve_request(self, partial_id, base_dir):
        status = PipelineStatus(base_dir)
        result = status._resolve_request_id(partial_id)

        if result is None:
            print(f"No request matches prefix '{partial_id}'")
            sys.exit(1)

        if isinstance(result, list):
            print(f"Multiple requests match prefix '{partial_id}':")
            for r in result:
                print(f"  - {r}")
            print("Please specify further.")
            sys.exit(1)

        return result

    def _autofill_request(self, status: PipelineStatus):
        active = status.list_active()

        if not active:
            print("No active requests found.")
            sys.exit(1)

        if len(active) == 1:
            return active[0]

        print("Multiple active requests detected:")
        for r in active:
            print(f"  - {r}")
        print("Please specify a request_id.")
        sys.exit(1)

    # ------------------------------------------------------------
    # Command handlers
    # ------------------------------------------------------------
    def run(self):
        if len(self.args) < 1:
            print("Usage: pl run <request_id> [base_dir]")
            sys.exit(1)

        partial = self.args[0]
        base_dir = self._base_dir(self.args[1:])
        request_id = self._resolve_request(partial, base_dir)

        PipelineRunner.run_request(request_id, base_dir)

    def pause(self):
        status = PipelineStatus(DEFAULT_BASE)

        if len(self.args) >= 1:
            partial = self.args[0]
            base_dir = self._base_dir(self.args[1:])
            request_id = self._resolve_request(partial, base_dir)
        else:
            request_id = self._autofill_request(status)
            base_dir = DEFAULT_BASE

        PipelineControl(base_dir).pause(request_id)
        print(f"Paused request {request_id}")

    def resume(self):
        status = PipelineStatus(DEFAULT_BASE)

        if len(self.args) >= 1:
            partial = self.args[0]
            base_dir = self._base_dir(self.args[1:])
            request_id = self._resolve_request(partial, base_dir)
        else:
            request_id = self._autofill_request(status)
            base_dir = DEFAULT_BASE

        PipelineControl(base_dir).resume(request_id)
        print(f"Resumed request {request_id}")

    def stop(self):
        status = PipelineStatus(DEFAULT_BASE)

        if len(self.args) >= 1:
            partial = self.args[0]
            base_dir = self._base_dir(self.args[1:])
            request_id = self._resolve_request(partial, base_dir)
        else:
            request_id = self._autofill_request(status)
            base_dir = DEFAULT_BASE

        PipelineControl(base_dir).stop(request_id)
        print(f"Stopped request {request_id}")

    def status(self):
        status = PipelineStatus(DEFAULT_BASE)

        if len(self.args) == 0:
            print(status.summary())
            return

        partial = self.args[0]
        base_dir = self._base_dir(self.args[1:])
        request_id = self._resolve_request(partial, base_dir)

        print(status.get_status(request_id))

    def list(self):
        base_dir = self._base_dir(self.args)
        status = PipelineStatus(base_dir)
        active = status.list_active()

        if not active:
            print("No active pipelines.")
        else:
            print("Active pipelines:")
            for r in active:
                print(f"  - {r}")

    # ------------------------------------------------------------
    # NEW: pl logs <request_id>
    # ------------------------------------------------------------
    def logs(self):
        if len(self.args) < 1:
            print("Usage: pl logs <request_id>")
            sys.exit(1)

        partial = self.args[0]
        base_dir = self._base_dir(self.args[1:])
        request_id = self._resolve_request(partial, base_dir)

        req_dir = os.path.join(base_dir, "requests", request_id)
        jobs_dir = os.path.join(req_dir, "jobs")

        if not os.path.exists(jobs_dir):
            print(f"No jobs directory for request {request_id}")
            sys.exit(1)

        job_ids = sorted(os.listdir(jobs_dir))
        if not job_ids:
            print(f"No jobs found for request {request_id}")
            sys.exit(1)

        latest_job = job_ids[-1]
        log_path = os.path.join(jobs_dir, latest_job, "stage.log")

        if not os.path.exists(log_path):
            print(f"No stage.log found for job {latest_job}")
            sys.exit(1)

        print(f"--- Logs for {request_id} / {latest_job} ---\n")
        with open(log_path) as f:
            print(f.read())

    # ------------------------------------------------------------
    # NEW: pl clean <request_id>
    # ------------------------------------------------------------
    def clean(self):
        if len(self.args) < 1:
            print("Usage: pl clean <request_id>")
            sys.exit(1)

        partial = self.args[0]
        base_dir = self._base_dir(self.args[1:])
        request_id = self._resolve_request(partial, base_dir)

        req_dir = os.path.join(base_dir, "requests", request_id)

        if not os.path.exists(req_dir):
            print(f"Request directory not found: {req_dir}")
            sys.exit(1)

        # Safety: do not delete active requests
        status = PipelineStatus(base_dir)
        if request_id in status.list_active():
            print(f"Cannot clean active request {request_id}. Stop it first.")
            sys.exit(1)

        import shutil
        shutil.rmtree(req_dir)

        print(f"Cleaned request {request_id}")

    # ------------------------------------------------------------
    # Dispatcher
    # ------------------------------------------------------------
    def dispatch(self):
        if not self.cmd:
            print("Commands: run, pause, resume, stop, status, list, logs, clean")
            sys.exit(1)

        if self.cmd == "run":
            return self.run()
        if self.cmd == "pause":
            return self.pause()
        if self.cmd == "resume":
            return self.resume()
        if self.cmd == "stop":
            return self.stop()
        if self.cmd == "status":
            return self.status()
        if self.cmd == "list":
            return self.list()
        if self.cmd == "logs":
            return self.logs()
        if self.cmd == "clean":
            return self.clean()

        print(f"Unknown command: {self.cmd}")
        sys.exit(1)


def main():
    cli = PipelineCLI(sys.argv)
    cli.dispatch()


if __name__ == "__main__":
    main()
