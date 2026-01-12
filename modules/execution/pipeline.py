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
    Class-oriented command-line interface for pipeline control.
    Commands:
        run <request_id>
        pause [request_id]
        resume [request_id]
        stop [request_id]
        status [request_id]
        list
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
            print("Usage: pipeline run <request_id> [base_dir]")
            sys.exit(1)

        request_id = self.args[0]
        base_dir = self._base_dir(self.args[1:])
        PipelineRunner.run_request(request_id, base_dir)

    def pause(self):
        status = PipelineStatus(DEFAULT_BASE)

        if len(self.args) >= 1:
            request_id = self.args[0]
            base_dir = self._base_dir(self.args[1:])
        else:
            request_id = self._autofill_request(status)
            base_dir = DEFAULT_BASE

        PipelineControl(base_dir).pause(request_id)
        print(f"Paused request {request_id}")

    def resume(self):
        status = PipelineStatus(DEFAULT_BASE)

        if len(self.args) >= 1:
            request_id = self.args[0]
            base_dir = self._base_dir(self.args[1:])
        else:
            request_id = self._autofill_request(status)
            base_dir = DEFAULT_BASE

        PipelineControl(base_dir).resume(request_id)
        print(f"Resumed request {request_id}")

    def stop(self):
        status = PipelineStatus(DEFAULT_BASE)

        if len(self.args) >= 1:
            request_id = self.args[0]
            base_dir = self._base_dir(self.args[1:])
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

        request_id = self.args[0]
        base_dir = self._base_dir(self.args[1:])
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
    # Dispatcher
    # ------------------------------------------------------------
    def dispatch(self):
        if not self.cmd:
            print("Commands: run, pause, resume, stop, status, list")
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

        print(f"Unknown command: {self.cmd}")
        sys.exit(1)


def main():
    cli = PipelineCLI(sys.argv)
    cli.dispatch()


if __name__ == "__main__":
    main()
