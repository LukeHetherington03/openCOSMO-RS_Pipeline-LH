#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modules/cli/cli.py

Entry point for the  pl  command-line interface.

Usage
-----
    pl q  <subcommand>    Queue control and inspection
    pl r  <subcommand>    Request management
    pl env <subcommand>   Environment validation

All subcommand routing is delegated to the relevant Commands class.
See each module for full subcommand documentation.
"""

import sys

from modules.cli.q_commands import QueueCommands
from modules.cli.r_commands import RequestCommands
from modules.cli.env_commands import EnvCommands


class PipelineCLI:

    def __init__(self, argv):
        self.argv = argv
        self.cmd  = argv[1] if len(argv) > 1 else None
        self.args = argv[2:]

    def dispatch(self):
        if not self.cmd:
            return self._help()

        if self.cmd == "q":
            return QueueCommands.dispatch(self.args)

        if self.cmd == "r":
            return RequestCommands.dispatch(self.args)

        if self.cmd == "env":
            return EnvCommands.dispatch(self.args)

        print(f"Unknown command: {self.cmd}")
        self._help()

    def _help(self):
        print("""
Usage: pl <command>

Queue:
    pl q start                  Start the queue worker
    pl q stop                   Graceful stop — cancels pending work,
                                  lets running items finish checkpoints,
                                  propagates SIGTERM to all child processes
    pl q kill                   Immediate stop — SIGKILL to entire process
                                  tree; in-flight items re-run on resume
    pl q status                 Show worker state and queue counts
    pl q list                   List pending and running requests
    pl q cancel <id>            Remove a pending request from the queue
    pl q reprio <id> <prio>     Change priority of a pending request

Requests:
    pl r list                   List all requests
    pl r status <id>            Show pipeline state for a request
    pl r info <id>              Detailed request information
    pl r logs <id>              Print stage logs
    pl r note <id> "text"       Attach a note to a request
    pl r tag <id> tag1 tag2     Tag a request
    pl r pin <id>               Pin a request
    pl r publish <id>           Mark a request as published
    pl r archive <id>           Archive a request
    pl r trash <id>             Move a request to trash
    pl r advance [<id>]         Force current pool to drain — SIGTERMs all running
                                  COSMO-RS subprocesses so they fail cleanly and the
                                  stage completes. Requires <id> to skip confirmation.

Environment:
    pl env check                Run full environment validation
    pl env software             Validate external executables
    pl env resources            Validate openCOSMO paths + constant files
    pl env chemistry            Validate chemistry JSON files
    pl env pip                  Validate required Python packages
    pl env table                Pretty table output
""")


def main():
    PipelineCLI(sys.argv).dispatch()


if __name__ == "__main__":
    main()