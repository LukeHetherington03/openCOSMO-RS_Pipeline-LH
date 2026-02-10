#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

from modules.cli.q_commands import QueueCommands
from modules.cli.r_commands import RequestCommands


class PipelineCLI:

    def __init__(self, argv):
        self.argv = argv
        self.cmd = argv[1] if len(argv) > 1 else None
        self.args = argv[2:]

    def dispatch(self):
        if not self.cmd:
            return self._help()

        if self.cmd == "q":
            return QueueCommands.dispatch(self.args)

        if self.cmd == "r":
            return RequestCommands.dispatch(self.args)

        print(f"Unknown command: {self.cmd}")
        self._help()

    def _help(self):
        print("""
Usage: pl <command>

Queue:
    pl q start
    pl q stop
    pl q drain
    pl q status
    pl q list
    pl q cancel <id>
    pl q reprio <id> <prio>

Requests:
    pl r list
    pl r status <id>
    pl r info <id>
    pl r logs <id>
    pl r note <id> "text"
    pl r tag <id> tag1 tag2
    pl r pin <id>
    pl r publish <id>
    pl r archive <id>
    pl r trash <id>
""")


def main():
    PipelineCLI(sys.argv).dispatch()


if __name__ == "__main__":
    main()
