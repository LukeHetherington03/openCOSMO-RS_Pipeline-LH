#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import signal

from modules.execution.queue import QueueManager


class QueueCommands:

    @staticmethod
    def dispatch(args):
        if not args:
            print("Usage: pl q <start|stop|drain|status|list|cancel|reprio>")
            return

        sub = args[0]

        if sub == "start": return QueueCommands.start()
        if sub == "stop": return QueueCommands.stop()
        if sub == "drain": return QueueCommands.drain()
        if sub == "status": return QueueCommands.status()
        if sub == "list": return QueueCommands.list()
        if sub == "cancel": return QueueCommands.cancel(args[1:])
        if sub == "reprio": return QueueCommands.reprio(args[1:])

        print(f"Unknown queue command: {sub}")

    # ------------------------------------------------------------
    # Worker control
    # ------------------------------------------------------------
    @staticmethod
    def start():
        if QueueManager.is_worker_running():
            print(f"Worker already running (PID {QueueManager.read_pid()})")
            return

        subprocess.Popen(
            [sys.executable, "-m", "modules.execution.worker"],
            stdout=open(os.devnull, "w"),
            stderr=open(os.devnull, "w")
        )

        print("Queue worker started.")

    @staticmethod
    def stop():
        pid = QueueManager.read_pid()
        if pid:
            try:
                os.kill(pid, signal.SIGTERM)
                print(f"Stopped worker (PID {pid})")
            except ProcessLookupError:
                print("Worker not running.")
        QueueManager.clear_pid()

    @staticmethod
    def drain():
        QueueManager.enable_drain()
        print("Drain mode enabled.")

    # ------------------------------------------------------------
    # Queue inspection
    # ------------------------------------------------------------
    @staticmethod
    def status():
        stats = QueueManager.stats()
        print(f"Worker running: {QueueManager.is_worker_running()}")
        print(f"Pending:   {stats['pending']}")
        print(f"Running:   {stats['running']}")
        print(f"Completed: {stats['completed']}")
        print(f"Failed:    {stats['failed']}")

    @staticmethod
    def list():
        print("Pending:")
        for rid in QueueManager.list_pending():
            print(f"  - {rid}")

        print("\nRunning:")
        for rid in QueueManager.list_running():
            print(f"  - {rid}")

    # ------------------------------------------------------------
    # Queue operations
    # ------------------------------------------------------------
    @staticmethod
    def cancel(args):
        rid = args[0]
        QueueManager.cancel(rid)
        print(f"Cancelled {rid}")

    @staticmethod
    def reprio(args):
        rid, prio = args[0], int(args[1])
        QueueManager.reprioritise(rid, prio)
        print(f"Updated priority for {rid} → {prio}")
