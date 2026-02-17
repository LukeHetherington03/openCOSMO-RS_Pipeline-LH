#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import subprocess
import signal
import json

from modules.execution.queue import QueueManager

CONFIG_PATH = os.path.abspath("config/paths.json")

def load_base_dir():
    with open(CONFIG_PATH) as f:
        cfg = json.load(f)
    return cfg["base_dir"]


class QueueCommands:

    @staticmethod
    def dispatch(args):
        if not args:
            print("Usage: pl q <start|stop|status|list|cancel|reprio>")
            return

        sub = args[0]

        if sub == "start": return QueueCommands.start()
        if sub == "stop": return QueueCommands.stop()
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
        base_dir = load_base_dir()
        QueueManager.init_base_dir(base_dir)

        if QueueManager.is_worker_running():
            print(f"Worker already running (PID {QueueManager.read_pid()})")
            return

        subprocess.Popen(
            [sys.executable, "-m", "modules.execution.worker", base_dir],
            stdout=open(os.devnull, "w"),
            stderr=open(os.devnull, "w")
        )

        print("Queue worker started.")

    @staticmethod
    def stop():
        base_dir = load_base_dir()
        QueueManager.init_base_dir(base_dir)

        pid = QueueManager.read_pid()
        if pid:
            try:
                os.kill(pid, signal.SIGTERM)
                print(f"Stopping worker (PID {pid})...")
            except ProcessLookupError:
                print("Worker not running.")
            QueueManager.clear_pid()
            print("Worker stopped successfully.")
        else:
            print("No worker is currently running.")

    # ------------------------------------------------------------
    # Queue inspection
    # ------------------------------------------------------------
    @staticmethod
    def status():
        base_dir = load_base_dir()
        QueueManager.init_base_dir(base_dir)

        pending = len(os.listdir(QueueManager.DIRS["pending"]))
        running = len(os.listdir(QueueManager.DIRS["running"]))
        completed = len(os.listdir(QueueManager.DIRS["completed"]))
        failed = len(os.listdir(QueueManager.DIRS["failed"]))

        print(f"Worker running: {QueueManager.is_worker_running()}")
        print(f"Pending:   {pending}")
        print(f"Running:   {running}")
        print(f"Completed: {completed}")
        print(f"Failed:    {failed}")

    @staticmethod
    def list():
        base_dir = load_base_dir()
        QueueManager.init_base_dir(base_dir)

        print("Pending:")
        for f in os.listdir(QueueManager.DIRS["pending"]):
            print("  -", f[:-5])

        print("\nRunning:")
        for f in os.listdir(QueueManager.DIRS["running"]):
            print("  -", f[:-5])

    # ------------------------------------------------------------
    # Queue operations
    # ------------------------------------------------------------
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
