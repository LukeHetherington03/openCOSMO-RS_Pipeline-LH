#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import json

from modules.request_tools.view import RequestView
from modules.request_tools.annotation import RequestAnnotation

CONFIG_PATH = os.path.abspath("config/paths.json")

def load_base_dir():
    with open(CONFIG_PATH) as f:
        cfg = json.load(f)
    return cfg["base_dir"]


class RequestCommands:

    @staticmethod
    def dispatch(args):
        if not args:
            print("Usage: pl r <list|status|info|logs|note|tag|pin|publish|archive|trash>")
            return

        sub = args[0]
        rest = args[1:]

        if sub == "list": return RequestCommands.list(rest)
        if sub == "status": return RequestCommands.status(rest)
        if sub == "info": return RequestCommands.info(rest)
        if sub == "logs": return RequestCommands.logs(rest)

        if sub == "note": return RequestCommands.note(rest)
        if sub == "tag": return RequestCommands.tag(rest)
        if sub == "pin": return RequestCommands.pin(rest)
        if sub == "publish": return RequestCommands.publish(rest)
        if sub == "archive": return RequestCommands.archive(rest)
        if sub == "trash": return RequestCommands.trash(rest)

        print(f"Unknown request command: {sub}")

    # ------------------------------------------------------------
    # Read-only operations
    # ------------------------------------------------------------
    @staticmethod
    def list(args):
        base_dir = load_base_dir()
        for rid in RequestView.list_all(base_dir):
            rv = RequestView.load(base_dir, rid)
            ps = rv.pipeline_state
            print(f"{rid} | {ps['state']} | {ps['stage']} | {rv.data.get('state', {}).get('user_status', '')}")

    @staticmethod
    def status(args):
        if not args:
            print("Usage: pl r status <id>")
            return
        base_dir = load_base_dir()
        rid = args[0]
        rv = RequestView.load(base_dir, rid)
        print(rv.pretty_status())

    @staticmethod
    def info(args):
        if not args:
            print("Usage: pl r info <id>")
            return
        base_dir = load_base_dir()
        rid = args[0]
        rv = RequestView.load(base_dir, rid)
        print(rv.pretty_info())

    @staticmethod
    def logs(args):
        if not args:
            print("Usage: pl r logs <id>")
            return
        base_dir = load_base_dir()
        rid = args[0]
        rv = RequestView.load(base_dir, rid)
        print(rv.latest_logs())

    # ------------------------------------------------------------
    # Annotation operations
    # ------------------------------------------------------------
    @staticmethod
    def note(args):
        if len(args) < 2:
            print("Usage: pl r note <id> \"text\"")
            return
        base_dir = load_base_dir()
        rid = args[0]
        text = " ".join(args[1:])
        RequestAnnotation(base_dir).note(rid, text)
        print(f"Added note to {rid}")

    @staticmethod
    def tag(args):
        if len(args) < 2:
            print("Usage: pl r tag <id> tag1 tag2 ...")
            return
        base_dir = load_base_dir()
        rid = args[0]
        tags = args[1:]
        RequestAnnotation(base_dir).tag(rid, tags)
        print(f"Added tags to {rid}")

    @staticmethod
    def pin(args):
        base_dir = load_base_dir()
        rid = args[0]
        RequestAnnotation(base_dir).pin(rid)
        print(f"Pinned {rid}")

    @staticmethod
    def publish(args):
        base_dir = load_base_dir()
        rid = args[0]
        RequestAnnotation(base_dir).publish(rid)
        print(f"Published {rid}")

    @staticmethod
    def archive(args):
        base_dir = load_base_dir()
        rid = args[0]
        RequestAnnotation(base_dir).archive(rid)
        print(f"Archived {rid}")

    @staticmethod
    def trash(args):
        base_dir = load_base_dir()
        rid = args[0]
        RequestAnnotation(base_dir).trash(rid)
        print(f"Trashed {rid}")
