#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os

from modules.request_tools.view import RequestView
from modules.request_tools.annotation import RequestAnnotation


BASE_DIR = os.path.abspath("pipeline_data")


class RequestCommands:

    # ------------------------------------------------------------
    # Dispatcher
    # ------------------------------------------------------------
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
        for rid in RequestView.list_all(BASE_DIR):
            rv = RequestView.load(BASE_DIR, rid)
            ps = rv.pipeline_state
            print(f"{rid} | {ps['state']} | {ps['stage']} | {rv.data.get('state', {}).get('user_status', '')}")

    @staticmethod
    def status(args):
        if not args:
            print("Usage: pl r status <id>")
            return
        rid = args[0]
        rv = RequestView.load(BASE_DIR, rid)
        print(rv.pretty_status())

    @staticmethod
    def info(args):
        if not args:
            print("Usage: pl r info <id>")
            return
        rid = args[0]
        rv = RequestView.load(BASE_DIR, rid)
        print(rv.pretty_info())

    @staticmethod
    def logs(args):
        if not args:
            print("Usage: pl r logs <id>")
            return
        rid = args[0]
        rv = RequestView.load(BASE_DIR, rid)
        print(rv.latest_logs())

    # ------------------------------------------------------------
    # Annotation operations
    # ------------------------------------------------------------
    @staticmethod
    def note(args):
        if len(args) < 2:
            print("Usage: pl r note <id> \"text\"")
            return
        rid = args[0]
        text = " ".join(args[1:])
        RequestAnnotation(BASE_DIR).note(rid, text)
        print(f"Added note to {rid}")

    @staticmethod
    def tag(args):
        if len(args) < 2:
            print("Usage: pl r tag <id> tag1 tag2 ...")
            return
        rid = args[0]
        tags = args[1:]
        RequestAnnotation(BASE_DIR).tag(rid, tags)
        print(f"Added tags to {rid}")

    @staticmethod
    def pin(args):
        rid = args[0]
        RequestAnnotation(BASE_DIR).pin(rid)
        print(f"Pinned {rid}")

    @staticmethod
    def publish(args):
        rid = args[0]
        RequestAnnotation(BASE_DIR).publish(rid)
        print(f"Published {rid}")

    @staticmethod
    def archive(args):
        rid = args[0]
        RequestAnnotation(BASE_DIR).archive(rid)
        print(f"Archived {rid}")

    @staticmethod
    def trash(args):
        rid = args[0]
        RequestAnnotation(BASE_DIR).trash(rid)
        print(f"Trashed {rid}")
