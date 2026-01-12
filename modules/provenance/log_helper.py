import os
from datetime import datetime

class LogHelper:
    """
    Unified logging helper for Request and Job logs.
    Provides timestamped logging with optional indentation.
    """

    @staticmethod
    def timestamp():
        return datetime.now().isoformat(timespec="seconds")

    @staticmethod
    def write(path, message, indent=0, echo=False):
        os.makedirs(os.path.dirname(path), exist_ok=True)

        prefix = " " * indent
        line = f"[{LogHelper.timestamp()}] {prefix}{message}"

        with open(path, "a") as f:
            f.write(line + "\n")

        if echo:
            print(line)

    @staticmethod
    def header(path, title, echo=False):
        line = f"===== {title} ====="
        LogHelper.write(path, line, indent=0, echo=echo)

    @staticmethod
    def section(path, title, echo=False):
        line = f"--- {title} ---"
        LogHelper.write(path, line, indent=0, echo=echo)
