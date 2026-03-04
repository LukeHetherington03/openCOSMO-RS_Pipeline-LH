import os
from datetime import datetime


class LogHelper:
    """
    Unified logging helper for Request and Job logs.

    =========================================================================
    LOG LEVEL TAXONOMY
    =========================================================================

    [INFO]       General stage progress — loading files, counts, config values
    [CONFIG]     Values resolved from config/defaults/parameters at startup
    [RESOURCES]  Allocation decisions — cores, workers, timing
    [START]      Item beginning — emitted when worker picks up an item
    [COMPLETE]   Item finished successfully — with → detail block beneath
    [FAILED]     Item failed — with → error block beneath
    [RESUME]     Item skipped because checkpoint already exists
    [WARNING]    Non-fatal issue — missing file, fallback triggered, etc.
    [ERROR]      Stage-level error before fail/raise
    [POOL]       Pool lifecycle — started, drained
    [STAGE]      Stage lifecycle — started, complete, wall time
    [SKIP]       Item deliberately skipped (missing input, invalid data, etc.)
    [DEBUG]      Verbose detail — written to debug log only, not main log

    =========================================================================
    FORMAT
    =========================================================================

    Standard line:
        [2025-01-15T14:23:01] [LEVEL]     message

    Continuation (detail rows beneath COMPLETE / FAILED):
        [2025-01-15T14:23:01] [COMPLETE]  item_id  worker=X  pid=Y  elapsed=Z
                                          → key:    value
                                          → key:    value

    Level tags are padded to a fixed width so message columns align across
    all levels. Continuation → lines align to the same message column.

    =========================================================================
    USAGE
    =========================================================================

    # Standard levels
    LogHelper.info(path, "Loaded 36 entries")
    LogHelper.config(path, "enable_fallback=True  strict=False")
    LogHelper.resources(path, "cores_per_item=3  n_workers=10")
    LogHelper.warning(path, "Missing XYZ for ABCDEF_c003 — skipping")
    LogHelper.error(path, "Stage input not found")

    # Pool lifecycle
    LogHelper.pool_started(path, workers=10, items=36, cores_per_item=3)
    LogHelper.pool_drained(path, successful=34, failed=2, wall=182.4)

    # Stage lifecycle
    LogHelper.stage_started(path, "OrcacosmoStage")
    LogHelper.stage_complete(path, "OrcacosmoStage", items=36, ok=34,
                             failed=2, wall=182.4)

    # Item events
    LogHelper.item_start(path, "ABCDEF_c001",
                         worker="ForkPoolWorker-1", pid=84318)

    LogHelper.item_complete(path, "ABCDEF_c001",
                            worker="ForkPoolWorker-1", pid=84318,
                            elapsed=38.2,
                            details=[
                                ("method",   "TZVPD"),
                                ("fallback", "no"),
                                ("output",   "orcacosmo_outputs/ABCDEF_c001.orcacosmo"),
                            ])

    LogHelper.item_failed(path, "ABCDEF_c002",
                          worker="ForkPoolWorker-2", pid=84319,
                          elapsed=41.1,
                          error="TZVPD and TZVP both failed — cpcm missing")

    LogHelper.item_resume(path, "ABCDEF_c001")

    LogHelper.item_skip(path, "ABCDEF_c003",
                        reason="Missing XYZ file")

    # Legacy compat (still works — delegates to info())
    LogHelper.write(path, "some message")
    LogHelper.header(path, "Section title")
    LogHelper.section(path, "Sub-section")
    """

    # ── Layout constants ──────────────────────────────────────────────────────

    # Width of the bracketed level tag including trailing space,
    # e.g. "[COMPLETE]  " = 12 chars.  All tags are padded to this.
    _LEVEL_WIDTH = 12

    # Timestamp width: "[2025-01-15T14:23:01] " = 22 chars
    _TS_WIDTH = 22

    # Full prefix width used to align continuation → lines
    _INDENT = _TS_WIDTH + _LEVEL_WIDTH   # 34 chars

    # ── Internals ─────────────────────────────────────────────────────────────

    @staticmethod
    def _timestamp() -> str:
        return datetime.now().isoformat(timespec="seconds")

    @staticmethod
    def _tag(level: str) -> str:
        """Return a fixed-width bracketed tag, e.g. '[COMPLETE]  '."""
        bracket = f"[{level}]"
        return bracket.ljust(LogHelper._LEVEL_WIDTH)

    @staticmethod
    def _write_line(path: str, level: str, message: str, echo: bool = False):
        """Write a single timestamped log line."""
        os.makedirs(os.path.dirname(path), exist_ok=True)
        line = f"[{LogHelper._timestamp()}] {LogHelper._tag(level)}{message}"
        with open(path, "a") as f:
            f.write(line + "\n")
        if echo:
            print(line)

    @staticmethod
    def _write_continuation(path: str, key: str, value: str, echo: bool = False):
        """
        Write a continuation detail line aligned under the message column.

            <indent spaces>→ key:    value

        Key is padded to 10 chars so values align across rows.
        """
        os.makedirs(os.path.dirname(path), exist_ok=True)
        key_col = f"{key}:".ljust(10)
        line = f"{' ' * LogHelper._INDENT}→ {key_col} {value}"
        with open(path, "a") as f:
            f.write(line + "\n")
        if echo:
            print(line)

    @staticmethod
    def _append_raw(path: str, text: str):
        """Append a raw pre-formatted line (no timestamp, no level tag)."""
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "a") as f:
            f.write(text + "\n")

    # ── Standard levels ───────────────────────────────────────────────────────

    @staticmethod
    def info(path: str, message: str, echo: bool = False):
        LogHelper._write_line(path, "INFO", message, echo)

    @staticmethod
    def config(path: str, message: str, echo: bool = False):
        LogHelper._write_line(path, "CONFIG", message, echo)

    @staticmethod
    def resources(path: str, message: str, echo: bool = False):
        LogHelper._write_line(path, "RESOURCES", message, echo)

    @staticmethod
    def warning(path: str, message: str, echo: bool = False):
        LogHelper._write_line(path, "WARNING", message, echo)

    @staticmethod
    def error(path: str, message: str, echo: bool = False):
        LogHelper._write_line(path, "ERROR", message, echo)

    @staticmethod
    def debug(path: str, message: str, echo: bool = False):
        """Write to a separate debug log only — not the main stage log."""
        LogHelper._write_line(path, "DEBUG", message, echo)

    # ── Pool lifecycle ────────────────────────────────────────────────────────

    @staticmethod
    def pool_started(
        path: str,
        workers: int,
        items: int,
        cores_per_item: int,
        echo: bool = False,
    ):
        msg = (
            f"started  "
            f"workers={workers}  "
            f"items={items}  "
            f"cores_per_item={cores_per_item}"
        )
        LogHelper._write_line(path, "POOL", msg, echo)

    @staticmethod
    def pool_drained(
        path: str,
        successful: int,
        failed: int,
        wall: float,
        echo: bool = False,
    ):
        msg = (
            f"drained  "
            f"successful={successful}  "
            f"failed={failed}  "
            f"wall={wall:.1f}s"
        )
        LogHelper._write_line(path, "POOL", msg, echo)

    # ── Stage lifecycle ───────────────────────────────────────────────────────

    @staticmethod
    def stage_started(path: str, stage_name: str, echo: bool = False):
        LogHelper._write_line(path, "STAGE", f"{stage_name} started", echo)

    @staticmethod
    def stage_complete(
        path: str,
        stage_name: str,
        items: int,
        ok: int,
        failed: int,
        wall: float,
        echo: bool = False,
    ):
        msg = (
            f"{stage_name} complete  "
            f"items={items}  "
            f"ok={ok}  "
            f"failed={failed}  "
            f"wall={wall:.1f}s"
        )
        LogHelper._write_line(path, "STAGE", msg, echo)

    @staticmethod
    def stage_failed(
        path: str,
        stage_name: str,
        reason: str,
        echo: bool = False,
    ):
        LogHelper._write_line(
            path, "STAGE", f"{stage_name} FAILED  reason={reason}", echo
        )

    # ── Item events ───────────────────────────────────────────────────────────

    @staticmethod
    def item_start(
        path: str,
        item_id: str,
        worker: str,
        pid: int,
        echo: bool = False,
    ):
        msg = f"{item_id}  worker={worker}  pid={pid}"
        LogHelper._write_line(path, "START", msg, echo)

    @staticmethod
    def item_complete(
        path: str,
        item_id: str,
        worker: str,
        pid: int,
        elapsed: float,
        details: list,          # list of (key, value) tuples — stage-specific
        echo: bool = False,
    ):
        """
        Write a COMPLETE line followed by indented → detail rows.

        details = [
            ("method",   "TZVPD"),
            ("fallback", "no"),
            ("output",   "orcacosmo_outputs/ABCDEF_c001.orcacosmo"),
        ]
        """
        msg = f"{item_id}  worker={worker}  pid={pid}  elapsed={elapsed:.1f}s"
        LogHelper._write_line(path, "COMPLETE", msg, echo)
        for key, value in details:
            LogHelper._write_continuation(path, key, str(value), echo)

    @staticmethod
    def item_failed(
        path: str,
        item_id: str,
        worker: str,
        pid: int,
        elapsed: float,
        error: str,
        echo: bool = False,
    ):
        """Write a FAILED line followed by an indented → error row."""
        msg = f"{item_id}  worker={worker}  pid={pid}  elapsed={elapsed:.1f}s"
        LogHelper._write_line(path, "FAILED", msg, echo)
        LogHelper._write_continuation(path, "error", error, echo)

    @staticmethod
    def item_resume(path: str, item_id: str, echo: bool = False):
        """Item skipped because a valid checkpoint already exists."""
        LogHelper._write_line(
            path, "RESUME", f"{item_id}  checkpoint exists — skipping", echo
        )

    @staticmethod
    def item_skip(
        path: str,
        item_id: str,
        reason: str,
        echo: bool = False,
    ):
        """Item deliberately skipped due to missing/invalid input."""
        LogHelper._write_line(path, "SKIP", f"{item_id}  reason={reason}", echo)

    # ── Dividers (for human readability in long logs) ─────────────────────────

    @staticmethod
    def divider(path: str, echo: bool = False):
        """Write a plain divider line — no timestamp, no level."""
        LogHelper._append_raw(path, "─" * 80)
        if echo:
            print("─" * 80)

    # ── Legacy compatibility ──────────────────────────────────────────────────
    # These delegate to info() / divider() so existing callers keep working
    # without changes. New code should use the typed methods above.

    @staticmethod
    def write(path: str, message: str, indent: int = 0, echo: bool = False):
        """Legacy: timestamped write with optional indent. Delegates to info()."""
        prefix = " " * indent
        LogHelper.info(path, f"{prefix}{message}", echo)

    @staticmethod
    def header(path: str, title: str, echo: bool = False):
        """Legacy: section header. Kept for backward compat."""
        LogHelper.divider(path, echo)
        LogHelper._write_line(path, "STAGE", title, echo)
        LogHelper.divider(path, echo)

    @staticmethod
    def section(path: str, title: str, echo: bool = False):
        """Legacy: sub-section marker. Kept for backward compat."""
        LogHelper._write_line(path, "INFO", f"── {title} ──", echo)