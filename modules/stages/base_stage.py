#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modules/stages/base_stage.py

Base class for all pipeline stages.

=========================================================================
RESPONSIBILITIES
=========================================================================

  Lifecycle
    run()           Owns mark_running / mark_complete / mark_failed.
                    Stages implement execute() only — never call
                    mark_complete() themselves.

  Resource attributes (injected by PipelineRunner before job.run())
    self.n_workers        pool size
    self.cores_per_item   per-item core budget (passed to tools)
    self.budget_cores     total core budget for this job

  Checkpoint infrastructure
    self.checkpoints_dir
    _write_checkpoint(item_id, result_dict)
    _load_checkpoint(item_id)  -> dict | None
    _load_all_checkpoints()    -> list[dict]
    _clear_checkpoints()

  Parallel execution
    run_parallel(worker_fn, build_args_fn)
      - worker_fn        module-level pure function, picklable
                         must return a dict containing:
                           "item_id"     str
                           "status"      "ok" | "failed"
                           "log_details" list of (key, value) tuples
                           "error"       str (required when status="failed")
      - build_args_fn    stage method: item -> dict
                         should call self._base_args(item) and extend it

  Standard args skeleton
    _base_args(item) -> dict
      Always contains:
        item_id, stage, outputs_dir, workdir,
        cores_per_item, log_path, checkpoints_dir

  Logging
    All [STAGE] / [POOL] / [START] / [COMPLETE] / [FAILED] / [RESUME]
    events are emitted by BaseStage automatically.
    Stages emit [CONFIG] / [INFO] / [WARNING] / [SKIP] themselves.

  Debug context
    Written to stage_logs/stage_context.log only — never the main stage log.

=========================================================================
STRICT MODE
=========================================================================

  Strict mode is stage-specific. BaseStage reads self.strict_mode
  as a fallback (default False) but never sets it — each stage's
  _load_stage_config() sets self.strict_mode using its own hierarchy:

    pipeline_spec args > request-level strict > stage defaults

  BaseStage.run_parallel() honours self.strict_mode for pool abort.

=========================================================================
CHECKPOINT CONTRACT
=========================================================================

  Every worker function must return a dict with at minimum:

    {
        "item_id":     str,           # matches the item passed to worker
        "status":      "ok",          # or "failed"
        "log_details": [              # list of (key, value) for [COMPLETE]
            ("key", "value"),
            ...
        ],
        # "error": str                # required when status == "failed"
        # ... stage-specific payload
    }

  BaseStage strips "log_details" before writing the checkpoint to disk.
  The checkpoint on disk contains everything else.

"""

import json
import os
import shutil
import time
from concurrent.futures import (
    CancelledError,
    BrokenExecutor,
    ProcessPoolExecutor,
    as_completed,
)
from datetime import datetime
from pathlib import Path

from modules.utils.atomic_write import AtomicWriter
from modules.utils.log_helper import LogHelper


# ─────────────────────────────────────────────────────────────────────────────
# Module-level executor slot  [new 03/03/2026]
# ─────────────────────────────────────────────────────────────────────────────
# QueueWorker's SIGTERM handler calls _cancel_active_executor().
# run_parallel() registers and clears this slot around the pool context manager
# so the signal handler can always reach the live executor.

_active_executor: ProcessPoolExecutor | None = None


def _cancel_active_executor() -> bool:
    """
    Cancel any live ProcessPoolExecutor.

    Called from QueueWorker._handle_sigterm().  Cancels all pending futures
    and requests the pool to shut down without waiting.  Running futures
    complete their own work (writing checkpoints) before their worker
    processes are terminated via the process group signal.

    Returns True if an executor was active, False if idle.
    Safe to call from a signal handler — no locks, no allocations.
    """
    global _active_executor
    ex = _active_executor
    if ex is not None:
        try:
            ex.shutdown(wait=False, cancel_futures=True)
        except Exception:
            pass   # never raise from a signal handler
        return True
    return False


class BaseStage:

    def __init__(self, job):
        self.job        = job
        self.request    = job.request
        self.request_id = job.request_id
        self.job_id     = job.job_id

        self.inputs_dir  = job.inputs_dir
        self.outputs_dir = job.outputs_dir
        self.parameters  = job.parameters
        self.config      = job.config

        self._stage_output = None

        # ── Resource attributes ───────────────────────────────────────
        # Injected by PipelineRunner.allocate() before job.run().
        # Always present — stages can rely on these unconditionally.
        self.n_workers       = int(self.parameters.get("n_workers",       1))
        self.cores_per_item  = int(self.parameters.get("cores_per_item",  1))
        self.budget_cores    = int(self.parameters.get("budget_cores",    1))

        # ── Strict mode default ───────────────────────────────────────
        # Stages override this in _load_stage_config() using their own
        # hierarchy.  BaseStage never forces a value — this is the
        # last-resort fallback only.
        self.strict_mode = False

        # ── Checkpoint directory ──────────────────────────────────────
        self.checkpoints_dir = os.path.join(self.outputs_dir, "checkpoints")

        # ── Log paths ─────────────────────────────────────────────────
        # Both log files live under stage_logs/ within the job directory.
        #   stage.log         — main stage log (events, items, errors)
        #   stage_context.log — context snapshot written once at startup
        stage_logs_dir   = os.path.join(job.job_dir, "stage_logs")
        os.makedirs(stage_logs_dir, exist_ok=True)
        self._stage_log  = os.path.join(stage_logs_dir, "stage.log")
        self._debug_log  = os.path.join(stage_logs_dir, "stage_context.log")

    # =========================================================================
    # Lifecycle — run() owns everything; stages implement execute() only
    # =========================================================================

    def run(self):
        self._write_debug_context()

        stage_name = self.__class__.__name__
        LogHelper.stage_started(self._stage_log, stage_name)

        self.job.mark_running()
        t_start = time.perf_counter()

        try:
            self.execute()

            wall = time.perf_counter() - t_start

            # Count outcomes from job state
            n_items  = len(self.job.items)
            n_ok     = len(self.job.completed_items)
            n_failed = len(self.job.failed_items)

            LogHelper.stage_complete(
                self._stage_log, stage_name,
                items=n_items, ok=n_ok, failed=n_failed, wall=wall,
            )

            self.job.mark_complete()

            self.request.update_pipeline_state(
                last_completed_stage  = self.job.stage,
                last_completed_job_id = self.job.job_id,
            )

        except Exception as e:
            wall = time.perf_counter() - t_start
            LogHelper.stage_failed(self._stage_log, stage_name, reason=str(e))

            self.job.mark_failed(str(e))

            self.request.update_pipeline_state(
                state          = "failed",
                halt_reason    = str(e),
                current_stage  = self.job.stage,
                current_job_id = self.job.job_id,
            )

            raise

    # ── Stages must override this ─────────────────────────────────────────────

    def execute(self):
        raise NotImplementedError("Stage must implement execute()")

    # =========================================================================
    # Parallel execution
    # =========================================================================

    def run_parallel(self, worker_fn, build_args_fn):
        """
        Generic parallel executor with checkpoint/resume.

        Registers the live ProcessPoolExecutor in the module-level
        _active_executor slot so QueueWorker's SIGTERM handler can cancel
        pending futures without waiting.  Running futures complete their
        own checkpoint writes before their processes are terminated via
        the process group signal.

        On CancelledError / BrokenExecutor the method returns cleanly
        rather than raising — the stage proceeds to _assemble_output with
        whatever checkpoints were written.

        Parameters
        ----------
        worker_fn : callable
            Module-level pure function (picklable).
            Signature: worker_fn(args: dict) -> dict
            Must include "item_id", "status", "log_details" in return dict.
            Must include "error" when status == "failed".

        build_args_fn : callable
            Stage method that builds the args dict for one item.
            Signature: build_args_fn(item: str) -> dict
            Should call self._base_args(item) and extend the result.
        """
        global _active_executor

        os.makedirs(self.checkpoints_dir, exist_ok=True)

        # ── Resume: skip already-checkpointed items ───────────────────
        still_pending = []
        for item in list(self.job.pending_items):
            existing = self._load_checkpoint(item)
            if existing is not None:
                LogHelper.item_resume(self._stage_log, item)
                self.job.update_progress(item, success=True)
            else:
                still_pending.append(item)

        if not still_pending:
            LogHelper.info(
                self._stage_log,
                "All items already checkpointed — nothing to run"
            )
            return

        n_workers = min(self.n_workers, len(still_pending))

        LogHelper.pool_started(
            self._stage_log,
            workers        = n_workers,
            items          = len(still_pending),
            cores_per_item = self.cores_per_item,
        )

        t_pool_start  = time.perf_counter()
        n_ok = n_failed = 0
        interrupted     = False
        n_total         = len(still_pending)

        pool = ProcessPoolExecutor(max_workers=n_workers)
        _active_executor = pool     # register for SIGTERM handler

        try:
            # Submit all pending items
            futures: dict = {}
            for submit_idx, item in enumerate(still_pending, start=1):
                args            = build_args_fn(item)
                future          = pool.submit(worker_fn, args)
                futures[future] = (item, submit_idx, time.perf_counter())
                LogHelper.item_start(
                    self._stage_log,
                    item_id = f"[{submit_idx}/{n_total}] {item}",
                    worker  = "submitted",
                    pid     = 0,          # PID comes back with the result
                )

            # Collect results as they complete
            n_done = 0
            for future in as_completed(futures):
                item, submit_idx, t_submit = futures[future]
                elapsed = time.perf_counter() - t_submit
                n_done += 1
                pos_tag = f"[{n_done}/{n_total}]"

                try:
                    exc = future.exception()
                except (CancelledError, BrokenExecutor):
                    # Pool was shut down externally (pl q stop path)
                    interrupted = True
                    LogHelper.warning(
                        self._stage_log,
                        f"[POOL] {item} — future cancelled (external shutdown)"
                    )
                    continue

                if exc:
                    # ── Worker raised an exception ────────────────────
                    n_failed += 1
                    LogHelper.item_failed(
                        self._stage_log,
                        item_id = f"{pos_tag} {item}",
                        worker  = "unknown",
                        pid     = 0,
                        elapsed = elapsed,
                        error   = str(exc),
                    )
                    self.job.update_progress(item, success=False)

                    if self.strict_mode:
                        pool.shutdown(wait=False, cancel_futures=True)
                        self.fail(
                            f"Strict mode — aborting after failure of {item}"
                        )

                else:
                    # ── Worker returned a result dict ─────────────────
                    result = future.result()

                    status      = result.get("status", "ok")
                    log_details = result.pop("log_details", [])
                    worker_name = result.pop("worker_name", "unknown")
                    worker_pid  = result.pop("worker_pid",  0)

                    if status == "failed":
                        # Worker returned a controlled failure
                        n_failed += 1
                        LogHelper.item_failed(
                            self._stage_log,
                            item_id = f"{pos_tag} {item}",
                            worker  = worker_name,
                            pid     = worker_pid,
                            elapsed = elapsed,
                            error   = result.get("error", "unknown error"),
                        )
                        self.job.update_progress(item, success=False)

                        if self.strict_mode:
                            pool.shutdown(wait=False, cancel_futures=True)
                            self.fail(
                                f"Strict mode — aborting after failure of {item}"
                            )

                    else:
                        # Success — write checkpoint then log
                        n_ok += 1
                        self._write_checkpoint(item, result)

                        LogHelper.item_complete(
                            self._stage_log,
                            item_id = f"{pos_tag} {item}",
                            worker  = worker_name,
                            pid     = worker_pid,
                            elapsed = elapsed,
                            details = log_details,
                        )
                        self.job.update_progress(item, success=True)

        except (CancelledError, BrokenExecutor):
            interrupted = True
            LogHelper.warning(
                self._stage_log, "[POOL] Pool interrupted by external shutdown"
            )

        finally:
            _active_executor = None     # always clear the slot
            try:
                pool.shutdown(wait=False, cancel_futures=True)
            except Exception:
                pass

        wall_pool = time.perf_counter() - t_pool_start

        if interrupted:
            LogHelper.warning(
                self._stage_log,
                f"[POOL] Interrupted — {n_ok} ok, {n_failed} failed before "
                f"shutdown. Checkpoints written for completed items."
            )
        else:
            LogHelper.pool_drained(
                self._stage_log,
                successful = n_ok,
                failed     = n_failed,
                wall       = wall_pool,
            )

    # =========================================================================
    # Standard args skeleton
    # =========================================================================

    def _base_args(self, item: str) -> dict:
        """
        Standard args fields present for every stage's worker.
        Stage _build_args() calls this and extends the result.

        Fields
        ------
        item_id         str   the item identifier
        stage           str   stage name for logging inside worker
        outputs_dir     str   stage outputs directory
        workdir         str   per-item scratch directory
        cores_per_item  int   core budget for this item's tool calls
        log_path        str   stage log path (worker appends via LogHelper)
        checkpoints_dir str   where worker writes its checkpoint
        """
        return {
            "item_id":        item,
            "stage":          self.job.stage,
            "outputs_dir":    self.outputs_dir,
            "workdir":        os.path.join(self.outputs_dir, "workdirs", item),
            "cores_per_item": self.cores_per_item,
            "log_path":       self._stage_log,
            "checkpoints_dir":self.checkpoints_dir,
        }

    # =========================================================================
    # Checkpoint infrastructure
    # =========================================================================

    def _checkpoint_path(self, item_id: str) -> str:
        return os.path.join(self.checkpoints_dir, f"{item_id}.json")

    def _write_checkpoint(self, item_id: str, result: dict):
        """Write result dict atomically. log_details already stripped."""
        os.makedirs(self.checkpoints_dir, exist_ok=True)
        with AtomicWriter(self._checkpoint_path(item_id)) as f:
            json.dump(result, f, indent=2)

    def _load_checkpoint(self, item_id: str) -> dict | None:
        """Return checkpoint dict if it exists and is valid, else None."""
        path = self._checkpoint_path(item_id)
        if not os.path.exists(path):
            return None
        try:
            with open(path) as f:
                data = json.load(f)
            # Treat as missing if it recorded a failure
            if data.get("status") == "failed":
                return None
            return data
        except Exception:
            return None  # corrupt — rerun

    def _load_all_checkpoints(self) -> list[dict]:
        """Load all successful checkpoint dicts, sorted by item_id."""
        results = []
        if not os.path.exists(self.checkpoints_dir):
            return results
        for fname in sorted(os.listdir(self.checkpoints_dir)):
            if not fname.endswith(".json"):
                continue
            path = os.path.join(self.checkpoints_dir, fname)
            try:
                with open(path) as f:
                    data = json.load(f)
                if data.get("status") != "failed":
                    results.append(data)
            except Exception:
                LogHelper.warning(
                    self._stage_log,
                    f"Could not read checkpoint {fname} — skipping"
                )
        return results

    def _clear_checkpoints(self):
        """Remove all checkpoints (fresh run)."""
        if os.path.exists(self.checkpoints_dir):
            shutil.rmtree(self.checkpoints_dir)
        os.makedirs(self.checkpoints_dir, exist_ok=True)

    # =========================================================================
    # Stage input / output helpers
    # =========================================================================

    def get_stage_input(self) -> str:
        stage_input = self.parameters.get("stage_input")
        if not stage_input:
            self.fail("Stage requires 'stage_input' but none was provided.")
        return stage_input

    def set_stage_output(self, filename: str) -> str:
        self._stage_output = self.output_path(filename)
        return self._stage_output

    def get_stage_output(self) -> str:
        return self._stage_output

    # =========================================================================
    # Logging helpers — delegates to LogHelper with stage log path
    # =========================================================================

    def log(self, msg: str, indent: int = 0):
        """General info log."""
        prefix = " " * indent
        LogHelper.info(self._stage_log, f"{prefix}{msg}")

    def log_info(self, msg: str):
        LogHelper.info(self._stage_log, msg)

    def log_config(self, msg: str):
        LogHelper.config(self._stage_log, msg)

    def log_resources(self, msg: str):
        LogHelper.resources(self._stage_log, msg)

    def log_warning(self, msg: str):
        LogHelper.warning(self._stage_log, msg)

    def log_error(self, msg: str):
        LogHelper.error(self._stage_log, msg)

    def log_skip(self, item_id: str, reason: str):
        LogHelper.item_skip(self._stage_log, item_id, reason)

    def log_header(self, title: str):
        LogHelper.header(self._stage_log, title)

    def log_section(self, title: str):
        LogHelper.section(self._stage_log, title)

    def log_debug(self, msg: str):
        """Write to debug log only — never the main stage log."""
        LogHelper.debug(self._debug_log, msg)

    # =========================================================================
    # Item tracking — delegates to job
    # =========================================================================

    def set_items(self, items, sort_by_complexity: bool = False):
        """
        Register items with the job, optionally sorting hardest-first.

        Parameters
        ----------
        items              : list of item IDs
        sort_by_complexity : if True, sort descending by rotatable_bonds
                             from molecule metadata before registering.
                             Reads config.constant_files.metadata_dir.
                             Silently no-ops if metadata_dir is absent or
                             any individual metadata file is unreadable.
        """
        if sort_by_complexity:
            items = self._sort_by_complexity(items)
        self.job.set_items(items)

    def _sort_by_complexity(self, items: list) -> list:
        """
        Sort items descending by rotatable_bonds from molecule metadata.

        Works for both plain inchi_keys (orcacosmo) and inchi_key_confNNN
        lookup_ids (optimisation) -- _confNNN is stripped before lookup.

        Degradation:
          - metadata_dir missing/not a directory -> return items unchanged
          - individual file missing or corrupt   -> score = 0, sorts last
          - any unexpected exception             -> return items unchanged
        """
        metadata_dir = (
            (self.config or {})
            .get("constant_files", {})
            .get("metadata_dir")
        )
        if not metadata_dir or not os.path.isdir(metadata_dir):
            return items

        def _score(item_id: str) -> int:
            inchi_key = item_id.split("_conf")[0]
            meta_path = os.path.join(metadata_dir, f"{inchi_key}.json")
            try:
                with open(meta_path) as f:
                    meta = json.load(f)
                return int(meta.get("rotatable_bonds", 0))
            except Exception:
                return 0

        try:
            return sorted(items, key=_score, reverse=True)
        except Exception:
            return items

    def update_progress(self, item, success: bool = True):
        self.job.update_progress(item, success)

    # =========================================================================
    # Path helpers
    # =========================================================================

    def input_path(self, *parts) -> str:
        return self.job.input_path(*parts)

    def output_path(self, *parts) -> str:
        return self.job.output_path(*parts)

    def require_file(self, path: str, description: str = "required file") -> str:
        if not os.path.exists(path):
            self.fail(f"{description} not found: {path}")
        return path

    # =========================================================================
    # Error helper
    # =========================================================================

    def fail(self, message: str):
        LogHelper.error(self._stage_log, message)
        raise RuntimeError(message)

    # =========================================================================
    # Debug context — written to logs/stage_context.log only
    # =========================================================================

    def _write_debug_context(self):
        """
        Write a structured context snapshot to logs/stage_context.log.

        Layout
        ------
        Header banner  — stage name, timestamp, request/job IDs
        RESOURCES      — available (from parameters.resources) and
                         allocated (injected by PipelineRunner)
        PARAMETERS     — non-path scalar values from self.parameters
        PATHS          — self.parameters["config"] JSON dump (all paths)
        """
        try:
            SEP   = "=" * 68
            stage = self.job.stage.upper()
            ts    = datetime.now().isoformat(timespec="milliseconds") + "Z"
            lines = []

            # ── Header ────────────────────────────────────────────────────
            lines += [
                SEP,
                f"{stage} STAGE — CONTEXT SNAPSHOT",
                f"Timestamp  : {ts}",
                f"Request ID : {self.request_id}",
                f"Job ID     : {self.job_id}",
                SEP,
                "",
            ]

            # ── Resources ─────────────────────────────────────────────────
            resources = self.parameters.get("resources") or {}
            avail_parts = []
            for k, v in resources.items():
                avail_parts.append(f"{k}={v}")
            avail_str = "  ".join(avail_parts) if avail_parts else "not specified"

            lines += [
                "RESOURCES",
                f"  Available  : {avail_str}",
                f"  Allocated  : cores_per_item={self.cores_per_item}"
                f"  n_workers={self.n_workers}"
                f"  budget_cores={self.budget_cores}",
                "",
            ]

            # ── Parameters (non-path scalars) ─────────────────────────────
            # Show everything except "config" (paths) and "resources"
            # (already shown above).  Scalars only — skip nested dicts/lists
            # that aren't useful here.
            SKIP_KEYS = {"config", "resources", "cores_per_item", "n_workers", "budget_cores"}
            params_rows = []
            for k, v in self.parameters.items():
                if k in SKIP_KEYS:
                    continue
                if isinstance(v, (dict, list)):
                    # Show type hint rather than raw dump
                    params_rows.append((k, f"<{type(v).__name__}>"))
                else:
                    params_rows.append((k, str(v)))

            if params_rows:
                col = max(len(k) for k, _ in params_rows) + 2
                lines.append("PARAMETERS")
                for k, v in params_rows:
                    lines.append(f"  {k:<{col}}: {v}")
                lines.append("")

            # ── Paths config JSON ──────────────────────────────────────────
            config = self.parameters.get("config")
            if config:
                lines.append("PATHS CONFIG")
                lines.append(json.dumps(config, indent=2))
                lines.append("")

            with open(self._debug_log, "w") as f:
                f.write("\n".join(lines) + "\n")

        except Exception:
            pass  # never let context logging sink a stage

    # =========================================================================
    # Legacy compat
    # =========================================================================

    def strict(self, section: str) -> bool:
        """Legacy helper — kept for stages not yet migrated."""
        return bool(self.config.get(section, {}).get("strict", False))