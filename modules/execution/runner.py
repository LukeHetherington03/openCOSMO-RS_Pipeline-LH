#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modules/execution/runner.py

Minimal, pure executor for a Request.

Responsibilities:
  - load the Request
  - run each Job sequentially
  - allocate resources before each job and inject into parameters
  - write logs to the request directory
  - update pipeline_state.json via Request/Job APIs

No:
  - nohup / scheduler lock / pause/resume/stop / control.json / CLI entrypoint
"""

import json
import os
import shutil
import time
from datetime import datetime

from modules.build.request_manager import Request
from modules.build.job_manager import Job
from modules.execution.resource_allocator import ResourceAllocator

_VERBOSE: bool = False   # set True by start_worker(-v) or direct mode in main.py


def _progress(msg: str) -> None:
    """Print a timestamped progress line to stdout when verbose mode is on."""
    if _VERBOSE:
        ts = datetime.now().strftime("%H:%M:%S")
        print(f"[{ts}] {msg}", flush=True)


# Extra output files to copy to final_outputs/ per last stage (beyond canonical JSON)
_STAGE_EXTRA_FILES = {
    "solubility":   ["solubility_results.csv", "solubility_human_summary.txt"],
    "orcacosmo":    ["orcacosmo_summary.csv"],
    "optimisation": ["optimisation_summary.csv"],
    "pruning":      ["pruning_summary.csv"],
}


class PipelineRunner:

    @classmethod
    def run_request(cls, request_id: str, base_dir: str):
        """Entry point used exclusively by QueueWorker."""
        req = Request.load(base_dir, request_id)
        _progress(f"Request {request_id}  starting  ({len(req.jobs)} jobs)")
        req.log_header("Pipeline started")
        cls._log_pipeline_header(req)
        cls._execute_pipeline(req)
        _progress(f"Request {request_id}  complete")

    # =========================================================================
    # Pipeline header — echoed once at the very top of the request log
    # =========================================================================

    @classmethod
    def _log_pipeline_header(cls, req: "Request"):
        """Echo title, pipeline version, and planned stage sequence."""
        config  = req.parameters.get("config", {})
        version = config.get("pipeline_version", "unknown")
        title   = req.parameters.get("title", req.request_id)

        stages = []
        for job_id in req.jobs:
            try:
                stages.append(Job.load(req, job_id).stage)
            except Exception:
                stages.append("?")

        req.log(f"Title          : {title}")
        req.log(f"Pipeline ver.  : {version}")
        req.log(f"Stage sequence : {' -> '.join(stages)}")

    # =========================================================================
    # Main execution loop
    # =========================================================================

    @classmethod
    def _execute_pipeline(cls, req: "Request"):
        """Sequential job execution with per-job resource allocation."""
        if not req.jobs:
            raise RuntimeError(f"Request {req.request_id} has no jobs registered.")

        config     = req.parameters.get("config", {})
        parameters = req.parameters
        allocator  = ResourceAllocator(config=config, parameters=parameters)

        t_pipeline_start = time.monotonic()
        completed_jobs   = []   # (job_id, stage) in execution order
        n_jobs           = len(req.jobs)

        start_index   = cls._find_resume_index(req)
        current_index = start_index
        current_job   = Job.load(req, req.jobs[current_index])

        while True:
            # -- Resource allocation ------------------------------------------
            n_items = cls._precount_items(current_job) or 1
            result  = allocator.allocate(stage=current_job.stage, n_items=n_items)

            req.log(f"[RESOURCES] stage={current_job.stage}  n_items={n_items}")
            for line in result.detail_lines():
                req.log(line)
            if result.warning:
                req.log(f"[WARNING] {result.warning}")

            current_job.parameters.update(result.as_params())

            # -- Execute -------------------------------------------------------
            _progress(
                f"[{current_index + 1}/{n_jobs}] {current_job.stage}  starting"
                f"  ({n_items} items · {result.n_workers} workers"
                f" × {result.cores_per_item} cores)"
            )
            req.log(f"Running job: {current_job.job_id}")
            current_job.run()
            req.log(f"Completed job: {current_job.job_id}")
            try:
                with open(current_job.job_state_path) as _f:
                    _st = json.load(_f)
                _wall = cls._fmt_duration(cls._job_wall_seconds(_st))
                _ok   = _st.get("n_ok", "?")
                _fail = _st.get("n_failed", 0)
                _fstr = f"  {_fail} failed" if _fail else ""
                _progress(
                    f"[{current_index + 1}/{n_jobs}] {current_job.stage}  done"
                    f"  ({_ok} ok{_fstr} · {_wall})"
                )
            except Exception:
                _progress(f"[{current_index + 1}/{n_jobs}] {current_job.stage}  done")

            completed_jobs.append((current_job.job_id, current_job.stage))

            next_job = req.create_next_job_by_index(current_index)
            if next_job is None:
                t_total = time.monotonic() - t_pipeline_start
                cls._log_completion_summary(req, completed_jobs, t_total)
                cls._log_pipeline_warnings_summary(req, completed_jobs)
                cls._write_final_outputs(req, completed_jobs)
                req.log_header("Pipeline complete")
                return

            current_index += 1
            current_job    = next_job

    # =========================================================================
    # Completion summary
    # =========================================================================

    @classmethod
    def _log_completion_summary(cls, req: "Request", completed_jobs: list, t_total: float):
        """
        Emit a stage-by-stage summary table and failure detail block
        immediately before the final Pipeline complete banner.

        Reads job_state.json for wall time and item counts. Falls back
        gracefully if any state file is missing or malformed.

        Example output:

          PIPELINE SUMMARY
          Stage              Items    OK  Failed      Wall  Items/hr
          ------------------------------------------------------------
          cleaning               1     1       -      0.8s      4.5k
          generation             2     2       -     47.1s       153
          pruning              400   400       -      0.9s     1600k
          optimisation           4     4       -    8m 12s        29
          orcacosmo              4     4       -   1h 14m          3
          solubility             2     2       -     12.3s       586
          ------------------------------------------------------------
          Total                                    1h 23m

        Followed by a FAILED ITEMS block if any items failed.
        """
        SEP = "-" * 62

        rows = []
        any_failures = False

        for job_id, stage in completed_jobs:
            try:
                job = Job.load(req, job_id)
                with open(job.job_state_path) as f:
                    state = json.load(f)
                n_items  = state.get("n_items",  0)
                n_ok     = state.get("n_ok",     0)
                n_failed = state.get("n_failed", 0)
                wall_s   = cls._job_wall_seconds(state)
            except Exception:
                n_items = n_ok = n_failed = 0
                wall_s  = None

            if n_failed:
                any_failures = True

            rows.append({
                "stage":    stage,
                "job_id":   job_id,
                "n_items":  n_items,
                "n_ok":     n_ok,
                "n_failed": n_failed,
                "wall_s":   wall_s,
            })

        hdr = (
            f"  {'Stage':<18} {'Items':>6}  {'OK':>4}  {'Failed':>6}"
            f"  {'Wall':>8}  {'Items/hr':>9}"
        )

        req.log("")
        req.log("PIPELINE SUMMARY")
        req.log(SEP)
        req.log(hdr)
        req.log(SEP)

        for r in rows:
            failed_str = str(r["n_failed"]) if r["n_failed"] else "-"
            req.log(
                f"  {r['stage']:<18} {r['n_items']:>6}  {r['n_ok']:>4}"
                f"  {failed_str:>6}  {cls._fmt_duration(r['wall_s']):>8}"
                f"  {cls._fmt_rate(r['n_items'], r['wall_s']):>9}"
            )

        req.log(SEP)
        req.log(f"  {'Total':<18} {'':>6}  {'':>4}  {'':>6}  {cls._fmt_duration(t_total):>8}")
        req.log("")

        # -- Failure detail block (only present when something failed) ---------
        if any_failures:
            req.log("FAILED ITEMS")
            req.log(SEP)
            for r in rows:
                if not r["n_failed"]:
                    continue
                req.log(f"  {r['stage']} ({r['n_failed']} failed):")
                try:
                    job = Job.load(req, r["job_id"])
                    for item_id in job.failed_items:
                        reason = job.failure_reasons.get(item_id, "no reason recorded")
                        req.log(f"    x  {item_id}  --  {reason}")
                except Exception:
                    req.log("    (could not load failure details)")
            req.log("")


    # =========================================================================
    # Pipeline warnings summary
    # =========================================================================

    @classmethod
    def _log_pipeline_warnings_summary(cls, req: "Request", completed_jobs: list):
        """
        Aggregate [WARNING] lines from all stage logs AND pipeline-level
        warnings from request.log into a single PIPELINE WARNINGS block,
        written to request.log and saved as final_outputs/pipeline_warnings.txt.
        """
        SEP = "-" * 62
        all_sections = {}   # section_name -> [message_str, ...]

        # ── Pipeline-level warnings (resource allocator, etc.) ────────────────
        pipeline_warnings = []
        try:
            with open(req.request_log_path) as f:
                for line in f:
                    if "] [INFO] " in line:
                        parts = line.split("] [INFO] ", 1)
                        msg = parts[1].strip() if len(parts) == 2 else ""
                        if msg.startswith("[WARNING] "):
                            pipeline_warnings.append(msg[len("[WARNING] "):])
        except Exception:
            pass
        if pipeline_warnings:
            all_sections["pipeline"] = pipeline_warnings

        # ── Stage-level warnings ───────────────────────────────────────────────
        for job_id, stage in completed_jobs:
            job      = Job.load(req, job_id)
            log_path = os.path.join(job.job_dir, "stage_logs", "stage.log")
            if not os.path.exists(log_path):
                continue
            warnings = []
            try:
                with open(log_path) as f:
                    for line in f:
                        if "] [WARNING] " in line:
                            parts = line.split("] [WARNING] ", 1)
                            msg = parts[1].rstrip() if len(parts) == 2 else line.rstrip()
                            warnings.append(msg)
            except Exception:
                pass
            if warnings:
                all_sections[stage] = warnings

        # Deduplicate within each section (identical messages from repeated allocations etc.)
        for section in all_sections:
            seen = []
            for msg in all_sections[section]:
                if msg not in seen:
                    seen.append(msg)
            all_sections[section] = seen

        total = sum(len(v) for v in all_sections.values())

        # ── Build block ────────────────────────────────────────────────────────
        out_lines = ["", "PIPELINE WARNINGS", SEP]
        if not all_sections:
            out_lines.append("  No warnings — pipeline completed cleanly.")
        else:
            out_lines.append(f"  {total} warning(s) across {len(all_sections)} section(s):")
            out_lines.append("")
            for section, warnings in all_sections.items():
                out_lines.append(f"  [{section.upper()}]  ({len(warnings)} warning(s))")
                for w in warnings:
                    out_lines.append(f"    ! {w}")
                out_lines.append("")
        out_lines += [SEP, ""]

        # ── Write to request.log ───────────────────────────────────────────────
        for line in out_lines:
            req.log(line)

        # ── Save to final_outputs/pipeline_warnings.txt ───────────────────────
        final_dir = os.path.join(req.request_dir, "final_outputs")
        os.makedirs(final_dir, exist_ok=True)
        try:
            with open(os.path.join(final_dir, "pipeline_warnings.txt"), "w") as f:
                f.write("\n".join(out_lines) + "\n")
        except Exception:
            pass

    # =========================================================================
    # Final outputs directory
    # =========================================================================

    @classmethod
    def _write_final_outputs(cls, req: "Request", completed_jobs: list):
        """
        Collect key outputs into {request_dir}/final_outputs/.

        Layout:
          molecule_metadata/       — per-molecule JSON from cleaning stage
          {last_stage_output.*}    — canonical + extra files from last stage
          pipeline_warnings.txt    — already written by _log_pipeline_warnings_summary
        """
        final_dir = os.path.join(req.request_dir, "final_outputs")
        os.makedirs(final_dir, exist_ok=True)

        cleaning_job = None
        last_job_id, last_stage = None, None
        for job_id, stage in completed_jobs:
            if stage == "cleaning":
                cleaning_job = Job.load(req, job_id)
            last_job_id, last_stage = job_id, stage

        # ── Last stage outputs ────────────────────────────────────────────────
        if last_job_id and last_stage != "cleaning":
            try:
                last_job  = Job.load(req, last_job_id)
                stage_out = last_job.get_output_file()
                if stage_out and os.path.exists(stage_out):
                    shutil.copy(stage_out,
                                os.path.join(final_dir, os.path.basename(stage_out)))
                for fname in _STAGE_EXTRA_FILES.get(last_stage, []):
                    p = os.path.join(last_job.outputs_dir, fname)
                    if os.path.exists(p):
                        shutil.copy(p, os.path.join(final_dir, fname))
            except Exception:
                pass

        # ── molecule_metadata/ from cleaning stage ────────────────────────────
        if cleaning_job:
            meta_dir = os.path.join(cleaning_job.outputs_dir, "molecule_metadata")
            if os.path.isdir(meta_dir):
                dest = os.path.join(final_dir, "molecule_metadata")
                if os.path.exists(dest):
                    shutil.rmtree(dest)
                shutil.copytree(meta_dir, dest)

        req.log(f"[OUTPUT] final_outputs/ → {final_dir}")

    # =========================================================================
    # Private helpers
    # =========================================================================

    @classmethod
    def _job_wall_seconds(cls, state: dict):
        """Return wall-clock seconds from a job_state dict, or None."""
        from datetime import datetime, timezone
        start = state.get("start_timestamp")
        stop  = state.get("stop_timestamp")
        if not start or not stop:
            return None
        try:
            fmt = "%Y-%m-%dT%H:%M:%S.%fZ"
            t0  = datetime.strptime(start, fmt).replace(tzinfo=timezone.utc)
            t1  = datetime.strptime(stop,  fmt).replace(tzinfo=timezone.utc)
            return (t1 - t0).total_seconds()
        except Exception:
            return None

    @classmethod
    def _fmt_duration(cls, seconds) -> str:
        """0.8s  /  47.1s  /  8m 12s  /  1h 14m"""
        if seconds is None:
            return "?"
        s = int(seconds)
        if s < 60:
            return f"{seconds:.1f}s"
        m, s = divmod(s, 60)
        if m < 60:
            return f"{m}m {s:02d}s"
        h, m = divmod(m, 60)
        return f"{h}h {m:02d}m"

    @classmethod
    def _fmt_rate(cls, n_items: int, wall_s) -> str:
        """Items/hr as compact string: 29 / 1.4k / 2.1M"""
        if not wall_s or wall_s <= 0 or not n_items:
            return "?"
        rate = n_items / (wall_s / 3600)
        if rate >= 1_000_000:
            return f"{rate / 1_000_000:.1f}M"
        if rate >= 1_000:
            return f"{rate / 1_000:.1f}k"
        return f"{int(rate)}"

    @classmethod
    def _precount_items(cls, job) -> int:
        """
        Best-effort item count before the stage runs.
        Uses job.items if populated (resume), else reads stage_input.
        Returns 0 if count cannot be determined.
        """
        if job.items:
            return len(job.items)

        stage_input = job.parameters.get("stage_input")
        if not stage_input or not os.path.exists(stage_input):
            return 0

        try:
            with open(stage_input) as f:
                data = json.load(f)
            if isinstance(data, (list, dict)):
                return len(data)
        except Exception:
            pass

        try:
            with open(stage_input) as f:
                lines = [l for l in f if l.strip() and not l.startswith("#")]
            return max(0, len(lines) - 1)
        except Exception:
            pass

        return 0

    @classmethod
    def _find_resume_index(cls, req: "Request") -> int:
        """
        Determine which job index to start from.
          - No jobs completed   -> start at 0
          - Some completed      -> resume from next unfinished job
        """
        for i, job_id in enumerate(req.jobs):
            job        = Job.load(req, job_id)
            state_path = job.job_state_path

            if not os.path.exists(state_path):
                return i

            with open(state_path) as f:
                state = json.load(f)

            if state.get("status") != "completed":
                return i

        return len(req.jobs) - 1