#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modules/execution/resource_allocator.py

Computes core allocation for pipeline stage execution.

=========================================================================
CONFIGURATION  (resource_allocation.json)
=========================================================================

    {
      "total_cores": 40,    physical cores on the node (validation only)
      "max_cores":   32     pipeline budget — what stages are allowed to use
    }

    total_cores  Used only to validate that max_cores does not exceed it.
                 If absent, validation is skipped.

    max_cores    The actual allocation budget.
                 If absent, the emergency fallback is used (see below).

=========================================================================
CORE BUDGET HIERARCHY
=========================================================================

    1. parameters["resources"]["cores"]
       Runtime override supplied when the request is created.
       Must be <= max_cores (validated, not silently clamped).

    2. resource_allocation.json  "max_cores"
       Node-level budget set by the operator.

    3. Emergency fallback  (only when max_cores is absent)
       Samples free cores via psutil, applies headroom rule:
         a. subtract 2  (guaranteed minimum headroom)
         b. round down to nearest valid multiple:
              <= 16 cores -> nearest multiple of 2
              >  16 cores -> nearest multiple of 4

       Examples (subtract 2 first, then round):
         free=40 -> 38 -> mult-of-4 -> 36
         free=41 -> 39 -> mult-of-4 -> 36
         free=42 -> 40 -> mult-of-4 -> 40
         free=18 -> 16 -> mult-of-2 -> 16
         free=15 -> 13 -> mult-of-2 -> 12
         free=4  ->  2 -> mult-of-2 ->  2
         free=2  ->  0 -> floor to  ->  1  (always at least 1)

=========================================================================
FREE CORE CHECK  (advisory only when max_cores IS set)
=========================================================================

    When max_cores is configured, psutil is used only as a warning check.
    If sampled free cores < max_cores the stage logs a warning but proceeds
    with max_cores unchanged — the user declared this budget intentionally.

    When max_cores is NOT configured, psutil determines the budget (case 3
    above) and the headroom rule is applied to the result.

=========================================================================
ALLOCATION LOGIC
=========================================================================

    Given a budget of `effective_cores` and `n_items`:

        cores_per_item = floor(budget / n_items)

    This is the highest value where every item can run simultaneously
    in a single pool pass, capped by a per-stage maximum.

        n_workers = min(n_items, floor(budget / cores_per_item))

=========================================================================
STAGE CEILINGS  (STAGE_MAX_CORES)
=========================================================================

    orcacosmo    4   ORCA MPI overhead dominates above ~4 procs
    generation   4   CREST scales well to ~4
    optimisation 4   May use ORCA for opt; give headroom
    solubility   1   Pure Python openCOSMO-RS
    pruning      1   In-memory sort/filter
    cleaning     1   RDKit/pandas, GIL-bound

=========================================================================
WHAT PROCESSES THE RESULT
=========================================================================

    PipelineRunner._execute_pipeline
        Calls ResourceAllocator.allocate(stage, n_items) once per job.
        Injects result via AllocationResult.as_params() into
        job.parameters before job.run().

    BaseStage.__init__
        Reads injected fields from self.parameters and exposes them as
        instance attributes (self.n_workers, self.cores_per_item, etc.).
        Stages use these attributes directly — they never call the
        allocator themselves.

=========================================================================
USAGE
=========================================================================

    allocator = ResourceAllocator(
        config=job.config,
        parameters=req.parameters,
    )

    result = allocator.allocate(stage="orcacosmo", n_items=36)

    print(result.summary())
    # cores_per_item=3  n_workers=10  budget=32  free=38

    if result.warning:
        log("[WARNING] " + result.warning)

    # Inject into job parameters (done by PipelineRunner)
    job.parameters.update(result.as_params())

"""

import json
import os
from dataclasses import dataclass, field

try:
    import psutil
    _PSUTIL_AVAILABLE = True
except ImportError:
    _PSUTIL_AVAILABLE = False


# ── Stage configuration ───────────────────────────────────────────────────────

# Maximum cores_per_item per stage — reflects tool scaling behaviour
STAGE_MAX_CORES: dict[str, int] = {
    "orcacosmo":    8,
    "generation":   4,
    "optimisation": 8,
    "solubility":   1,
    "pruning":      1,
    "cleaning":     1,
}

# Stages where multi-core per item is meaningful
MULTI_CORE_STAGES: set[str] = {"orcacosmo", "generation", "optimisation"}


# ── Headroom rule ─────────────────────────────────────────────────────────────

def _apply_headroom(free_cores: int) -> int:
    """
    Apply the headroom rule to a raw free-core count.

    Steps:
      1. Subtract 2  (guaranteed minimum headroom)
      2. Round down to nearest valid multiple:
           <= 16 -> multiple of 2
           >  16 -> multiple of 4
      3. Floor at 1  (never return 0)

    Examples:
      40 -> 38 -> 36   (mult-of-4)
      41 -> 39 -> 36
      42 -> 40 -> 40
      18 -> 16 -> 16   (mult-of-2, on boundary)
      15 -> 13 -> 12
       4 ->  2 ->  2
       2 ->  0 ->  1   (floor)
       1 -> -1 ->  1   (floor)
    """
    after_headroom = free_cores - 2

    if after_headroom <= 1:
        return 1

    if after_headroom <= 16:
        return (after_headroom // 2) * 2
    else:
        return (after_headroom // 4) * 4


# ── Result dataclass ──────────────────────────────────────────────────────────

@dataclass
class AllocationResult:
    stage:          str
    n_items:        int
    cores_per_item: int
    n_workers:      int
    budget:         int                       # cores the pipeline is allowed to use
    free_cores:     int                       # sampled from psutil (0 if not sampled)
    budget_source:  str       = field(default="unknown")
    budget_reason:  str       = field(default="")
    stage_ceiling:  int       = field(default=1)
    warning:        str|None  = field(default=None)
    audit_trail:    list      = field(default_factory=list)

    def summary(self) -> str:
        """
        One-line string suitable for a [RESOURCES] log entry.

            cores_per_item=3  n_workers=10  budget=32  free=38
        """
        parts = [
            f"cores_per_item={self.cores_per_item}",
            f"n_workers={self.n_workers}",
            f"budget={self.budget}",
            f"free={self.free_cores}" if self.free_cores > 0 else "free=unknown",
        ]
        if self.warning:
            parts.append(f"WARN: {self.warning}")
        return "  ".join(parts)

    def detail_lines(self) -> list[str]:
        """
        Multi-line audit breakdown for the stage log.
        Written by PipelineRunner under a [RESOURCES] block.

        Example output:
          Budget source : runtime_override  (parameters[resources][cpus]=20)
          Stage ceiling : 4  (orcacosmo — ORCA MPI overhead dominates above ~4 procs)
          Items         : 2
          cores_per_item: 3  (floor(20 / 2) = 10, capped to stage ceiling 4 → 3... wait
                              floor(20 / 2) = 10, capped to 4)
          n_workers     : 5  (floor(20 / 4) = 5, capped to n_items 2 → 2)
          free_cores    : 38  (advisory — node appears healthy)
        """
        lines = []
        lines.append(f"  Budget source  : {self.budget_source}  ({self.budget_reason})")
        lines.append(f"  Budget         : {self.budget} cores")
        lines.append(f"  Stage ceiling  : {self.stage_ceiling} cores/item  ({self.stage} ceiling)")
        lines.append(f"  Items          : {self.n_items}")
        raw_cpi = max(1, self.budget // self.n_items) if self.n_items > 0 else 1
        lines.append(
            f"  cores_per_item : {self.cores_per_item}"
            f"  (floor({self.budget}/{self.n_items})={raw_cpi}"
            f", capped to ceiling {self.stage_ceiling})"
        )
        raw_workers = self.budget // self.cores_per_item if self.cores_per_item > 0 else 1
        lines.append(
            f"  n_workers      : {self.n_workers}"
            f"  (floor({self.budget}/{self.cores_per_item})={raw_workers}"
            f", capped to n_items {self.n_items})"
        )
        if self.free_cores > 0:
            health = "healthy" if self.free_cores >= self.budget else "UNDER LOAD"
            lines.append(f"  free_cores     : {self.free_cores}  (advisory — node {health})")
        else:
            lines.append(f"  free_cores     : unknown  (psutil not available or not sampled)")
        if self.warning:
            lines.append(f"  WARNING        : {self.warning}")
        return lines

    def as_params(self) -> dict:
        """
        Return a dict for injection into job.parameters.
        BaseStage reads these fields from self.parameters.
        """
        return {
            "cores_per_item": self.cores_per_item,
            "n_workers":      self.n_workers,
            "budget_cores":   self.budget,
        }


# ── Main class ────────────────────────────────────────────────────────────────

class ResourceAllocator:
    """
    Resolves core budget and computes per-stage allocation.

    Parameters
    ----------
    config : dict
        Pipeline config dict (from paths.json / job.config).
        Used to locate resource_allocation.json.

    parameters : dict
        Request/job parameters dict.
        Checked for {"resources": {"cores": N}} runtime override.
    """

    def __init__(self, config: dict, parameters: dict):
        self.config     = config or {}
        self.parameters = parameters or {}

        self._node_config = self._load_node_config()
        self._validate_node_config()
        self._budget      = self._resolve_budget()

    # ── Public API ────────────────────────────────────────────────────────────

    def allocate(self, stage: str, n_items: int) -> AllocationResult:
        """
        Compute allocation for a given stage and item count.

        Parameters
        ----------
        stage   : str   Stage name, e.g. "orcacosmo"
        n_items : int   Number of items to process

        Returns
        -------
        AllocationResult
        """
        if n_items <= 0:
            return AllocationResult(
                stage=stage, n_items=0,
                cores_per_item=1, n_workers=0,
                budget=self._budget, free_cores=0,
                warning=None,
            )

        # Advisory free-core check
        free_cores, warning = self._check_free_cores()

        # Cores per item and worker count
        cores_per_item = self._compute_cores_per_item(stage, n_items)
        n_workers      = min(n_items, max(1, self._budget // cores_per_item))

        # Budget source from trail (populated by _resolve_budget)
        trail         = getattr(self, "_budget_trail", [])
        budget_source = trail[0][0] if trail else "unknown"
        budget_reason = trail[0][1] if trail else ""
        stage_ceiling = STAGE_MAX_CORES.get(stage, 1)

        return AllocationResult(
            stage          = stage,
            n_items        = n_items,
            cores_per_item = cores_per_item,
            n_workers      = n_workers,
            budget         = self._budget,
            free_cores     = free_cores,
            budget_source  = budget_source,
            budget_reason  = budget_reason,
            stage_ceiling  = stage_ceiling,
            warning        = warning,
            audit_trail    = trail,
        )

    @property
    def budget(self) -> int:
        return self._budget

    @property
    def node_config(self) -> dict:
        return self._node_config

    # ── Config loading + validation ───────────────────────────────────────────

    def _load_node_config(self) -> dict:
        """Load resource_allocation.json via config."""
        ra_path = (
            self.config
            .get("resource_allocation", {})
            .get("defaults")
        )
        if ra_path and os.path.exists(ra_path):
            try:
                with open(ra_path) as f:
                    return json.load(f)
            except Exception as e:
                raise RuntimeError(
                    f"Failed to parse resource_allocation.json: {e}"
                ) from e
        return {}

    def _validate_node_config(self):
        """
        Validate that max_cores does not exceed total_cores.
        Raises ValueError at startup — not silently clamped.
        """
        total = self._node_config.get("total_cores")
        max_c = self._node_config.get("max_cores")

        if total is None or max_c is None:
            return  # nothing to validate

        try:
            total = int(total)
            max_c = int(max_c)
        except (TypeError, ValueError) as e:
            raise ValueError(
                f"resource_allocation.json: total_cores and max_cores "
                f"must be integers. Got total={total!r}, max={max_c!r}"
            ) from e

        if max_c > total:
            raise ValueError(
                f"resource_allocation.json: max_cores ({max_c}) "
                f"exceeds total_cores ({total}). "
                f"max_cores must be <= total_cores."
            )

    # ── Budget resolution ─────────────────────────────────────────────────────

    def _resolve_budget(self) -> int:
        """
        Resolve core budget.

        Hierarchy:
          1. parameters["resources"]["cpus"]  (or "cores")  runtime override
          2. node_config["max_cores"]          operator config
          3. psutil / os.cpu_count() fallback  headroom rule applied

        Decision trail is stored in self._budget_trail for logging.
        """
        self._budget_trail: list[tuple[str, str]] = []  # (source, reason)

        resources = self.parameters.get("resources") or {}

        # Accept both "cpus" (documented API) and "cores" (legacy compat)
        raw_cpus  = resources.get("cpus") or resources.get("cores")

        # 1. Runtime override
        if raw_cpus is not None:
            try:
                v = int(raw_cpus)
                if v > 0:
                    max_c = self._node_config.get("max_cores")
                    if max_c is not None and v > int(max_c):
                        raise ValueError(
                            f"parameters resources.cpus ({v}) exceeds "
                            f"resource_allocation.json max_cores ({max_c}). "
                            f"Lower your request or raise max_cores."
                        )
                    self._budget_trail.append((
                        "runtime_override",
                        f"parameters[resources][cpus]={v}"
                    ))
                    return v
            except (TypeError, ValueError) as e:
                if "exceeds" in str(e):
                    raise
                self._budget_trail.append((
                    "runtime_override_skipped",
                    f"non-integer value {raw_cpus!r}"
                ))

        # 2. Operator config
        try:
            v = int(self._node_config.get("max_cores") or 0)
            if v > 0:
                self._budget_trail.append((
                    "node_config",
                    f"resource_allocation.json max_cores={v}"
                ))
                return v
            else:
                self._budget_trail.append((
                    "node_config_skipped",
                    f"max_cores absent or zero"
                ))
        except (TypeError, ValueError):
            self._budget_trail.append((
                "node_config_skipped",
                "max_cores not parseable as int"
            ))

        # 3. Emergency fallback
        self._budget_trail.append((
            "fallback",
            "neither runtime override nor max_cores available"
        ))
        return self._fallback_budget()

    def _fallback_budget(self) -> int:
        """
        Emergency fallback when max_cores is not configured.
        Samples free cores via psutil and applies the headroom rule.
        Falls back to os.cpu_count() if psutil is unavailable.
        """
        if _PSUTIL_AVAILABLE:
            try:
                usage      = psutil.cpu_percent(interval=1, percpu=True)
                free_cores = sum(1 for u in usage if u < 20.0)
                return _apply_headroom(free_cores)
            except Exception:
                pass

        # psutil unavailable
        physical = os.cpu_count() or 2
        return _apply_headroom(physical)

    # ── Advisory free core check ──────────────────────────────────────────────

    def _check_free_cores(self) -> tuple[int, str | None]:
        """
        Advisory check of actual CPU availability.

        When max_cores IS configured:
            Sample psutil. Warn if free < budget but do not change budget.
            Returns (free_cores, warning_or_None).

        When max_cores is NOT configured:
            Free cores were already used to set the budget in
            _fallback_budget(). No second sample needed.
            Returns (0, None).
        """
        max_configured = bool(self._node_config.get("max_cores"))

        if not max_configured:
            return 0, None

        if not _PSUTIL_AVAILABLE:
            return 0, None

        try:
            usage      = psutil.cpu_percent(interval=1, percpu=True)
            free_cores = sum(1 for u in usage if u < 20.0)
        except Exception:
            return 0, None

        warning = None
        if free_cores < self._budget:
            warning = (
                f"only {free_cores} cores appear free "
                f"but budget is {self._budget} — "
                f"node may be under load"
            )

        return free_cores, warning

    # ── Allocation computation ────────────────────────────────────────────────

    def _compute_cores_per_item(self, stage: str, n_items: int) -> int:
        """
        Find highest cores_per_item such that all items fit in one pool pass.

            cores_per_item = floor(budget / n_items)

        Capped by STAGE_MAX_CORES.
        Stages not in MULTI_CORE_STAGES always return 1.
        """
        if stage not in MULTI_CORE_STAGES:
            return 1

        stage_max = STAGE_MAX_CORES.get(stage, 1)
        max_cpi   = max(1, self._budget // n_items)
        return min(max_cpi, stage_max)