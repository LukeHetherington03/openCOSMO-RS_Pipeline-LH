# Developer Guide

openCOSMO-RS Pipeline

---

## Contents

1. [Architecture overview](#1-architecture-overview)
2. [Code change / pycache note](#2-code-change--pycache-note)
3. [Adding a new stage](#3-adding-a-new-stage)
4. [Adding a new generation backend](#4-adding-a-new-generation-backend)
5. [Adding a new optimisation engine](#5-adding-a-new-optimisation-engine)
6. [BaseStage contract](#6-basestage-contract)
7. [Worker function contract](#7-worker-function-contract)
8. [Logging conventions](#8-logging-conventions)
9. [Job and Request data model](#9-job-and-request-data-model)
10. [Pipeline runner and resource allocation](#10-pipeline-runner-and-resource-allocation)
11. [Queue system](#11-queue-system)

---

## 1. Architecture overview

```
CLI (pl)
  └── QueueCommands / RequestCommands / EnvCommands
        └── PipelineRunner.submit()  →  queue (SQLite)

QueueWorker (background process)
  └── PipelineRunner.run()
        ├── ResourceAllocator       allocates cores/workers per stage
        ├── Job.create_new()        creates job directory + state
        └── Job.run()
              └── {Stage}Stage(job).run()
                    └── execute()   stage-specific logic
                          └── run_parallel(worker_fn, build_args_fn)
                                └── ProcessPoolExecutor
                                      └── module-level worker function
```

### Key directories

```
modules/
    build/
        job_manager.py          Job class — directory layout, state, lifecycle
        request_manager.py      Request class — top-level container
    cli/
        cli.py                  Entry point; delegates to q/r/env command classes
        q_commands.py           Queue control commands
        r_commands.py           Request inspection commands
        env_commands.py         Environment validation
    execution/
        runner.py               PipelineRunner — orchestrates stages end-to-end
        resource_allocator.py   Decides n_workers and cores_per_item per stage
        queue_worker.py         Background worker process
    stages/
        base_stage.py           BaseStage — lifecycle, parallel execution, logging
        cleaning_stage.py
        generation_stage.py
        pruning_stage.py
        optimisation_stage.py
        orcacosmo_stage.py
        solubility_stage.py
    utils/
        log_helper.py           Structured log writer (all levels)
        atomic_write.py         Atomic file writer for checkpoints

config/
    paths.json                  Software/resource paths
    resource_allocation.json    Core budget ceiling
    generation_defaults.json
    pruning_defaults.json
    optimisation_defaults.json
    optimisation_engines.json   Available optimisation engine definitions
    orcacosmo_defaults.json
    solubility_defaults.json

CONSTANT_FILES/
    molecule_metadata/          Pre-computed per-molecule metadata (InChIKey keyed)
    solvents/                   Per-solvent COSMO files
    solvent_lists/              Named solvent lists (JSON)
    chemistry/                  Functional group definitions etc.
```

### Data flow between stages

Each stage produces exactly one **canonical output file** and passes it to the next stage via `stage_input` in the job parameters.

| Stage | Canonical output |
|-------|-----------------|
| cleaning | `cleaned.csv` |
| generation | `energies.json` |
| pruning | `energies.json` |
| optimisation | `energies.json` |
| orcacosmo | `orcacosmo_summary.json` |
| solubility | `solubility_results.json` |

The `energies.json` format is a **ConformerSet JSON** — a list of conformer records that accumulates provenance and optimisation history through the pipeline.

---

## 2. Code change / pycache note

The queue worker is a long-running background process. Python caches compiled bytecode in `__pycache__/`. If you change source files while the worker is running, the worker continues executing the old cached version.

**You must restart the worker after any code change:**

```bash
pl q stop    # wait for running items to finish
pl q start
```

Or for an immediate restart (in-flight items will re-run from checkpoints on resume):

```bash
pl q kill
pl q start
```

---

## 3. Adding a new stage

1. **Create** `modules/stages/<name>_stage.py` with a class named `<Name>Stage` extending `BaseStage`.

2. **Implement `execute()`** — this is the only method stages must define. Do not call `mark_running()`, `mark_complete()`, or `mark_failed()` — BaseStage.run() owns the lifecycle.

3. **Register the canonical output** in `Job.STAGE_OUTPUTS` in `modules/build/job_manager.py`:
   ```python
   STAGE_OUTPUTS = {
       ...
       "<name>": "<name>_output.json",
   }
   ```

4. **Add config defaults** (optional) in `config/<name>_defaults.json` and reference it in `config/paths.json`.

5. **Add extra output files** to `_STAGE_EXTRA_FILES` in `modules/execution/runner.py` if the stage produces auxiliary outputs (CSV, txt) that should be copied to `final_outputs/`:
   ```python
   _STAGE_EXTRA_FILES = {
       ...
       "<name>": ["<name>_results.csv", "<name>_summary.txt"],
   }
   ```

6. The stage can then be referenced by name in `pipeline_sequence` in a `request.json`.

### Minimal stage skeleton

```python
class MyStage(BaseStage):

    def execute(self):
        self._load_stage_config()
        # ... read self.parameters["stage_input"]
        # ... process items
        # ... write canonical output to self.outputs_dir
```

---

## 4. Adding a new generation backend

Generation backends are module-level functions (picklable for multiprocessing) dispatched through a single entry point.

### Steps

1. **Register** the backend name in `GENERATION_BACKENDS`:
   ```python
   GENERATION_BACKENDS = {
       "rdkit":     "_worker_rdkit",
       "crest":     "_worker_crest",
       "openbabel": "_worker_openbabel",
       "mybackend": "_worker_mybackend",   # add here
   }
   ```

2. **Add a dispatch branch** in `_generation_worker()`:
   ```python
   def _generation_worker(args):
       backend = args["engine"]
       if backend == "mybackend":
           return _worker_mybackend(args)
       ...
   ```

3. **Implement the worker function** at module level (not a method):
   ```python
   def _worker_mybackend(args: dict) -> dict:
       """
       args contains (from _base_args + _build_args):
           item_id, inchi_key, smiles, charge, multiplicity,
           seed, cores_per_item, num_confs,
           xyz_out_dir, workdir, log_path, checkpoints_dir,
           timestamp, stage, outputs_dir
           # plus any tool paths you added in _build_args()
       """
       # ... generate conformers ...

       return {
           "item_id":     args["item_id"],
           "status":      "ok",           # or "failed"
           "log_details": [               # shown in [COMPLETE] log line
               ("n_confs", str(n_generated)),
               ("method",  "mybackend"),
           ],
           "worker_name": multiprocessing.current_process().name,
           "worker_pid":  os.getpid(),
           "conformers":  [c.to_dict() for c in conformer_list],
           # "error": "message"  # required when status == "failed"
       }
   ```

4. **Add resource defaults** in `config/generation_defaults.json`:
   ```json
   {
     "cores_per_item": {
       "mybackend": 2
     },
     "parallel": {
       "mybackend": true
     }
   }
   ```

5. **Resolve executable paths** in `GenerationStage._build_args()` if your backend needs a tool binary.

6. **Add version detection** in `GenerationStage._get_tool_versions()` if desired.

---

## 5. Adding a new optimisation engine

Optimisation engines are defined declaratively in `config/optimisation_engines.json`. No Python code change is required for engines that fit an existing family (`orca`, `xtb`, `gxtb`, `forcefield`).

### Engine families

| Family | Required fields | Description |
|--------|----------------|-------------|
| `orca` | `method`, `basis`, `opt`/`sp`, `cpcm`, `alpb` | DFT via ORCA |
| `xtb` | `gfn`, `level` | GFN-xTB via XTB binary |
| `gxtb` | `gfn`, `level`, `opt_flag` | GFN-xTB via gXTB binary |
| `forcefield` | `forcefield` | RDKit force field (MMFF94, UFF) |

### Example: add a new gXTB tight engine

```json
"gxtb_opt_tight_gfn1": {
    "family": "gxtb",
    "gfn": 1,
    "level": "tight",
    "opt_flag": "--opt tight"
}
```

### Adding a new engine family

If none of the existing families fit your backend:

1. Add the engine definition to `optimisation_engines.json` with `"family": "myfamily"`.
2. Add a dispatch branch in `OptimisationStage._run_engine()` (or the equivalent backend runner).
3. Implement the runner as a module-level function following the worker function contract (see section 7).

---

## 6. BaseStage contract

All stage classes inherit from `BaseStage` and must implement `execute()`.

### Attributes available to all stages

| Attribute | Source | Description |
|-----------|--------|-------------|
| `self.job` | Constructor | The `Job` instance |
| `self.parameters` | `job.parameters` | Runtime parameters (stage_args merged with defaults) |
| `self.config` | `job.config` | Global `config/paths.json` contents |
| `self.inputs_dir` | `job.inputs_dir` | `jobs/<job_id>/inputs/` |
| `self.outputs_dir` | `job.outputs_dir` | `jobs/<job_id>/outputs/` |
| `self.n_workers` | Injected | Number of parallel workers for this stage |
| `self.cores_per_item` | Injected | CPU cores per item (passed to tools) |
| `self.budget_cores` | Injected | Total core budget (`n_workers * cores_per_item`) |
| `self.checkpoints_dir` | BaseStage | `outputs/checkpoints/` |
| `self.strict_mode` | Stage config | Abort pool on first failure |

### Lifecycle (do not override `run()`)

```
BaseStage.run()
  _write_debug_context()      writes stage_context.log once at startup
  stage_started log
  job.mark_running()
  → self.execute()            STAGES IMPLEMENT THIS ONLY
  stage_complete log
  job.mark_complete()
  request.update_pipeline_state()
```

If `execute()` raises, `run()` catches it, logs `[STAGE] FAILED`, calls `job.mark_failed()`, and re-raises.

### Checkpoint API

```python
self._write_checkpoint(item_id, result_dict)  # write atomically
self._load_checkpoint(item_id)                # -> dict | None
self._load_all_checkpoints()                  # -> list[dict]
self._clear_checkpoints()                     # delete all checkpoints
```

Checkpoints allow stages to resume after a crash or stop without re-running completed items.

### Logging API

```python
self.log(msg)              # [INFO] to stage log
self.log_config(msg)       # [CONFIG] to stage log
self.log_info(msg)         # [INFO] to stage log (same as log())
self.log_warning(msg)      # [WARNING] — also captured in pipeline warnings summary
self.log_error(msg)        # [ERROR] to stage log
self.log_skip(item_id, reason)   # [SKIP] for deliberately skipped items
```

Only `log_warning()` lines appear in the pipeline warnings summary and `final_outputs/pipeline_warnings.txt`. Use it for anything a user needs to review (defaults, missing inputs, fallbacks).

---

## 7. Worker function contract

Worker functions are module-level functions (not methods) that run in child processes via `ProcessPoolExecutor`. They must be picklable.

### Required return dict

```python
{
    "item_id":     str,       # must match the item passed in
    "status":      str,       # "ok" or "failed"
    "log_details": list,      # list of (key, value) tuples shown in [COMPLETE] line
    # Required when status == "failed":
    "error":       str,
    # Stage-specific payload (written to checkpoint):
    ...
}
```

`BaseStage.run_parallel()` strips `log_details` before writing the checkpoint to disk. Everything else is preserved in the checkpoint.

### What BaseStage does for each result

- **ok**: logs `[COMPLETE]`, writes checkpoint, calls `job.update_progress(item, success=True)`
- **failed**: logs `[FAILED]` with `error`, calls `job.update_progress(item, success=False)`
- **strict_mode on + failed**: cancels all pending futures, shuts down the pool

### `_base_args(item)` — standard args

Every `build_args_fn` should call `self._base_args(item)` and extend the result:

```python
def _build_args(self, item):
    args = self._base_args(item)
    args.update({
        "my_param": self.my_param,
        "tool_exe": self.config["mytool"]["executable"],
    })
    return args
```

`_base_args` always provides: `item_id`, `stage`, `outputs_dir`, `workdir`, `cores_per_item`, `log_path`, `checkpoints_dir`.

---

## 8. Logging conventions

All logging goes through `LogHelper` in `modules/utils/log_helper.py`.

### Level taxonomy

| Level | Used for |
|-------|---------|
| `[INFO]` | General progress — file counts, stage stages |
| `[CONFIG]` | Resolved config values at startup |
| `[RESOURCES]` | Core/worker allocation decisions |
| `[START]` | Item picked up by a worker |
| `[COMPLETE]` | Item finished — with detail block |
| `[FAILED]` | Item failed — with error block |
| `[RESUME]` | Item skipped (checkpoint exists) |
| `[WARNING]` | Non-fatal issue — defaults, fallbacks, missing optional data |
| `[ERROR]` | Stage-level error before raise |
| `[POOL]` | Pool started / drained |
| `[STAGE]` | Stage started / complete / failed |
| `[SKIP]` | Item deliberately skipped |
| `[DEBUG]` | Verbose detail (debug log only, not main log) |

### Log files per job

```
jobs/<job_id>/
    stage_logs/
        stage.log          main stage log — all levels above
        stage_context.log  debug context snapshot (startup only)
```

The root-level `stage.log` (directly in `job_dir`) is **not** used — always read from `stage_logs/`.

### Pipeline warnings

The pipeline scanner reads `[WARNING]` lines from stage logs after the pipeline completes and collates them into:
- A `[PIPELINE WARNINGS]` block in `request.log`
- `final_outputs/pipeline_warnings.txt`

Warnings are grouped by section (stage name extracted from the message) and deduplicated within each section.

---

## 9. Job and Request data model

### Request

```
pipeline_data/requests/<request_id>/
    request.json          submitted pipeline spec
    request_state.json    pipeline state (stage, status, timestamps)
    request.log           high-level request log
    jobs/
        <job_id>/
            job_state.json      job status, item counts, parameters
            inputs/
                parameters.json   merged parameters for this job
                <stage_input>     canonical output copied from previous stage
            outputs/
                <canonical_output>
                checkpoints/
                ...
            stage_logs/
                stage.log
                stage_context.log
    final_outputs/
        ...
```

### Job state fields (`job_state.json`)

| Field | Description |
|-------|-------------|
| `job_id` | Unique job identifier (`J-<timestamp>-<stage>`) |
| `stage` | Stage name |
| `status` | `initialised` / `running` / `completed` / `failed` |
| `items` | All item IDs for this job |
| `pending_items` | Items not yet processed |
| `completed_items` | Successfully processed items |
| `failed_items` | Failed items |
| `parameters` | Full merged parameters dict |
| `output_file` | Canonical output path (set on completion) |

### Pipeline state (`request_state.json`)

Tracks `state` (`running`/`completed`/`failed`), `current_stage`, `last_completed_stage`, `completed_jobs` (list of `[job_id, stage]` pairs), `halt_reason`.

---

## 10. Pipeline runner and resource allocation

`PipelineRunner` in `modules/execution/runner.py` is the central orchestrator.

### Resource allocation

Before each stage, `ResourceAllocator.allocate()` computes:

- `n_items` — number of items the stage will process
- `max_cores` — from `config/resource_allocation.json` (or `parameters.resources.cpus`)
- `cores_per_item` — from stage defaults or `parameters.cores_per_item`
- `n_workers` — `floor(max_cores / cores_per_item)`, capped at `n_items`

The computed values are injected into `parameters` as `n_workers`, `cores_per_item`, `budget_cores` before `job.run()` is called.

If `max_cores` cannot accommodate even one full `cores_per_item` budget, a warning is emitted and `cores_per_item` is reduced to fit.

### Stage chaining

After each stage completes, the canonical output file is copied into the next job's `inputs/` directory and registered as `stage_input` in that job's parameters. This copy is deterministic — there is no path guessing.

---

## 11. Queue system

The queue is backed by a SQLite database. Requests are submitted with a priority (default 5, lower = higher priority).

### Queue worker lifecycle

```
pl q start   →  QueueWorker starts as background process
                Polls queue every N seconds
                Picks highest-priority pending request
                Calls PipelineRunner.run(request)
                Loops

pl q stop    →  Sends SIGTERM to worker
                Worker cancels pending futures (in-flight items finish)
                SIGTERM propagated to child process groups
                Worker exits cleanly

pl q kill    →  Sends SIGKILL to worker and all children
                In-flight items will re-run from checkpoints on resume
```

### Resume

If a request was interrupted (stop/kill/crash), resubmitting it resumes from the last completed stage. Items that wrote checkpoints before interruption are skipped (they emit `[RESUME]`).
