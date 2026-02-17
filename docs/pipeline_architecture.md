# Pipeline Architecture  
openCOSMO‑RS Pipeline

This document describes the internal architecture of the openCOSMO‑RS Pipeline.  
It is intended for developers, maintainers, and advanced users who want to understand how the pipeline works internally, how stages communicate, and how reproducibility is guaranteed.

---

# 1. Architectural Overview

The pipeline is built around three core abstractions:

1. **Request**  
   Represents a full pipeline run. Immutable record of user intent.

2. **Job**  
   Represents a single stage execution within a Request.

3. **Stage**  
   The actual scientific computation (cleaning, generation, optimisation, etc.).

These components interact in a deterministic, config‑driven workflow:

```
Request → Job 0 → Job 1 → Job 2 → … → Job N
```

Each Job:
- Has a canonical input  
- Has a canonical output  
- Writes logs  
- Writes a job_state.json  
- Is resumable  

Each Stage:
- Implements scientific logic  
- Reads from canonical input  
- Writes canonical output  

---

# 2. Request

A Request is created when the user runs:

```
pl run
```

or programmatically via:

```python
Request.create_new(...)
```

A Request contains:

```
request.json
pipeline_state.json
request.log
jobs/
```

### request.json (immutable)
Contains:
- pipeline_spec  
- global config  
- stage arguments  
- dataset name  
- user  
- hostname  
- Git version  
- timestamps  

This file **must never be modified** after creation.

### pipeline_state.json (mutable)
Tracks:
- current stage  
- last completed stage  
- current job ID  
- request state (running, completed, failed)  
- halt reason (if any)  

This file is updated automatically by the pipeline.

---

# 3. Job

Each stage execution creates a Job directory:

```
jobs/J-<timestamp>-<stage>/
    inputs/
    outputs/
    job_state.json
    stage.log
```

### job_state.json
Tracks item‑level progress:

```
{
  "status": "running" | "completed" | "failed",
  "pending_items": [...],
  "completed_items": [...],
  "failed_items": [...],
  "created_at": "...",
  "completed_at": "...",
  "error_message": null
}
```

This enables:
- Resume  
- Partial progress  
- Deterministic recovery  
- Safe interruption (Ctrl+C)  

### Canonical Inputs & Outputs
Each stage declares its canonical output:

```python
self.set_stage_output("cleaned.csv")
```

The pipeline uses these canonical filenames to:
- Chain stages  
- Support continuation  
- Support reproducibility  

Canonical filenames are defined in:

```
Job.STAGE_OUTPUTS
```

---

# 4. Stage

Stages live in:

```
modules/stages/
```

Each stage subclasses `BaseStage` and implements:

```python
def execute(self):
    ...
```

Stages receive:

### `self.parameters`
Per‑stage arguments from `pipeline_spec`.

### `self.config`
Global configuration from `paths.json`.

### `self.job`
The Job object, containing:
- request_id  
- job_id  
- stage name  
- canonical paths  

### `self.input_path(name)`
Returns the canonical input path.

### `self.output_path(name)`
Returns the canonical output path.

### `self.update_progress(item)`
Moves an item from pending → completed.

### `self.fail(message)`
Marks the job as failed and stops execution.

---

# 5. Execution Flow

The pipeline is executed by:

```
PipelineRunner.run_request(request_id, base_dir)
```

Execution steps:

1. Load Request  
2. Determine resume index  
3. Load Job  
4. Run Stage  
5. Write canonical output  
6. Create next Job  
7. Repeat until final stage  

---

# 6. Resume Logic

Resume is fully deterministic.

### `_find_resume_index()` logic:

- If job_state.json missing → job never started → resume here  
- If status != completed → resume here  
- If all jobs completed → pipeline ends  

### Stage‑level resume
Stages use:

```
pending_items
completed_items
failed_items
```

to resume item‑level work.

### Example
If optimisation is interrupted:

```
pending_items = ["mol_27"]
completed_items = ["mol_1", ..., "mol_26"]
failed_items = []
```

Resume will:
- Skip cleaning  
- Skip generation  
- Skip pruning  
- Resume optimisation  
- Run only the pending item  

---

# 7. Continuation Logic

Continuation creates a **new Request** that uses the output of a previous Job as its input.

This is used for:
- Multi‑step workflows  
- Branching pipelines  
- Re‑running downstream stages with new parameters  

Continuation preserves:
- Provenance  
- Canonical outputs  
- Deterministic chaining  

---

# 8. Provenance

Every stage writes:

### stage.log
Detailed execution log.

### request.log
High‑level pipeline log.

### metadata provenance
Each metadata file includes:

```
pipeline_version
cleaning_timestamp
generated_by_request
generated_by_job
source_file
```

### Git version
Injected automatically:

```
main@abc1234
```

This ensures full reproducibility.

---

# 9. Deterministic Design Philosophy

The pipeline is intentionally:

### ✔ Config‑driven  
No hidden defaults.

### ✔ Deterministic  
Same input → same output.

### ✔ Canonical  
Each stage has a single canonical output.

### ✔ Reproducible  
Git version embedded in metadata.

### ✔ Resumable  
Item‑level resume for long stages.

### ✔ Extensible  
New stages can be added easily.

---

# 10. Adding a New Stage

To add a new stage:

1. Create `modules/stages/<name>_stage.py`
2. Subclass `BaseStage`
3. Implement `execute()`
4. Declare canonical output with `set_stage_output()`
5. Add stage to `Job.STAGE_OUTPUTS`
6. Add stage to `pipeline_spec`
7. Document the stage in `docs/stages/<name>.md`

---

# 11. Next Steps

Continue with:

- `docs/stages/cleaning.md`  
- `docs/stages/generation.md`  
- `docs/stages/pruning.md`  
- `docs/stages/optimisation.md`  
- `docs/stages/orcacosmo.md`  
- `docs/stages/solubility.md`  
- `docs/resume_and_continuation.md`  

