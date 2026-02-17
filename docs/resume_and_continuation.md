# Resume & Continuation  
openCOSMO‑RS Pipeline

This document explains how the openCOSMO‑RS Pipeline implements **resume** and **continuation**, two core mechanisms that ensure robustness, reproducibility, and efficient long‑running execution.

Resume and continuation are related but fundamentally different:

- **Resume** → Continue an *interrupted* Request  
- **Continuation** → Create a *new* Request using the output of a previous one  

Both mechanisms are deterministic and fully provenance‑tracked.

---

# 1. Overview

The pipeline is designed for:

- Long‑running scientific computations  
- HPC preemption  
- Partial failures  
- Multi‑stage workflows  
- Branching workflows  

To support these use cases, the pipeline provides:

### ✔ Resume  
Recover from interruption without repeating completed work.

### ✔ Continuation  
Start a new Request using the canonical output of a previous Request.

Both mechanisms rely on:

- Canonical inputs/outputs  
- Job‑level state (`job_state.json`)  
- Request‑level state (`pipeline_state.json`)  
- Deterministic stage chaining  

---

# 2. Resume

Resume is used when a Request was **interrupted** before completion.

Typical causes:

- User pressed Ctrl+C  
- HPC job timeout  
- Node failure  
- Worker shutdown  
- Unexpected error  

Resume ensures:

- No duplicated work  
- No corrupted state  
- No skipped items  
- Deterministic recovery  

---

# 3. How Resume Works

Resume is invoked via:

```
python3 main_resume.py
```

The script:

1. Loads the Request  
2. Reads `pipeline_state.json`  
3. Determines the stage to resume  
4. Loads the corresponding Job  
5. Reads `job_state.json`  
6. Resumes only the **pending** items  
7. Continues the pipeline normally  

### 3.1 pipeline_state.json

Tracks:

```
current_stage
current_job_id
last_completed_stage
state
```

### 3.2 job_state.json

Tracks:

```
pending_items
completed_items
failed_items
status
```

### 3.3 Resume Rules

The pipeline resumes at the **first incomplete job**:

- If `job_state.json` is missing → job never started → resume here  
- If `status != completed` → resume here  
- If all jobs completed → pipeline ends  

### 3.4 Item‑Level Resume

Stages that operate on multiple items (e.g., molecules) resume at the item level:

- Completed items are skipped  
- Pending items are processed  
- Failed items are retried or reported  

This makes resume efficient and safe.

---

# 4. Continuation

Continuation is used when the user wants to:

- Start a new Request  
- Use the output of a previous Request  
- Change parameters for downstream stages  
- Branch a workflow  
- Perform solubility‑only or ORCA‑only runs  
- Re‑run later stages without repeating earlier ones  

Continuation is **not** recovery — it is **workflow branching**.

---

# 5. How Continuation Works

Continuation is invoked via:

```
python3 main_continue.py
```

The script:

1. Loads the previous Request  
2. Identifies the canonical output of the final completed stage  
3. Creates a **new Request**  
4. Sets the new Request’s stage 0 input to the previous output  
5. Writes a new `request.json`  
6. Writes a new `pipeline_state.json`  
7. Creates a new Job 0  
8. Runs or queues the new Request  

### 5.1 New Request, New Provenance

Continuation creates a **completely new Request**, with:

- New request_id  
- New job lineage  
- New metadata  
- New logs  
- New pipeline_state.json  

But it also records:

```
"continued_from": "<previous_request_id>"
```

This preserves provenance.

---

# 6. Differences Between Resume and Continuation

| Feature | Resume | Continuation |
|--------|--------|--------------|
| Purpose | Recover from interruption | Start a new workflow branch |
| Creates new Request? | No | Yes |
| Repeats completed work? | No | No |
| Uses previous output? | Yes (same Request) | Yes (new Request) |
| Provenance | Same Request | New Request with link to old |
| Typical use | HPC preemption | Running solubility only, ORCA only, parameter sweeps |

---

# 7. Canonical Outputs Enable Both

Both resume and continuation rely on **canonical outputs**, defined in:

```
Job.STAGE_OUTPUTS
```

Examples:

- `cleaned.csv`
- `energies.json`
- `orcacosmo_summary.json`
- `solubility_results.json`

Because every stage produces a deterministic, canonical output:

- Resume knows exactly where to continue  
- Continuation knows exactly what to use as input  

This is the backbone of reproducibility.

---

# 8. Example Workflows

## 8.1 Resume After Interruption

User starts a run:

```
python3 main.py
```

HPC node times out during optimisation.

User resumes:

```
python3 main_resume.py
```

Pipeline resumes at:

- Stage: optimisation  
- Item: first pending molecule  

---

## 8.2 Continue to Solubility Only

User completes optimisation + ORCA COSMO.

Then runs:

```
python3 main_continue.py
```

The new Request:

- Uses `orcacosmo_summary.json` as input  
- Runs only the solubility stage  

---

## 8.3 Branching Workflow

User wants:

- One branch with XTB optimisation  
- One branch with ORCA optimisation  

Workflow:

```
python3 main.py              # branch A
python3 main_continue.py     # branch B
```

Branch B modifies the pipeline_spec to use ORCA.

---

# 9. Safety & Determinism

Resume and continuation are designed to be:

### ✔ Deterministic  
Same input → same output.

### ✔ Safe  
No overwriting of previous results.

### ✔ Reproducible  
Git version embedded in metadata.

### ✔ Transparent  
All decisions logged in:

- request.log  
- stage.log  
- job_state.json  
- pipeline_state.json  

---

# 10. Summary

- **Resume** recovers an interrupted Request  
- **Continuation** creates a new Request using previous output  
- Both rely on canonical outputs  
- Both preserve provenance  
- Both are deterministic  
- Both are essential for HPC and long‑running workflows  

---

# 11. Related Documents

- `docs/execution_and_queueing.md`  
- `docs/pipeline_architecture.md`  
- `docs/configuration.md`  
- `docs/stages/cleaning.md`  
- `docs/stages/optimisation.md`  

