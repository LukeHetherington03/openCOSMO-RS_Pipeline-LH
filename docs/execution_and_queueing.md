# Execution & Queueing  
openCOSMO‑RS Pipeline

This document explains how to **create**, **run**, **resume**, and **queue** pipeline Requests using the openCOSMO‑RS Pipeline execution system.  
It focuses only on the commands required for execution.  
A full CLI reference is provided separately in `docs/cli_reference.md`.

---

# 1. Overview

The pipeline supports two execution modes:

1. **Direct execution**  
   Runs immediately in the foreground.

2. **Queued execution**  
   Adds a Request to a persistent queue.  
   A worker process picks up jobs and executes them.

Both modes begin the same way:

### ✔ A Request is created using a Python entrypoint  
### ✔ Execution mode is chosen inside the script  
### ✔ Queue management is handled separately using the CLI  

This ensures all Requests are created consistently and reproducibly.

---

# 2. Creating Requests (Python Entry Points)

Requests are created using the main Python scripts:

### **Create a new Request**
```
python3 main.py
```

### **Resume an interrupted Request**
```
python3 main_resume.py
```

### **Continue from a previous Request**
```
python3 main_continue.py
```

Each script:

- Loads `config/paths.json`
- Injects the Git version
- Builds the Request object
- Writes `request.json`
- Writes `pipeline_state.json`
- Creates the first Job
- Chooses execution mode based on `USE_QUEUE`

This is the **canonical** way to start pipeline runs.

---

# 3. Choosing Execution Mode (Inside main.py)

Inside `main.py`:

```python
USE_QUEUE = False   # Set to True to enqueue instead of running directly
```

Users choose execution mode by editing this line.

## **Direct Execution**
```
USE_QUEUE = False
```

The pipeline runs immediately:

- No queue involved  
- No worker required  
- Ideal for development, debugging, and small datasets  

## **Queued Execution**
```
USE_QUEUE = True
```

The Request is added to the queue:

- Worker processes handle execution  
- Ideal for HPC clusters and long‑running workloads  
- Allows multiple Requests to be queued safely  

---

# 4. Queue Management (Required CLI Commands Only)

To use queued execution, users must manage the queue using the CLI.

Add the alias:

```bash
alias pl='python3 -m modules.cli.cli'
```

Reload:

```bash
source ~/.bashrc
```

### **Start the worker**
```
pl q start
```

### **Stop the worker**
```
pl q stop
```

### **Check queue status**
```
pl q status
```

### **List queued Requests**
```
pl q list
```

These are the only CLI commands required for execution.  
All other CLI commands are documented separately.

---

# 5. Inspecting Requests (Minimal Required CLI)

During execution, users may want to inspect Requests.

### **Check Request status**
```
pl r status <id>
```

### **View logs**
```
pl r logs <id>
```

These are the only Request commands needed during execution.  
Full Request management commands are documented separately.

---

# 6. Environment Validation (Minimal Required CLI)

Before running the pipeline, users should validate their environment:

### **Validate required Python packages**
```
pl env pip
```

### **Validate external executables**
```
pl env software
```

### **Validate openCOSMO paths + constant files**
```
pl env resources
```

These are the only environment commands required for execution.  
Full environment commands are documented separately.

---

# 7. Execution Flow (Direct or Queued)

Regardless of execution mode, the pipeline follows the same internal flow:

1. Load Request  
2. Determine resume index  
3. Load Job  
4. Run Stage  
5. Write canonical output  
6. Create next Job  
7. Update pipeline_state.json  
8. Continue until final stage  

Stages never:

- Guess inputs  
- Guess outputs  
- Write outside their job directory  

This ensures full reproducibility.

---

# 8. Resume Logic

Resume is based on:

### `pipeline_state.json`
Tracks:

- current_stage  
- current_job_id  
- last_completed_stage  

### `job_state.json`
Tracks:

- pending_items  
- completed_items  
- failed_items  
- status  

Resume rules:

- If job_state.json missing → job never started → resume here  
- If status != completed → resume here  
- If all jobs completed → pipeline ends  

This makes resume:

- Deterministic  
- Safe  
- Reproducible  

---

# 9. Continuation Logic

Continuation creates a **new Request** with:

- New request_id  
- New job lineage  
- Stage 0 input = previous job’s canonical output  

Continuation is used for:

- Multi‑step workflows  
- Parameter sweeps  
- Re‑running downstream stages  
- Solubility‑only runs  
- ORCA‑only runs  

---

# 10. Worker Process

A worker continuously processes queued Requests:

```
pl q start
```

The worker:

1. Waits for a queued Request  
2. Locks it  
3. Runs the pipeline  
4. Updates queue state  
5. Moves to the next Request  

Workers are:

- Safe  
- Multi‑instance aware  
- Crash‑resistant  
- Deterministic  

---

# 11. HPC Usage

The pipeline is HPC‑friendly:

- No environment modules required  
- No shared state  
- No race conditions  
- Workers can run on multiple nodes  
- Resume works across nodes  
- Continuation works across nodes  

Recommended:

- Run workers via SLURM or PBS  
- Use `pl q start` inside a job script  
- Use `USE_QUEUE=True` in main.py to submit work  

---

# 12. Summary

The execution system provides:

- Direct execution  
- Resume  
- Continuation  
- Queueing  
- Worker processes  
- Deterministic chaining  
- Full provenance  
- HPC‑safe operation  

Requests are always created via Python entrypoints.  
Queue management and inspection are handled via the CLI.

---

# 13. Next Steps

Continue with:

- `docs/cli_reference.md`  
