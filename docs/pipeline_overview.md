# Pipeline Overview  
openCOSMO‑RS Pipeline

The openCOSMO‑RS Pipeline is a fully reproducible, config‑driven workflow for computing molecular solubilities using modern conformer generation, quantum chemistry, and COSMO‑RS modelling.  
This document provides a high‑level overview of the pipeline architecture, execution model, and scientific flow.

---

# 1. Purpose of the Pipeline

The pipeline automates the full solubility‑prediction workflow:

1. **Data ingestion**  
2. **Conformer generation**  
3. **Conformer pruning**  
4. **Geometry optimisation**  
5. **ORCA COSMO surface generation**  
6. **COSMO‑RS solubility prediction**

It is designed to be:

- **Reproducible** — every run is fully provenance‑tracked  
- **Deterministic** — canonical inputs/outputs at every stage  
- **Config‑driven** — no hard‑coded paths or parameters  
- **Resumable** — safe recovery from interruption  
- **Queue‑friendly** — supports HPC and long‑running workloads  
- **Extensible** — new stages can be added easily  

---

# 2. Architectural Overview

The pipeline is built around three core abstractions:

### **Request**
A full pipeline run.  
Contains:

- `request.json` (immutable user intent)  
- `pipeline_state.json` (mutable execution state)  
- A sequence of Jobs  

### **Job**
A single stage execution.  
Contains:

- Canonical inputs  
- Canonical outputs  
- `job_state.json`  
- Logs  
- Resume information  

### **Stage**
A scientific computation implemented as a Python class.  
Stages inherit from `BaseStage`, which provides:

- Lifecycle management  
- Logging  
- Item‑level progress tracking  
- Canonical path helpers  
- Strict‑mode enforcement  

Stages never write outside their job directory and never guess inputs or outputs.

---

# 3. Execution Model

Requests are created using Python entrypoints:

```
python3 main.py
python3 main_resume.py
python3 main_continue.py
```

Execution mode is chosen **inside the script**:

- **Direct execution** — run immediately  
- **Queued execution** — submit to the queue  

Queue management is handled via the CLI:

```
pl q start
pl q status
pl q stop
```

Request inspection is also done via the CLI:

```
pl r status <id>
pl r logs <id>
```

This separation ensures reproducibility and clarity.

---

# 4. Data Flow Through the Pipeline

The pipeline processes data through six scientific stages.  
Each stage consumes the canonical output of the previous one.

---

## **Stage 1 — Cleaning**

Input: raw CSV files  
Output: `cleaned.csv` + metadata JSONs  

Responsibilities:

- Standardise headers  
- Validate required fields  
- Canonicalise SMILES  
- Generate InChIKeys  
- Compute physchem descriptors  
- Write molecule metadata  

This stage establishes the molecular identity and metadata used throughout the pipeline.

---

## **Stage 2 — Generation**

Input: `cleaned.csv`  
Output: `energies.json` + XYZ files  

Responsibilities:

- Generate 3D conformers (RDKit, CREST, OpenBabel)  
- Compute generation‑level energies  
- Write XYZ files  
- Record provenance  

This stage creates the initial conformer ensemble.

---

## **Stage 3 — Pruning**

Input: `energies.json`  
Output: pruned `energies.json`  

Responsibilities:

- Remove invalid conformers  
- Apply RMSD/energy/percentile/N‑based pruning  
- Produce a compact, diverse conformer set  

This stage reduces computational cost for optimisation.

---

## **Stage 4 — Optimisation**

Input: pruned `energies.json`  
Output: optimised `energies.json`  

Responsibilities:

- Geometry optimisation (ORCA, gXTB, XTB, forcefield)  
- Checkpointing and resume  
- Convergence detection  
- Energy validation  
- Provenance tracking  

This stage produces high‑quality geometries for COSMO.

---

## **Stage 5 — ORCACOSMO**

Input: optimised `energies.json`  
Output: `orcacosmo_summary.json` + `.orcacosmo` files  

Responsibilities:

- Run ORCA CPCM single‑point calculations  
- TZVPD → TZVP fallback  
- Parse `.log`, `.cpcm`, `.cpcm_corr`  
- Reconstruct `.orcacosmo` files  
- Write raw + parsed outputs  

This stage generates COSMO surfaces for COSMO‑RS.

---

## **Stage 6 — Solubility**

Input: `orcacosmo_summary.json`  
Output: `solubility_results.json`  

Responsibilities:

- Build COSMO‑RS mixture inputs  
- Run COSMO‑RS engine  
- Write per‑molecule solubility predictions  
- Produce human‑readable summaries  

This stage produces the final scientific output.

---

# 5. Provenance & Reproducibility

The pipeline embeds provenance at every level:

- Git version  
- Request ID  
- Job ID  
- Timestamps  
- Engine versions  
- Input file paths  
- Optimisation history  
- COSMO reconstruction metadata  

Every output can be traced back to its exact inputs and parameters.

---

# 6. Resume & Continuation

### **Resume**
Recover an interrupted Request:

- Uses `pipeline_state.json`  
- Uses `job_state.json`  
- Resumes at item level  

### **Continuation**
Start a new Request using the output of a previous one:

- New request_id  
- New job lineage  
- Same canonical input  

This enables branching workflows and downstream re‑runs.

---

# 7. Configuration Philosophy

The pipeline is entirely config‑driven:

- Executable paths  
- CONSTANT_FILES  
- COSMO radii  
- Optimisation engines  
- Defaults for solubility  

Users modify configuration, not code.

---

# 8. Directory Structure

Each Request has:

```
request.json
pipeline_state.json
request.log
jobs/
```

Each Job has:

```
inputs/
outputs/
job_state.json
stage.log
```

This structure is deterministic and reproducible.

---

# 9. Extending the Pipeline

To add a new stage:

1. Create `modules/stages/<name>_stage.py`  
2. Subclass `BaseStage`  
3. Implement `execute()`  
4. Declare canonical output  
5. Add to `pipeline_spec`  
6. Document in `docs/stages/<name>.md`  

The architecture is intentionally modular.

---

# 10. Summary

The openCOSMO‑RS Pipeline provides:

- A complete end‑to‑end solubility workflow  
- Deterministic, reproducible scientific computation  
- Robust resume and continuation  
- HPC‑friendly queueing  
- Rich provenance  
- Modular, extensible architecture  

It is designed for both production‑grade solubility prediction and research‑grade experimentation.

