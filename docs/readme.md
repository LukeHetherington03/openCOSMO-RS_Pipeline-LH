# openCOSMO‑RS Pipeline

A fully reproducible, deterministic, and extensible scientific workflow engine for:

- Molecular cleaning & metadata generation  
- Conformer generation  
- Conformer pruning  
- Geometry optimisation (XTB, gXTB, ORCA)  
- ORCA COSMO calculations  
- COSMO‑RS solubility prediction  

The pipeline is designed for **robustness**, **traceability**, and **scientific reproducibility**, with:

- Canonical stage inputs/outputs  
- Deterministic job chaining  
- Full provenance tracking  
- Resume & continuation support  
- Config‑driven execution  
- Git version embedding  

---

## 🚀 Quickstart

### 1. Prepare your environment

Install dependencies (RDKit, ORCA, XTB, gXTB, CREST, COSMO‑RS bindings).  
See `docs/installation.md` for full details.

### 2. Configure paths

Edit:

```
config/paths.json
```

to point to:

- ORCA executable  
- XTB/gXTB executables  
- CREST  
- COSMO‑RS Python + C++ bindings  
- CONSTANT_FILES directory  
- pipeline_data directory  

### 3. Run the pipeline

```bash
python3 -m modules.main
```

This will:

- Create a new Request  
- Create the first Job  
- Execute each stage sequentially  
- Produce a full provenance trail  

### 4. Resume a request

If a job was interrupted (Ctrl+C, crash, HPC preemption):

```bash
python3 -m modules.main_resume
```

The pipeline will:

- Skip completed stages  
- Resume the interrupted job  
- Continue normally  

### 5. Continue from a previous request

To start a new pipeline using the output of a previous job:

```bash
python3 -m modules.main_continue
```

---

## 📂 Directory Structure

```
pipeline_data/
    requests/
        R-<timestamp>-<title>/
            request.json
            pipeline_state.json
            request.log
            jobs/
                J-<timestamp>-<stage>/
                    inputs/
                    outputs/
                    job_state.json
                    stage.log
```

### Key files:

- **request.json** — immutable record of user intent  
- **pipeline_state.json** — current stage, last completed stage  
- **job_state.json** — item‑level progress (pending/completed/failed)  
- **stage.log** — detailed execution log  
- **canonical outputs** — e.g., `cleaned.csv`, `energies.json`, `orcacosmo_summary.json`  

---

## 🧠 Pipeline Stages

### 1. Cleaning

- Reads raw CSVs  
- Standardises headers  
- Canonicalises SMILES  
- Generates InChIKeys  
- Computes physchem descriptors  
- Writes molecule metadata (local + global)  
- Output: `cleaned.csv`

### 2. Generation

- Generates conformers (RDKit or CREST)  
- Output: `energies.json`

### 3. Pruning

- Selects top‑N conformers  
- Output: `energies.json`

### 4. Optimisation

- Runs ORCA, XTB, or gXTB  
- Checkpointed  
- Fully resumable  
- Output: `energies.json`

### 5. ORCA COSMO

- Writes ORCA input files  
- Runs TZVPD → fallback TZVP  
- Parses log/cpcm/cpcm_corr  
- Reconstructs `.orcacosmo` files  
- Output: `orcacosmo_summary.json`

### 6. Solubility

- Runs COSMO‑RS  
- Output: `solubility_results.json`

---

## 🔁 Resume & Continuation

### Resume (same request)

If a job is interrupted:

- `job_state.json` tracks pending items  
- `pipeline_state.json` tracks current stage  

Resume with:

```bash
python3 -m modules.main_resume
```

### Continue (new request)

To start a new pipeline using the output of a previous job:

```bash
python3 -m modules.main_continue
```

This creates a new Request with:

- New request ID  
- New job lineage  
- Stage 0 input = previous job’s canonical output  

---

## 🧬 Provenance & Reproducibility

Every metadata file includes:

- Git version (`main@abc1234`)  
- Cleaning timestamp  
- Request ID  
- Job ID  
- Source file  
- Pipeline version  

Every stage:

- Has a canonical input  
- Has a canonical output  
- Never guesses or auto‑detects  
- Is deterministic  

---

## 🛠 Adding a New Stage

See `docs/developer_guide.md` for full details.

In short:

1. Create `modules/stages/<name>_stage.py`
2. Subclass `BaseStage`
3. Implement:
   - `execute()`
   - `set_stage_output()`
   - `require_file()`
4. Add canonical output to `Job.STAGE_OUTPUTS`
5. Add stage to pipeline spec

---

## 🧪 Example

```python
pipeline_spec = [
    {"stage": "cleaning", "args": {"input_csv": "data.csv"}},
    {"stage": "generation", "args": {"engine": "rdkit"}},
    {"stage": "pruning", "args": {"n": 1}},
    {"stage": "optimisation", "args": {"engine": "xtb_opt_normal"}},
    {"stage": "orcacosmo", "args": {}},
    {"stage": "solubility", "args": {}},
]
```

---

## 📄 License

MIT or your preferred license.

---

## 👤 Authors

Luke Hetherington  
(openCOSMO‑RS Pipeline Architect)

---

