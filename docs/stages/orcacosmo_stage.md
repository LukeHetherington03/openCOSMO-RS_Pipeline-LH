# ORCACOSMO Stage  
openCOSMO‑RS Pipeline

The ORCACOSMO Stage performs ORCA CPCM single‑point calculations on optimised geometries to generate COSMO surface data required for COSMO‑RS solubility prediction.  
It is the most computationally intensive stage and includes robust fallback logic, strict validation, and detailed provenance tracking.

This document describes the ORCACOSMO Stage in detail for both users and developers.

---

# 1. Purpose

The ORCACOSMO Stage:

- Reads `energies.json` from the Optimisation Stage  
- Loads CPCM radii and chemistry constants  
- Runs ORCA CPCM single‑point calculations (TZVPD → fallback TZVP)  
- Parses `.log`, `.cpcm`, and `.cpcm_corr` files  
- Reconstructs `.orcacosmo` files using the orchestrator  
- Writes a canonical `orcacosmo_summary.json`  
- Writes raw outputs, parsed outputs, and a summary CSV  
- Tracks fallback usage, failures, and missing XYZs  

It is the final quantum‑chemical stage before solubility prediction.

---

# 2. Canonical Input & Output

### **Input**
The stage requires:

```
stage_input = energies.json
```

This file is copied into:

```
jobs/J-.../inputs/energies.json
```

All XYZ files referenced in the entries are also copied into:

```
jobs/J-.../inputs/<lookup_id>.xyz
```

### **Output**
The stage declares its canonical output:

```
orcacosmo_summary.json
```

and also writes:

```
orcacosmo_summary.csv
item_to_lookup_mapping.json
raw_outputs/<lookup_id>.{log,cpcm,cpcm_corr}
parsed_outputs/<lookup_id>.{log,cpcm,cpcm_corr}.json
orcacosmo_outputs/<lookup_id>.orcacosmo
```

These appear under:

```
jobs/J-.../outputs/
```

---

# 3. Parameters & Configuration

The stage reads configuration from:

```
config/paths.json
```

### Required configuration keys

| Key | Description |
|-----|-------------|
| `orca.executable` | Path to ORCA binary |
| `orca.home` | ORCA installation directory |
| `constant_files.chemistry_dir` | CPCM radii JSON |
| `constant_files.metadata_dir` | Molecule metadata |
| `enable_tzvp_fallback` | Whether to allow fallback |
| `cosmo.strict` | Strict mode toggle |

### Strict Mode

Strict mode enforces:

- Missing XYZ → fail  
- Missing CPCM files → fail  
- Both TZVPD and TZVP failing → fail  

Without strict mode, these become warnings and the molecule is skipped.

---

# 4. Execution Flow

The ORCACOSMO Stage follows this sequence:

1. **Load configuration**  
2. **Load CPCM radii**  
3. **Prepare output directories**  
4. **Load entries + copy XYZs**  
5. **Discover lookup IDs**  
6. **For each lookup ID:**  
   - Assign item number  
   - Prepare workdir  
   - Write ORCA input (TZVPD)  
   - Run ORCA  
   - Validate CPCM outputs  
   - If invalid → fallback to TZVP  
   - Copy raw outputs  
   - Parse `.log`, `.cpcm`, `.cpcm_corr`  
   - Write parsed JSONs  
   - Build bundle  
   - Reconstruct `.orcacosmo` file  
   - Append summary entry  
7. **Write final summary JSON + CSV**  
8. **Write item‑to‑lookup mapping**  
9. **Print warning summary**  

---

# 5. CPCM Radii

The stage loads:

```
CONSTANT_FILES/chemistry/cpcm_radii.json
```

This file contains:

- Element‑specific radii  
- `cut_area` parameter  

These values are injected into ORCA input files.

---

# 6. ORCA Input Generation

The stage writes two possible ORCA input files:

### **TZVPD (primary method)**

```
! CPCM BP86 def2-TZVPD SP
```

### **TZVP (fallback method)**

```
! CPCM BP86 def2-TZVP SP
```

Both include:

- `%MaxCore 2000`  
- `%base` directive  
- `%cpcm` block with radii and cut_area  
- `%elprop` block  
- `* xyzfile 0 1 <lookup_id>.xyz`  

These files are stored in:

```
outputs/orca_inputs/<lookup_id>_tzvpd.inp
outputs/orca_inputs/<lookup_id>_tzvp.inp
```

---

# 7. ORCA Execution

ORCA is executed via:

```
orca <input.inp>
```

with environment:

- `OMP_NUM_THREADS = max_procs`  
- `LD_LIBRARY_PATH` set to ORCA directory  
- MPI variables removed to avoid conflicts  

Logs are written to:

```
tmp_exec/itemNNN/itemNNN.log
```

---

# 8. Fallback Logic (TZVPD → TZVP)

The stage first attempts:

```
TZVPD
```

If CPCM validation fails:

- If fallback disabled → fail  
- If fallback enabled → run TZVP  

Validation checks:

- `.cpcm` exists and >100 bytes  
- `.cpcm_corr` exists and >100 bytes  
- `.log` exists and >1000 bytes  

If both TZVPD and TZVP fail → stage fails.

---

# 9. Raw Output Copying

After ORCA finishes, raw outputs are copied to:

```
raw_outputs/<lookup_id>.log
raw_outputs/<lookup_id>.cpcm
raw_outputs/<lookup_id>.cpcm_corr
```

Missing files are logged.

---

# 10. Parsing

The stage uses three parsers:

- `OrcaLogParser`
- `OrcaCpcmParser`
- `OrcaCpcmCorrParser`

Parsed JSONs are written to:

```
parsed_outputs/<lookup_id>.log.json
parsed_outputs/<lookup_id>.cpcm.json
parsed_outputs/<lookup_id>.cpcm_corr.json
```

---

# 11. Bundle Construction

Each lookup ID produces a bundle:

```
{
  "meta": {...},
  "paths": {...},
  "optimisation_entry": {...}
}
```

This bundle is written to:

```
parsed_outputs/<lookup_id>_bundle.json
```

It includes:

- Lookup ID  
- InChIKey  
- SMILES  
- Method used (TZVPD/TZVP)  
- Fallback flag  
- ORCA version  
- CPCM radii source  
- Paths to parsed JSONs  
- Original optimisation entry  

---

# 12. COSMO Reconstruction

The stage uses:

```
OrcaCosmoOrchestrator(bundle)
```

to reconstruct the `.orcacosmo` file.

Output is written to:

```
orcacosmo_outputs/<lookup_id>.orcacosmo
```

This file is consumed by the Solubility Stage.

---

# 13. Summary Outputs

The stage writes:

### **Canonical output**
```
orcacosmo_summary.json
```

### **CSV summary**
```
orcacosmo_summary.csv
```

Columns include:

- lookup_id  
- inchi_key  
- smiles  
- energy  
- method_used  
- fallback_triggered  
- elapsed_seconds  
- orcacosmo_path  
- item_number  

### **Item mapping**
```
item_to_lookup_mapping.json
```

---

# 14. Warning System

The stage tracks:

- Successful items  
- Failed items  
- Fallback usage  
- Missing XYZ files  

A summary is printed at the end.

---

# 15. Failure Modes

The stage fails if:

- ORCA executable missing  
- CPCM radii missing  
- XYZ missing (strict mode)  
- Both TZVPD and TZVP fail  
- Parser errors (strict mode)  
- Missing CPCM files (strict mode)  

Failures update:

```
pipeline_state.json
job_state.json
```

and halt the pipeline.

---

# 16. Minimal Example

### pipeline_spec entry

```python
{
    "stage": "orcacosmo",
    "args": {
        "resources": {"cpus": 8}
    }
}
```

### Running the stage

```
python3 main.py
```

---

# 17. Summary

The ORCACOSMO Stage:

- Runs ORCA CPCM single‑point calculations  
- Supports TZVPD with automatic TZVP fallback  
- Validates CPCM outputs  
- Parses log and CPCM files  
- Reconstructs `.orcacosmo` files  
- Provides full provenance and summary reporting  
- Supports strict mode for reproducibility  

It is the final quantum‑chemical stage before solubility prediction.

