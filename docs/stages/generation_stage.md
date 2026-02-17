# Generation Stage  
openCOSMO‑RS Pipeline

The Generation Stage produces 3D conformers for each molecule in the cleaned dataset.  
It supports multiple backends (RDKit, CREST, OpenBabel), applies strict validation, and writes a canonical `energies.json` file along with per‑conformer XYZ files.

This document describes the Generation Stage in detail for both users and developers.

---

# 1. Purpose

The Generation Stage:

- Reads `cleaned.csv` from the Cleaning Stage  
- Validates molecules (SMILES + InChIKey)  
- Determines charge and spin from metadata  
- Generates 3D conformers using the selected backend  
- Computes generation‑level energies  
- Writes per‑conformer XYZ files  
- Produces a canonical `energies.json`  
- Produces a `summary.csv`  
- Records provenance for each conformer  

It is the foundation for pruning and optimisation.

---

# 2. Canonical Input & Output

### **Input**
The stage requires:

```
stage_input = cleaned.csv
```

This is enforced via:

```python
stage_input = self.get_stage_input()
self.require_file(stage_input, "cleaned dataset")
```

The file is copied into:

```
jobs/J-.../inputs/cleaned.csv
```

for reproducibility.

### **Output**
The stage declares its canonical output:

```
energies.json
```

and also writes:

```
summary.csv
xyz/<inchi_key>_confNNN.xyz
```

These appear under:

```
jobs/J-.../outputs/
```

---

# 3. Parameters

The Generation Stage reads the following parameters:

| Parameter | Type | Required | Description |
|----------|------|----------|-------------|
| `engine` | str | optional | `"rdkit"` (default), `"crest"`, `"openbabel"` |
| `n` | int | optional | Number of conformers per molecule |
| `seed` | int | optional | RNG seed (default: 42) |
| `threads` | int | optional | Thread count for RDKit/CREST |
| `charge` | int | optional | Override charge for all molecules |
| `spin` | int | optional | Override spin for all molecules |
| `gfn` | str | optional | `"gfn0"` or `"gfn2"` for CREST |
| `level` | str | optional | `"fast"` or `"normal"` (CREST) |

If `n` is not provided, a heuristic is used.

---

# 4. Strict Mode

Strict mode is enabled via:

```json
"generation": { "strict": true }
```

Strict mode enforces:

- Missing InChIKey → fail  
- Missing SMILES → fail  
- Invalid SMILES → fail  
- Backend failures → fail  

Without strict mode, these become warnings and the molecule is skipped.

---

# 5. Execution Flow

The Generation Stage follows this sequence:

1. **Load cleaned.csv**  
2. **Load metadata** (charge, spin, provenance)  
3. **Validate rows**  
4. **Determine conformer count**  
5. **Dispatch backend**  
6. **Write XYZ files**  
7. **Write energies.json**  
8. **Write summary.csv**  
9. **Write warning summary**  

---

# 6. Metadata Usage

Metadata is loaded from:

```
CONSTANT_FILES/molecule_metadata/<inchi_key>.json
```

Charge and spin are resolved with precedence:

1. Stage parameters (`charge`, `spin`)  
2. Metadata JSON  
3. Defaults: charge=0, spin=1  

This ensures consistent behaviour across backends.

---

# 7. Conformer Count Heuristic

If `n` is not provided, the stage uses the O’Boyle 2011 heuristic:

| Rotatable Bonds | Conformers |
|-----------------|------------|
| ≤ 7 | 50 |
| 8–12 | 200 |
| ≥ 13 | 300 |

The stage computes this per molecule and uses the **maximum** as a safe global value.

---

# 8. Backends

The stage supports three backends:

```
"rdkit"
"crest"
"openbabel"
```

Each backend is documented below.

---

# 9. RDKit Backend (Fully Functional)

The RDKit backend:

- Generates 3D conformers using ETKDGv3  
- Adds hydrogens  
- Uses MMFF94 (or UFF fallback) for energy evaluation  
- Writes XYZ files  
- Records provenance  

### Key features

- Fast  
- Pure Python  
- No external executables  
- Deterministic with seed  
- Good for initial conformer sets  

### Failure modes

- RDKit cannot parse SMILES  
- Embedding fails  
- Energy evaluation fails  

Strict mode converts warnings into errors.

---

# 10. CREST Backend (Fully Functional)

The CREST backend:

- Uses RDKit to generate a single seed conformer  
- Writes `input.xyz`  
- Runs CREST via subprocess  
- Parses `crest_conformers.xyz`  
- Parses `crest.energies`  
- Writes per‑conformer XYZ files  
- Records provenance including:
  - CREST version  
  - GFN level  
  - Charge/spin  

### GFN Level Resolution

Priority:

1. `gfn` parameter (`"gfn0"` or `"gfn2"`)  
2. `level` parameter (`"fast"` → gfn0, `"normal"` → gfn2`)  
3. Default: gfn2  

### Failure modes

- CREST executable missing  
- CREST returns non‑zero exit code  
- Missing CREST output files  
- No conformers parsed  

Strict mode converts warnings into errors.

---

# 11. OpenBabel Backend (NOT YET FUNCTIONAL)

The OpenBabel backend is present in the code but **not yet working**.

### Current state

- Backend is partially implemented  
- Requires `openbabel-wheel`  
- Code is truncated  
- Not tested  
- Not recommended for production  

### Documentation note

This backend is included for completeness but should be considered **experimental**.

---

# 12. Provenance

Each conformer record includes:

- Backend name  
- Backend version (RDKit/CREST/OpenBabel)  
- Seed  
- Threads  
- Charge/spin  
- Source row index  
- Timestamp  
- Forcefield (RDKit)  
- GFN level (CREST)  

This ensures full reproducibility.

---

# 13. Warning System

The stage tracks:

- Embedding failures  
- Optimisation failures  
- Zero conformers  
- CREST failures  
- OpenBabel failures  

At the end, a summary is printed.

---

# 14. Failure Modes

The stage fails if:

- Stage input missing  
- No valid molecules  
- Backend not found  
- Backend produces zero conformers  
- Strict mode violations  

Failures update:

```
pipeline_state.json
job_state.json
```

and halt the pipeline.

---

# 15. Minimal Example

### pipeline_spec entry

```python
{
    "stage": "generation",
    "args": {
        "engine": "rdkit",
        "n": 50,
        "seed": 42,
        "threads": 4
    }
}
```

### Running the stage

```
python3 main.py
```

---

# 16. Summary

The Generation Stage:

- Validates molecules  
- Loads metadata  
- Generates conformers  
- Computes energies  
- Writes XYZ files  
- Produces `energies.json`  
- Supports RDKit, CREST, and experimental OpenBabel  
- Provides deterministic, reproducible conformer sets  

It is the foundation for pruning and optimisation.

