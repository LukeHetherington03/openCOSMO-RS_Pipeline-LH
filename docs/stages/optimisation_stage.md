# Optimisation Stage  
openCOSMO‑RS Pipeline

The Optimisation Stage performs geometry optimisation on each conformer using one of several supported quantum‑chemical or semi‑empirical engines.  
It consumes the conformer set produced by the Generation or Pruning stages and produces a fully optimised `energies.json` with detailed provenance and timing information.

This document describes the Optimisation Stage in detail for both users and developers.

---

# 1. Purpose

The Optimisation Stage:

- Reads `energies.json` from Generation or Pruning  
- Loads engine defaults and engine registry from config  
- Performs geometry optimisation using ORCA, gXTB, XTB, or a forcefield placeholder  
- Tracks convergence, energies, and timing  
- Writes per‑iteration XYZ files and logs  
- Supports resume and checkpointing  
- Produces a canonical optimised `energies.json`  
- Writes an optional human‑readable summary  

It is the final geometry refinement step before ORCA COSMO or solubility prediction.

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

for reproducibility.

### **Output**
The stage declares its canonical output:

```
energies.json
```

and also writes:

```
xyz/<lookup_id>_optN.xyz
log/<lookup_id>_optN.log
checkpoints/<lookup_id>.json
optimisation_summary.csv  (optional)
```

These appear under:

```
jobs/J-.../outputs/
```

---

# 3. Parameters

The Optimisation Stage merges:

1. Engine defaults (`optimisation.defaults`)  
2. Engine registry (`optimisation.engines`)  
3. User overrides (`pipeline_spec`)  

Supported parameters include:

| Parameter | Type | Description |
|----------|------|-------------|
| `engine` | str | Engine name from registry |
| `level` | str | `"loose"`, `"normal"`, `"tight"`, `"vtight"` |
| `max_iter` | int | Max optimisation iterations |
| `gfn` | int | GFN level for XTB/gXTB |
| `keep_scratch` | bool | Keep scratch directories |
| `resume` | bool | Resume from checkpoints |
| `global_fail_threshold` | float | Abort if failure ratio exceeds threshold |

---

# 4. Strict Mode

Strict mode is enabled via:

```json
"optimisation": { "strict": true }
```

Strict mode enforces:

- Missing XYZ → fail  
- Backend failures → fail  
- All conformers unusable → fail  

Without strict mode, these become warnings and the conformer is skipped.

---

# 5. Execution Flow

The Optimisation Stage follows this sequence:

1. **Load energies.json**  
2. **Load engine defaults and registry**  
3. **Merge parameters**  
4. **Determine engine family**  
5. **Resume or start fresh**  
6. **Prepare XYZ inputs**  
7. **Optimise each conformer**  
8. **Write checkpoints**  
9. **Merge checkpoints into final energies.json**  
10. **Write summary**  

---

# 6. Resume & Checkpointing

The stage supports robust resume:

- Each conformer has a unique `lookup_id`  
- Each optimisation iteration writes a checkpoint:  
  ```
  outputs/checkpoints/<lookup_id>.json
  ```
- If the stage is interrupted, resume picks up only pending conformers  
- Checkpoints are merged at the end into the final `energies.json`  

This ensures deterministic recovery and safe HPC operation.

---

# 7. Input Preparation

The stage copies all XYZ files into:

```
inputs/xyz/
```

Missing XYZ files are logged and skipped.  
If no valid conformers remain, the stage fails.

---

# 8. Engine Selection

The engine is selected from:

```
optimisation.engines
```

Each engine entry defines:

- `family` (orca, gxtb, xtb, forcefield)  
- Default optimisation level  
- Method/basis (ORCA)  
- GFN level (XTB/gXTB)  
- CPCM/ALPB settings (ORCA)  

The stage then dispatches to the appropriate backend.

---

# 9. Backends

The stage supports four backend families.

---

## 9.1 ORCA Backend (Fully Functional)

The ORCA backend:

- Writes a full `.inp` file  
- Supports:
  - LooseOpt / Opt / TightOpt / VeryTightOpt  
  - CPCM (default or custom radii)  
  - ALPB solvation  
  - Arbitrary method/basis from engine registry  
- Runs ORCA via subprocess  
- Copies `<lookup>.xyz` to output directory  
- Parses logs using `ORCALogParser`  

### Failure modes

- ORCA executable missing  
- ORCA returns non‑zero exit code  
- Final XYZ missing  
- Parser fails  

Strict mode converts warnings into errors.

---

## 9.2 gXTB Backend (Fully Functional)

The gXTB backend:

- Uses XTB as the driver  
- Uses gXTB for gradients  
- Builds a driver string:  
  ```
  gxtb -grad -c xtbdriver.xyz
  ```
- Supports loose/normal/tight/vtight optimisation  
- Uses 80% of available CPU cores  
- Writes `xtbopt.xyz` as output  

### Failure modes

- Missing XTB or gXTB executables  
- Missing `xtbopt.xyz`  
- Non‑zero exit code  

---

## 9.3 Classic XTB Backend (Fully Functional)

The classic XTB backend:

- Runs XTB directly  
- Supports GFN levels  
- Iteration count depends on optimisation level  
- Writes `xtbopt.xyz`  

### Failure modes

- Missing XTB executable  
- Missing output XYZ  
- Non‑zero exit code  

---

## 9.4 Forcefield Backend (Placeholder)

A simple placeholder backend:

- Copies input XYZ to output  
- Writes a minimal log  
- Does not perform optimisation  

Useful for debugging or pipeline testing.

---

# 10. Log Parsing

The stage uses:

```
ORCALogParser
XTBLogParser
GxTBLogParser
```

Each parser extracts:

- Convergence status  
- Final energy  
- Additional metadata (if available)  

If parsing fails, the stage logs a warning and marks the conformer as partial.

---

# 11. Status Determination

A conformer is considered:

- **converged** → parser says converged AND XYZ exists  
- **partial** → parser says converged but XYZ missing  
- **failed** → parser says not converged  

This status is recorded in the optimisation history.

---

# 12. Energy Validation

Energy is considered invalid if:

- `None`  
- `NaN`  
- `abs(energy) > 1e6`  

Invalid energies are replaced with `None`.

---

# 13. Optimisation History

Each conformer accumulates a history:

```
{
  "stage": "optimisation",
  "engine": "...",
  "level": "...",
  "status": "...",
  "energy": ...,
  "xyz_path": "...",
  "log_path": "...",
  "elapsed_seconds": ...,
  "timing": {...},
  "backend_meta": {...},
  "timestamp": "..."
}
```

This provides full provenance for downstream analysis.

---

# 14. Global Failure Threshold

If the fraction of failed conformers exceeds:

```
global_fail_threshold
```

(default: 0.8)

the stage aborts early.

This prevents wasting compute on hopeless cases