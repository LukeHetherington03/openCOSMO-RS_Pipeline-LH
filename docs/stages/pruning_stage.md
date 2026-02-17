# Pruning Stage  
openCOSMO‑RS Pipeline

The Pruning Stage reduces the number of conformers per molecule by applying configurable pruning rules.  
It consumes the conformer set produced by the Generation or Optimisation stages and outputs a smaller, cleaner, and more physically meaningful conformer ensemble.

This document describes the Pruning Stage in detail for both users and developers.

---

# 1. Purpose

The Pruning Stage:

- Reads `energies.json` from Generation or Optimisation  
- Groups conformers by molecule (InChIKey)  
- Removes invalid or missing‑energy conformers  
- Applies user‑defined pruning rules  
- Writes a canonical pruned `energies.json`  
- Writes a `pruning_summary.csv`  
- Tracks warnings and strict‑mode failures  

It is the final step before geometry optimisation (if pruning follows generation) or before ORCA COSMO (if pruning follows optimisation).

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
pruning_summary.csv
```

These appear under:

```
jobs/J-.../outputs/
```

---

# 3. Parameters

The Pruning Stage supports a rich set of pruning parameters.  
All are optional — the user may specify any combination.

| Parameter | Type | Description |
|----------|------|-------------|
| `keep_all` | bool | Skip pruning entirely |
| `rmsd_threshold` | float | RMSD‑based pruning (Å) |
| `energy_window` | float | Keep conformers within ΔE of minimum |
| `max_energy` | float | Keep conformers with E ≤ max_energy |
| `percentile` | float | Keep conformers below percentile cutoff |
| `n` | int | Keep lowest‑energy N conformers |
| `n_high` | int | Keep highest‑energy N conformers |
| `n_start` | int | Slice start index (paired with `n`) |

Parameters are applied **in the order defined in the registry**, allowing complex pruning pipelines.

---

# 4. Strict Mode

Strict mode is enabled via:

```json
"pruning": { "strict": true }
```

Strict mode enforces:

- Missing energies → fail  
- All conformers pruned → fail  
- RMSD parsing failures → fail  

Without strict mode, these become warnings and the molecule is skipped.

---

# 5. Execution Flow

The Pruning Stage follows this sequence:

1. **Load energies.json**  
2. **Group conformers by molecule**  
3. **Track items for resume**  
4. **For each molecule:**  
   - Remove missing‑energy conformers  
   - Apply pruning rules  
   - Sort survivors  
   - Record summary row  
5. **Write pruned energies.json**  
6. **Write pruning_summary.csv**  
7. **Write warning summary**  

---

# 6. Dynamic Pruning Registry

The stage uses a dynamic registry:

```python
PRUNING_METHODS = {
    "rmsd_threshold": "_prune_rmsd",
    "energy_window": "_prune_energy_window",
    "max_energy": "_prune_max_energy",
    "percentile": "_prune_percentile",
    "n": "_prune_keep_lowest_n",
    "n_high": "_prune_keep_highest_n",
    "n_start": "_prune_slice",
}
```

This allows:

- Flexible pruning pipelines  
- Arbitrary combinations of methods  
- Clear provenance in `pruning_summary.csv`  

Each method is documented below.

---

# 7. Missing‑Energy Filtering (Always Applied)

Before any pruning rules, the stage removes conformers with:

- `energy = None`
- `energy = NaN`

These are logged and counted.

If **all** conformers are missing energy:

- Warning is recorded  
- Strict mode → fail  

---

# 8. Pruning Methods

Each method receives the current list of survivors and returns a new list.

### 8.1 RMSD Threshold Pruning

```
rmsd_threshold = float
```

Removes conformers whose RMSD to any already‑kept conformer is below the threshold.

- Uses RDKit’s `GetBestRMS`  
- Loads XYZ files via `Chem.MolFromXYZFile`  
- Keeps lowest‑energy conformers first  

### 8.2 Energy Window Pruning

```
energy_window = float  # default units: kcal
```

Keeps conformers with:

```
E - E_min <= energy_window
```

Energy units may be:

- hartree  
- kcal  
- kJ  

### 8.3 Max Energy Pruning

```
max_energy = float
```

Keeps conformers with:

```
E <= max_energy
```

### 8.4 Percentile Pruning

```
percentile = float
```

Keeps conformers with:

```
E <= percentile_cutoff
```

where cutoff is computed from the energy distribution.

### 8.5 Keep Lowest N

```
n = int
```

Keeps the lowest‑energy N conformers.

### 8.6 Keep Highest N

```
n_high = int
```

Keeps the highest‑energy N conformers.

### 8.7 Slice Pruning

```
n_start = int
n = int
```

Keeps a slice of the sorted conformer list:

```
survivors = conformers[n_start : n_start + n]
```

Supports negative indexing.

---

# 9. Summary Output

The stage writes:

```
pruning_summary.csv
```

with columns:

- Molecule  
- Total Conformers  
- With Valid Energy  
- Removed (Missing Energy)  
- Final Count  
- Pruning Steps Applied  
- RMSD Threshold (Å)  
- Energy Window  
- Max Energy  
- Percentile Cutoff  
- Keep Lowest N  
- Keep Highest N  
- Slice Start  

This provides a complete audit trail of pruning decisions.

---

# 10. Warning System

The stage tracks:

- Molecules with all conformers missing energy  
- Molecules with all conformers pruned  

At the end, a summary is printed.

---

# 11. Failure Modes

The stage fails if:

- Stage input missing  
- All conformers missing energy (strict mode)  
- All conformers pruned (strict mode)  
- RMSD parsing fails (strict mode)  
- Invalid pruning parameters  

Failures update:

```
pipeline_state.json
job_state.json
```

and halt the pipeline.

---

# 12. Minimal Example

### pipeline_spec entry

```python
{
    "stage": "pruning",
    "args": {
        "rmsd_threshold": 0.5,
        "energy_window": 3.0,
        "n": 10
    }
}
```

### Running the stage

```
python3 main.py
```

---

# 13. Summary

The Pruning Stage:

- Removes invalid conformers  
- Applies flexible pruning rules  
- Produces a clean, compact conformer set  
- Writes a canonical `energies.json`  
- Provides full provenance and summary reporting  
- Supports strict mode for reproducibility  

It is the final step before geometry optimisation or ORCA COSMO.

