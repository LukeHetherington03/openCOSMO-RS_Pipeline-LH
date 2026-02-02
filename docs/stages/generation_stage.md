# Stage: Generation

## Purpose
Generate 3D conformers for each molecule in the cleaned dataset and produce a canonical conformer store (`energies.json`) plus XYZ geometries.

This stage:
- reads molecule‑level data from CleaningStage,
- canonicalises SMILES,
- generates 3D conformers,
- assigns conformer IDs,
- writes XYZ files,
- computes energies (if backend supports),
- records provenance for each conformer.

No molecule‑level metadata is modified here.

---

## Inputs

### Required
- `clean_dataset.csv` from CleaningStage  
  (auto‑detected from `inputs/` if not explicitly provided)

### Required columns
- `inchi_key`
- `smiles` (canonical or raw; canonicalised here)

### Optional columns (passed through)
- `mol_name`
- `mol_name_iupac`
- `charge`

---

## Parameters

| Parameter | Type | Default | Description |
|----------|------|---------|-------------|
| `engine` | string | `"rdkit"` | Conformer generation backend |
| `num_confs` | int | `5` | Number of conformers per molecule |
| `seed` | int | `42` | Random seed for reproducibility |
| `strict` | bool | `false` | If true, any failure aborts the stage |

---

## Processing

### 1. Load `clean_dataset.csv`
- Validate required columns.
- Canonicalise SMILES using RDKit.

### 2. For each molecule
- Generate conformers using the selected backend.
- Optimise geometry (UFF for RDKit).
- Compute energy if backend supports it.
- Assign:
  - `conf_num` (0‑indexed)
  - `lookup_id = "{inchi_key}_conf{conf_num:03d}"`

### 3. Write outputs
- `.xyz` files for each conformer
- `summary.csv` (conformer‑level table)
- `energies.json` (canonical conformer store)
- `job_state.json` (stage provenance)

---

## Outputs

### `xyz/`
One XYZ file per conformer:

