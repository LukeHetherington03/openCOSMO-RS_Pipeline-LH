# Stage: Pruning

## Purpose
Reduce the number of conformers per molecule according to energy‑ and geometry‑based criteria, while preserving conformer identity and schema.

This stage:
- removes invalid or missing energies,
- applies pruning operations in a deterministic order,
- preserves all conformer metadata,
- writes a pruned `energies.json` and a molecule‑level pruning summary.

No geometry optimisation is performed here.

---

## Inputs

### Primary input
`energies.json` from GenerationStage or a previous PruningStage.

If `energies_file` is not provided, the stage auto‑detects:

- `inputs/energies.json`
- `outputs/energies.json`

### Expected schema (per conformer)
Each entry must match the canonical `ConformerRecord` schema:

- `lookup_id` — composite key `"{inchi_key}_conf{conf_num:03d}"`
- `inchi_key` — molecule identity
- `conf_num` — integer conformer index
- `xyz_path` — path to XYZ file
- `energy` — float or null
- `smiles` — canonical SMILES
- `provenance` — dict (backend, seed, timestamps, etc.)

---

## Arguments (`args`)

All pruning arguments are optional.  
If none are provided, only missing‑energy removal is applied.

### Energy‑based pruning
- `energy_window: float`  
  Keep conformers within ΔE of the minimum energy.

- `max_energy: float`  
  Drop conformers with `energy > max_energy`.

- `percentile: float`  
  Keep conformers with energy ≤ the given percentile cutoff.

- `n: int`  
  Keep the lowest N energies.  
  Also used with `n_start` for slicing.

- `n_high: int`  
  Keep the highest N energies.

### Geometry‑based pruning
- `rmsd_threshold: float`  
  RMSD clustering (currently a placeholder; logs and keeps all).

### Index‑based slicing
- `n_start: int`  
  Slice start index (supports negative indices).  
  Used together with `n`.

Examples:
- `n_start = 10, n = 5` → keep indices 10–14  
- `n_start = -1, n = 5` → keep last 5 conformers  

### Strict mode
`config["pruning"]["strict"] = True`

If a molecule has **all conformers missing energy**, the stage raises and fails.

---

## Processing

### 1. Group conformers by molecule (`inchi_key`)

### 2. Remove missing energies
Drop conformers where:
- `energy is None`, or
- `energy` is `NaN`.

If all conformers are removed:
- log an error,
- increment `molecules_all_missing_energy`,
- in strict mode → fail,
- in non‑strict mode → skip pruning for that molecule.

### 3. Apply pruning operations in deterministic order
If the corresponding arguments are present:

1. `rmsd_threshold` → RMSD pruning (placeholder)
2. `energy_window` → keep within ΔE of minimum
3. `max_energy` → keep `energy <= max_energy`
4. `percentile` → keep `energy <= percentile cutoff`
5. `n` → keep lowest N energies
6. `n_high` → keep highest N energies
7. `n_start + n` → slice by index (supports negative `n_start`)

### 4. Build pruning summary row
For each molecule:

- `inchi_key`
- `total_conformers`
- `valid_energy_conformers`
- `removed_missing_energy`
- `kept_after_pruning`
- all pruning args used

---

## Outputs

### `energies.json` (pruned)
Same schema as input, but with fewer conformers.  
All conformer fields (`lookup_id`, `inchi_key`, `conf_num`, `xyz_path`, `smiles`, `provenance`) are preserved.

### `pruning_summary.csv`
One row per molecule:

- `inchi_key`
- `total_conformers`
- `valid_energy_conformers`
- `removed_missing_energy`
- `kept_after_pruning`
- `rmsd_threshold`
- `energy_window`
- `max_energy`
- `percentile`
- `n`
- `n_high`
- `n_start`

---

## Guarantees

- Conformer identity is preserved:
  - `lookup_id`, `inchi_key`, `conf_num`, `xyz_path`, `smiles`, `provenance` remain unchanged.
- Only the **set** of conformers changes, not their content.
- Output `energies.json` remains fully compatible with OptimisationStage.
- Missing energies are always removed or explicitly reported.
- Pruning is deterministic and reproducible.

---

## Warnings

- Molecules with all conformers missing energy are counted and summarised.
- RMSD pruning is currently a placeholder and logs a warning when used.
- All warnings are summarised at the end of the stage.

