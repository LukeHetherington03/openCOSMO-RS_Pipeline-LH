# Arguments and Parameters Guide

openCOSMO-RS Pipeline

---

## Contents

1. [Pipeline spec (`request.json`)](#1-pipeline-spec-requestjson)
2. [Parameter hierarchy](#2-parameter-hierarchy)
3. [Global parameters](#3-global-parameters)
4. [Cleaning stage parameters](#4-cleaning-stage-parameters)
5. [Generation stage parameters](#5-generation-stage-parameters)
6. [Pruning stage parameters](#6-pruning-stage-parameters)
7. [Optimisation stage parameters](#7-optimisation-stage-parameters)
8. [OrcaCOSMO stage parameters](#8-orcacosmo-stage-parameters)
9. [Solubility stage parameters](#9-solubility-stage-parameters)
10. [Resource parameters](#10-resource-parameters)

---

## 1. Pipeline spec (`request.json`)

A `request.json` defines a complete pipeline run. The top-level fields are:

```json
{
  "request_name": "my_run",
  "input_csv":    "/path/to/molecules.csv",
  "pipeline_sequence": ["cleaning", "generation", "optimisation", "orcacosmo", "solubility"],
  "stage_args": [
    {},
    {"engine": "rdkit", "num_confs": 20},
    {"engine": "gxtb_opt_normal"},
    {"default_basis": "TZVP"},
    {}
  ]
}
```

### Top-level fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `request_name` | string | No | Human-readable label. Defaults to timestamp. |
| `input_csv` | string or list | Yes* | Path(s) to input CSV file(s). Also accepts `input_csvs`, `input_folder`, `input_folders`. |
| `pipeline_sequence` | list of strings | Yes | Ordered list of stage names to run. |
| `stage_args` | list of objects | Yes | Per-stage parameter overrides. **Must be the same length as `pipeline_sequence`** and aligned by index. |
| `strict` | bool | No | Global strict mode — propagates to all stages unless overridden per-stage. |
| `resources` | object | No | Runtime resource override (see [Resource parameters](#10-resource-parameters)). |
| `tags` | list of strings | No | Metadata tags for organising requests. |
| `notes` | string | No | Free-text note attached to the request. |

*`input_csv` / `input_csvs` / `input_folder` / `input_folders` are passed through to the cleaning stage.

### `pipeline_sequence` and `stage_args`

`pipeline_sequence` is an ordered list of stage names. `stage_args` is a parallel list — `stage_args[0]` applies to `pipeline_sequence[0]`, `stage_args[1]` to `pipeline_sequence[1]`, and so on.

Each element of `stage_args` is a dict of parameter overrides for that stage. Empty dicts `{}` are valid and use all defaults.

**Valid stage names:** `cleaning`, `generation`, `pruning`, `optimisation`, `orcacosmo`, `solubility`

**Common sequences:**

```json
["cleaning", "generation", "optimisation", "orcacosmo", "solubility"]
["cleaning", "generation", "pruning", "optimisation", "orcacosmo", "solubility"]
["cleaning", "generation"]
["orcacosmo", "solubility"]
```

You can run partial pipelines by starting mid-sequence. The first stage must receive appropriate input (or `stage_input` must be provided manually).

---

## 2. Parameter hierarchy

Each stage resolves its parameters in order of precedence (highest wins):

```
1. stage_args entry in request.json          (runtime override)
2. request-level "strict" / "resources"      (cross-stage overrides)
3. stage defaults JSON  (e.g. optimisation_defaults.json)
4. BaseStage / hard-coded fallback
```

For example, if `optimisation_defaults.json` has `"engine": "gxtb_opt_normal"` and the request's `stage_args` entry has `"engine": "orca_opt_fast"`, the request-level value wins.

---

## 3. Global parameters

These can be set at the top level of `request.json` and apply across stages.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `strict` | bool | `false` | If `true`, any per-item failure aborts the entire pool. Each stage can also override this individually. |

---

## 4. Cleaning stage parameters

The cleaning stage accepts raw CSV(s), canonicalises SMILES, generates InChIKeys, resolves charge/multiplicity/melting point, and writes `cleaned.csv` and per-molecule metadata.

### Input source parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `input_csv` | string or list | Path to a single CSV or list of CSV paths |
| `input_csvs` | string or list | Alias for `input_csv` |
| `input_folder` | string | Path to a directory — all `.csv` files within are loaded |
| `input_folders` | string or list | List of directories |

At least one of the above is required.

### Behaviour parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `strict` | bool | `false` | If `true`, any row that fails validation causes the stage to abort |

### Input CSV columns

| Column | Required | Description |
|--------|----------|-------------|
| `smiles` | Yes | SMILES string |
| `mol_name` | Recommended | Human-readable molecule name |
| `charge` | Optional | Formal charge (integer). If absent, RDKit or default (0) used per molecule. |
| `multiplicity` | Optional | Spin multiplicity. If absent, RDKit (closed-shell) or default (1) used. |
| `melting_temp_c` | Optional | Melting point in Celsius. Used for solubility saturation correction. |
| `melting_temp_source` | Optional | Source label for the melting point (e.g. `"ChemBook"`, `"user_input"`). |

If `charge` or `multiplicity` columns are absent, a `[WARNING]` is emitted in the log with a per-file breakdown of how many molecules received RDKit-computed vs default values.

### Outputs

| File | Description |
|------|-------------|
| `outputs/cleaned.csv` | Unified, canonicalised dataset for downstream stages |
| `outputs/molecule_metadata/<inchi_key>.json` | Per-molecule physchem descriptors, charge, multiplicity, melting point |

---

## 5. Generation stage parameters

Generates 3D conformers for each unique molecule in `cleaned.csv`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `engine` | string | `"rdkit"` | Conformer generation backend: `"rdkit"`, `"crest"`, `"openbabel"` |
| `num_confs` | int | engine-dependent | Maximum number of conformers to generate |
| `seed` | int | `42` | Random seed for reproducibility |
| `gfn` | string | `"gfn2"` | GFN version for CREST (e.g. `"gfn2"`) |
| `strict` | bool | `false` | Abort on first per-molecule failure |
| `global_fail_threshold` | float | `0.8` | If the fraction of failures exceeds this, the stage raises an error |

### Engine comparison

| Engine | Speed | Quality | Parallelism | Requires |
|--------|-------|---------|-------------|---------|
| `rdkit` | Fast | Medium (ETKDGv3 + MMFF energy) | Serial | RDKit only |
| `crest` | Slow | High (iMTD-GC conformer search) | Parallel | CREST binary |
| `openbabel` | Medium | Medium (genetic algorithm) | Serial | `obabel` on PATH |

### Outputs

| File | Description |
|------|-------------|
| `outputs/energies.json` | ConformerSet JSON — one entry per conformer |
| `outputs/xyz/<inchi_key>_conf<NNN>.xyz` | 3D coordinates |
| `outputs/summary.csv` | Per-molecule summary |

---

## 6. Pruning stage parameters

Reduces the conformer set by energy window and/or RMSD deduplication. Takes `energies.json` from generation and outputs a filtered `energies.json`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `energy_window` | float | `5.0` | Keep conformers within this energy window of the minimum |
| `energy_window_units` | string | `"kcal"` | Units for `energy_window`: `"kcal"` or `"kj"` |
| `rmsd_threshold` | float | `0.5` | RMSD threshold (Angstrom) for deduplication. Conformers more similar than this are merged (lowest energy kept). |
| `keep_all` | bool | `false` | If `true`, skip pruning and pass all conformers through unchanged |
| `max_energy` | float or null | `null` | Absolute energy ceiling (in native units). Conformers above this are discarded. |
| `percentile` | float or null | `null` | Keep only conformers in the lowest N-th percentile by energy |
| `n` | int or null | `null` | Keep at most N conformers per molecule (lowest energy) |
| `n_high` | int or null | `null` | Additionally keep at most N conformers from the high-energy tail |
| `n_start` | int or null | `null` | Keep at most N conformers from the start of the sorted list |
| `strict` | bool | `false` | Abort on first failure |

---

## 7. Optimisation stage parameters

Geometry-optimises each conformer using the selected engine. Takes `energies.json` and outputs a new `energies.json` with updated geometries and energies.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `engine` | string | `"gxtb_opt_normal"` | Optimisation engine key — must match a key in `optimisation_engines.json` |
| `max_iter` | int | `250` | Maximum optimisation iterations (passed to engine) |
| `keep_scratch` | bool | `true` | Retain per-conformer scratch directories after completion |
| `strict` | bool | `false` | Abort pool on first failure |
| `global_fail_threshold` | float | `0.8` | Raise error if failure fraction exceeds this |

### Available engines

#### gXTB engines (recommended for large batches)

| Engine key | Level | Notes |
|-----------|-------|-------|
| `gxtb_opt_loose` | loose | Fastest, least converged |
| `gxtb_opt_normal` | normal | **Default** — good balance |
| `gxtb_opt_tight` | tight | More thorough |
| `gxtb_opt_vtight` | vtight | Most thorough |

#### XTB engines

| Engine key | Level |
|-----------|-------|
| `xtb_opt_loose` | loose |
| `xtb_opt_normal` | normal |
| `xtb_opt_tight` | tight |
| `xtb_opt_vtight` | vtight |

#### ORCA engines

| Engine key | Method | Basis | Notes |
|-----------|--------|-------|-------|
| `orca_opt_fast` | BP86 | def2-TZVP(-f) | DFT optimisation, no solvation |
| `orca_opt_final` | BP86 | def2-TZVP | DFT optimisation, no solvation |
| `orca_opt_cpcm_fast` | BP86 | def2-TZVP(-f) | DFT + CPCM solvation |
| `orca_opt_cpcm_final` | BP86 | def2-TZVP | DFT + CPCM solvation |
| `orca_sp_fast` | BP86 | def2-TZVP | Single-point |
| `orca_sp_final` | BP86 | def2-TZVPD | Single-point with CPCM |
| `orca_xtb2_alpb_opt` | XTB2 | — | ORCA-driven xTB2 with ALPB solvation |

#### Force field engines (fast pre-optimisation)

| Engine key | Force field |
|-----------|------------|
| `forcefield_mmff` | MMFF94 |
| `forcefield_uff` | UFF |

### Outputs

| File | Description |
|------|-------------|
| `outputs/energies.json` | Updated ConformerSet with optimisation_history entries |
| `outputs/xyz/<lookup_id>_opt<N>.xyz` | Optimised geometry |
| `outputs/log/<lookup_id>_opt<N>.log` | Engine log |
| `outputs/optimisation_summary.csv` | Per-conformer summary |

---

## 8. OrcaCOSMO stage parameters

Runs ORCA CPCM COSMO single-point calculations to produce `.orcacosmo` files for openCOSMO-RS.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `default_basis` | string | `"TZVP"` | Primary basis set: `"TZVP"` or `"TZVPD"` |
| `fallback_basis` | string or false | `false` | Basis set to try if the primary fails: `"TZVP"`, `"TZVPD"`, or `false` (no fallback — fail immediately) |
| `functional` | string | `"BP86"` | DFT functional |
| `solvation` | string | `"CPCM"` | Solvation model keyword |
| `maxcore` | int | `2000` | ORCA `%MaxCore` in MB per process |
| `strict` | bool | `false` | Abort pool on first failure |

### Outputs

| File | Description |
|------|-------------|
| `outputs/orcacosmo_summary.json` | Per-conformer results including COSMO file paths |
| `outputs/orcacosmo_summary.csv` | Tabular summary |
| `outputs/orcacosmo_outputs/<lookup_id>.orcacosmo` | COSMO file for openCOSMO-RS |
| `outputs/raw_outputs/<lookup_id>.{log,cpcm,cpcm_corr}` | Raw ORCA outputs |

---

## 9. Solubility stage parameters

Runs openCOSMO-RS solubility calculations for each molecule against a list of solvents.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `solvent_list` | string | built-in default | Name of a solvent list file from `CONSTANT_FILES/solvent_lists/`, or path to a custom JSON solvent list |
| `temperature` | float | `298.15` | Calculation temperature in Kelvin |
| `calculation_type` | string | `"all"` | Type of calculation. `"all"` runs all solvents; other values stage-specific |
| `initial_solute_x` | float | `0` | Initial mole fraction guess for the solute |
| `saturation.enabled` | bool | `true` | Apply solubility saturation correction using melting point |
| `saturation.Gfus_mode` | string | `"MyrdalYalkowsky"` | Gibbs fusion energy model: `"MyrdalYalkowsky"` or `"Hfus"` |
| `saturation.Hfus` | string or float | `"N/A"` | Enthalpy of fusion (kJ/mol), used when `Gfus_mode = "Hfus"` |
| `saturation.SORcf` | float | `1.0` | Correction factor for the fusion model |
| `parallel` | bool | `true` | Run solvent calculations in parallel |
| `n_workers` | int | `4` | Number of parallel workers for solubility calculations |
| `strict` | bool | `false` | Abort on first failure |

### Solvent list format

A solvent list JSON file is a list of solvent names that must match entries in `CONSTANT_FILES/solvents/`:

```json
["water", "ethanol", "acetone", "dmso"]
```

If no solvent list is configured, the pipeline falls back to the default list (a broad set of common solvents). If that is also unavailable, an emergency water-only fallback is used and a warning is emitted.

### Outputs

| File | Description |
|------|-------------|
| `outputs/solubility_results.json` | Full machine-readable results per molecule/solvent |
| `outputs/solubility_results.csv` | Flat per-molecule/solvent table |
| `outputs/solubility_human_summary.txt` | Aligned human-readable table |

---

## 10. Resource parameters

Resource allocation can be overridden at the request level or per-stage.

### Global resource override

Set in the top-level `resources` dict in `request.json`:

```json
{
  "resources": {
    "cpus": 16
  }
}
```

| Parameter | Description |
|-----------|-------------|
| `cpus` | Total CPU cores available. Overrides `max_cores` in `config/resource_allocation.json`. |

### Per-stage resource override

Set in the relevant `stage_args` entry:

```json
{
  "cores_per_item": 4,
  "n_workers": 5
}
```

| Parameter | Description |
|-----------|-------------|
| `cores_per_item` | CPU cores per item. Overrides the stage default from `generation_defaults.json` / `optimisation_defaults.json` etc. |
| `n_workers` | Number of parallel workers. If not set, computed as `floor(max_cores / cores_per_item)`. |

### Resource allocation hierarchy

```
1. stage_args["cores_per_item"]          per-stage runtime override
2. stage_args["n_workers"]               per-stage runtime override
3. resources["cpus"]                     global runtime ceiling
4. config/resource_allocation.json       machine configuration
5. psutil (CPU count)                    automatic fallback
```

The allocator always ensures `n_workers * cores_per_item <= max_cores`. If the configured `cores_per_item` would exceed the budget even for a single worker, it is reduced automatically and a warning is emitted.
