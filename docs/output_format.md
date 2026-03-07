# Output Format Reference

openCOSMO-RS Pipeline

---

## Overview

Each pipeline stage writes a canonical output file. The final results are collected in `final_outputs/` at the end of a run. This document describes the schema for each output format.

---

## `final_outputs/` layout

```
final_outputs/
    solubility_results.json          machine-readable results (last stage = solubility)
    solubility_results.csv           flat per-molecule/combo table
    solubility_human_summary.txt     aligned human-readable table
    molecule_metadata/               per-molecule JSON descriptors (from cleaning)
        <InChIKey>.json
        ...
    pipeline_warnings.txt            all warnings from the run
```

For runs ending at an earlier stage, the equivalent stage outputs appear instead of the solubility files.

---

## `solubility_results.json`

A JSON array. One object per unique molecule.

```json
[
  {
    "inchi_key":                           "UHOVQNZJYSORNB-UHFFFAOYSA-N",
    "mol_name":                            "benzene",
    "smiles":                              "c1ccccc1",
    "melting_temp":                        5.5,
    "melting_temp_source":                 "literature",
    "experimental_solubility_mol_frac":    null,
    "aqsol_predicted_solubility_mol_frac": null,
    "conformer_history": [
      {
        "lookup_id":            "UHOVQNZJYSORNB-UHFFFAOYSA-N_c001",
        "conf_num":             1,
        "energy":               -232.451,
        "provenance":           { ... },
        "optimisation_history": [ ... ],
        "orcacosmo_history":    { ... },
        "orcacosmo_path":       "jobs/.../outputs/orcacosmo_outputs/UHOVQNZJYSORNB-UHFFFAOYSA-N_c001.orcacosmo"
      }
    ],
    "solvent_results": [
      {
        "combo_id":             "water_only",
        "combo_label":          "Pure water",
        "solvents":             [{"name": "water", "mol_frac": 1.0}],
        "temperature_K":        298.15,
        "predicted_solubility": -3.21,
        "n_solute_confs":       3,
        "n_solvent_confs":      {"water": 1},
        "mixture_inputs_path":  "jobs/.../outputs/results/.../mixture_inputs.txt",
        "raw_output_path":      "jobs/.../outputs/results/.../raw_output.txt",
        "error":                null
      }
    ]
  }
]
```

### Top-level fields

| Field | Type | Description |
|-------|------|-------------|
| `inchi_key` | string | Standard InChIKey identifier |
| `mol_name` | string | Molecule name from input CSV |
| `smiles` | string | Canonical SMILES |
| `melting_temp` | float / `"liquid"` / `"N/A"` | Melting point used for saturation correction (K); `"liquid"` if below 25°C; `"N/A"` if unavailable |
| `melting_temp_source` | string | Source of the melting point value |
| `experimental_solubility_mol_frac` | float / null | Experimental reference value if available |
| `aqsol_predicted_solubility_mol_frac` | float / null | AqSol ML predicted value if available |
| `conformer_history` | list | List of conformer records used in the calculation |
| `solvent_results` | list | One entry per solvent combo |

### `solvent_results` entry

| Field | Type | Description |
|-------|------|-------------|
| `combo_id` | string | Unique combo identifier from the solvent list |
| `combo_label` | string | Human-readable label |
| `solvents` | list | Solvent components and mole fractions |
| `temperature_K` | float | Calculation temperature |
| `predicted_solubility` | float / null | Predicted log10 solubility in mol fraction. `null` if the calculation failed. |
| `n_solute_confs` | int | Number of solute conformers used |
| `n_solvent_confs` | object | Number of conformers used per solvent component |
| `mixture_inputs_path` | string | Path to the mixture inputs file (COSMO-RS input) |
| `raw_output_path` | string | Path to the raw COSMO-RS output |
| `error` | string / null | Error message if the combo failed; `null` on success |

---

## `solubility_results.csv`

Flat table — one row per (molecule, solvent combo) combination.

Key columns:

| Column | Description |
|--------|-------------|
| `inchi_key` | Molecule InChIKey |
| `mol_name` | Molecule name |
| `smiles` | Canonical SMILES |
| `combo_id` | Solvent combo ID |
| `combo_label` | Solvent combo label |
| `temperature_K` | Calculation temperature |
| `predicted_solubility` | log10 solubility in mol fraction |
| `melting_temp` | Melting point used |
| `melting_temp_source` | Source of melting point |
| `n_solute_confs` | Number of solute conformers |
| `error` | Error message or empty |

---

## `energies.json` (ConformerSet)

Used as the canonical output of the generation, pruning, and optimisation stages. A JSON array — one object per conformer.

```json
[
  {
    "lookup_id":            "UHOVQNZJYSORNB-UHFFFAOYSA-N_c001",
    "inchi_key":            "UHOVQNZJYSORNB-UHFFFAOYSA-N",
    "conf_num":             1,
    "smiles":               "c1ccccc1",
    "energy":               -232.451,
    "xyz_path":             "jobs/.../outputs/xyz/UHOVQNZJYSORNB-UHFFFAOYSA-N_c001.xyz",
    "provenance": {
      "source":             "rdkit",
      "seed":               42,
      "gfn":                null,
      "timestamp":          "2025-01-15T14:23:01"
    },
    "optimisation_history": [
      {
        "stage":            "optimisation",
        "engine":           "gxtb_opt_normal",
        "level":            "normal",
        "status":           "converged",
        "energy":           -232.451,
        "xyz_path":         "jobs/.../outputs/xyz/UHOVQNZJYSORNB-UHFFFAOYSA-N_c001_opt1.xyz",
        "log_path":         "jobs/.../outputs/log/UHOVQNZJYSORNB-UHFFFAOYSA-N_c001_opt1.log",
        "elapsed_seconds":  12.4,
        "timing":           {"wall_seconds": 12.4, "cpu_seconds": 24.8},
        "backend_meta":     {"run_status": "converged", "version": "6.7.1", "command_str": "gxtb ..."},
        "timestamp":        "2025-01-15T14:23:01"
      }
    ]
  }
]
```

### Conformer fields

| Field | Type | Description |
|-------|------|-------------|
| `lookup_id` | string | `<inchi_key>_c<NNN>` — unique conformer identifier |
| `inchi_key` | string | Parent molecule InChIKey |
| `conf_num` | int / null | Conformer number (1-indexed) |
| `smiles` | string | Canonical SMILES |
| `energy` | float | Energy in the units of the last successful computation (Hartree for DFT/GFN, kcal/mol for force fields) |
| `xyz_path` | string | Path to the current best XYZ geometry file |
| `provenance` | object | Generation source metadata |
| `optimisation_history` | list | Ordered history of optimisation attempts |

### `optimisation_history` entry

| Field | Type | Description |
|-------|------|-------------|
| `stage` | string | Stage name (`"optimisation"`) |
| `engine` | string | Engine key used (e.g. `"gxtb_opt_normal"`) |
| `level` | string | Optimisation level |
| `status` | string | `"converged"`, `"partial"`, or `"failed"` |
| `energy` | float / null | Energy at this stage; `null` if failed |
| `xyz_path` | string / null | Path to XYZ at this stage |
| `log_path` | string | Path to engine log file |
| `elapsed_seconds` | float | Wall time for this step |
| `timing` | object | Detailed timing breakdown |
| `backend_meta` | object | Tool version, command string, run status |
| `timestamp` | string | ISO-8601 timestamp |

---

## `orcacosmo_summary.json`

Output of the orcacosmo stage. A JSON array — one object per conformer.

```json
[
  {
    "item_id":              "UHOVQNZJYSORNB-UHFFFAOYSA-N_c001",
    "lookup_id":            "UHOVQNZJYSORNB-UHFFFAOYSA-N_c001",
    "inchi_key":            "UHOVQNZJYSORNB-UHFFFAOYSA-N",
    "conf_num":             1,
    "smiles":               "c1ccccc1",
    "energy":               -232.451,
    "provenance":           { ... },
    "optimisation_history": [ ... ],
    "orcacosmo_history": {
      "method_used":        "TZVP",
      "fallback_triggered": false,
      "elapsed_default":    142.3,
      "elapsed_fallback":   null,
      "elapsed_total":      142.3,
      "orca_version":       "6.1.1",
      "raw_output_paths": {
        "log":      "jobs/.../outputs/raw_outputs/UHOVQNZJYSORNB-UHFFFAOYSA-N_c001.log",
        "cpcm":     "jobs/.../outputs/raw_outputs/UHOVQNZJYSORNB-UHFFFAOYSA-N_c001.cpcm",
        "cpcm_corr":"jobs/.../outputs/raw_outputs/UHOVQNZJYSORNB-UHFFFAOYSA-N_c001.cpcm_corr"
      }
    },
    "orcacosmo_path":       "jobs/.../outputs/orcacosmo_outputs/UHOVQNZJYSORNB-UHFFFAOYSA-N_c001.orcacosmo",
    "status":               "ok"
  }
]
```

### `orcacosmo_history` fields

| Field | Type | Description |
|-------|------|-------------|
| `method_used` | string | Basis set that succeeded (`"TZVP"` or `"TZVPD"`) |
| `fallback_triggered` | bool | Whether the fallback basis set was used |
| `elapsed_default` | float | Wall time for the default basis attempt (seconds) |
| `elapsed_fallback` | float / null | Wall time for the fallback attempt; `null` if not triggered |
| `elapsed_total` | float | Total wall time (seconds) |
| `orca_version` | string | ORCA version string |
| `raw_output_paths` | object | Paths to the raw ORCA output files |

---

## `molecule_metadata/<InChIKey>.json`

Written by the cleaning stage. One file per unique molecule.

```json
{
  "inchi_key":             "UHOVQNZJYSORNB-UHFFFAOYSA-N",
  "mol_name":              "benzene",
  "smiles":                "c1ccccc1",
  "canonical_smiles":      "c1ccccc1",
  "inchi":                 "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H",
  "charge":                0,
  "charge_source":         "user_input",
  "multiplicity":          1,
  "multiplicity_source":   "user_input",
  "melting_temp_k":        278.65,
  "melting_temp_source":   "literature",
  "molecular_weight":      78.11,
  "molecular_formula":     "C6H6",
  "n_heavy_atoms":         6,
  "n_rotatable_bonds":     0,
  "n_rings":               1,
  "n_aromatic_rings":      1,
  "hbd":                   0,
  "hba":                   0,
  "tpsa":                  0.0,
  "logp":                  1.9,
  "functional_groups":     ["aromatic"],
  "source_file":           "molecules.csv",
  "row_index":             0
}
```

### Source fields for charge/multiplicity/melting point

| `_source` value | Meaning |
|----------------|---------|
| `"user_input"` | Value was explicitly provided in the input CSV |
| `"rdkit_computed"` | Value was computed from the SMILES by RDKit |
| `"default"` | A fallback default was applied (RDKit could not determine the value) |

---

## `pipeline_warnings.txt`

Plain-text file collecting all `[WARNING]`-level log entries from the entire run, grouped by stage and deduplicated.

```
============================================================
 PIPELINE WARNINGS SUMMARY
============================================================

[PIPELINE] (1 warning(s))
  → Resource budget reduced: cores_per_item=8 exceeds available budget; using cores_per_item=6

[CLEANING] (2 warning(s))
  → [HEADER] 'charge' column absent — will use RDKit or default 0 per molecule
  → [META] UHOVQNZJYSORNB-UHFFFAOYSA-N: charge not fully user-provided (3 molecule(s)) — user_input=0, rdkit_computed=3, default=0

[ORCACOSMO] (1 warning(s))
  → UHOVQNZJYSORNB-UHFFFAOYSA-N_c001: TZVPD failed — retrying with TZVP fallback
```
