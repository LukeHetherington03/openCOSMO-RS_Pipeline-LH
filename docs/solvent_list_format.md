# Solvent List Format

openCOSMO-RS Pipeline

---

## Overview

The solubility stage calculates predicted solubility against one or more solvent environments (called **combos**). Each combo can be a pure solvent or a mixture of solvents at specified mole fractions.

Solvent combinations are defined in **solvent list JSON files** stored in `CONSTANT_FILES/solvent_lists/`.

---

## Solvent list JSON format

A solvent list is a JSON object with two fields:

```json
{
  "description": "Human-readable description of this list",
  "combos": {
    "<combo_id>": {
      "label": "Human-readable combo label",
      "solvents": [
        {"name": "<solvent_name>", "mol_frac": <fraction>},
        ...
      ]
    },
    ...
  }
}
```

### Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `description` | string | No | Description of the list. Stored in outputs for traceability. |
| `combos` | object | Yes | Dict of combo definitions, keyed by a unique combo ID. |

### Combo fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `label` | string | Yes | Human-readable label used in output tables and summaries. |
| `solvents` | list | Yes | List of solvent components. For a pure solvent, one entry with `mol_frac: 1.0`. For a mixture, multiple entries summing to 1.0. |

### Solvent entry fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `name` | string | Yes | Solvent name. Must match a directory name under `CONSTANT_FILES/solvents/`. |
| `mol_frac` | float | Yes | Mole fraction of this component in the solvent mixture. All `mol_frac` values in one combo must sum to 1.0. |

---

## Examples

### Pure solvent

```json
{
  "description": "Water only",
  "combos": {
    "water_only": {
      "label": "Pure water",
      "solvents": [
        {"name": "water", "mol_frac": 1.0}
      ]
    }
  }
}
```

### Multiple pure solvents

```json
{
  "description": "Common pharmaceutical solvents",
  "combos": {
    "water_only": {
      "label": "Pure water",
      "solvents": [{"name": "water", "mol_frac": 1.0}]
    },
    "glycerol_only": {
      "label": "Pure glycerol",
      "solvents": [{"name": "glycerol", "mol_frac": 1.0}]
    },
    "methanol_only": {
      "label": "Pure methanol",
      "solvents": [{"name": "methanol", "mol_frac": 1.0}]
    }
  }
}
```

### Binary mixtures

```json
{
  "description": "Water/glycerol mixtures",
  "combos": {
    "water_only": {
      "label": "Pure water",
      "solvents": [{"name": "water", "mol_frac": 1.0}]
    },
    "glycerol_only": {
      "label": "Pure glycerol",
      "solvents": [{"name": "glycerol", "mol_frac": 1.0}]
    },
    "50_50_water_glycerol": {
      "label": "50:50 water/glycerol",
      "solvents": [
        {"name": "water",    "mol_frac": 0.5},
        {"name": "glycerol", "mol_frac": 0.5}
      ]
    },
    "20_80_water_glycerol": {
      "label": "20:80 water/glycerol",
      "solvents": [
        {"name": "water",    "mol_frac": 0.2},
        {"name": "glycerol", "mol_frac": 0.8}
      ]
    }
  }
}
```

---

## Available solvents

The following solvents have pre-computed COSMO files in `CONSTANT_FILES/solvents/`:

| Solvent name | Notes |
|-------------|-------|
| `water` | |
| `methanol` | |
| `glycerol` | |
| `1,2-propanediol` | |
| `choline` | |
| `chloride` | |
| `betaine` | |
| `b-alanine` | |
| `citric_acid` | 2 conformers |
| `formic_acid` | |
| `fructose` | |
| `glucose` | |
| `l-tartaric_acid` | |
| `lactic_acid` | |
| `levulinic_acid` | 2 conformers |
| `malonic_acid` | |
| `menthol` | |
| `oxalic_acid` | |
| `proline` | |
| `sucrose` | |
| `thymol` | |
| `urea` | |
| `bitartrate` | |
| `2-phenyl_acetic_acid` | |
| `3-phenyl_propionic_acid` | |

The `name` field in the solvent list must exactly match the directory name under `CONSTANT_FILES/solvents/`.

---

## Adding a new solvent

To add a new solvent:

1. **Obtain the COSMO file** — run an ORCA CPCM calculation for the solvent molecule (the same calculation the `orcacosmo` stage runs for solutes). This produces a `.orcacosmo` file.

2. **Create the solvent directory:**
   ```
   CONSTANT_FILES/solvents/<solvent_name>/
   ```

3. **Place the COSMO file(s)** with the naming convention:
   ```
   <solvent_name>_c000.orcacosmo
   <solvent_name>_c001.orcacosmo   (if multiple conformers)
   ...
   ```

4. **Use the solvent name** in your solvent list JSON.

Solvents with multiple conformers are automatically handled — all conformer files found in the directory are used.

---

## Using a solvent list in a request

Specify the list by name (stem or filename) in the solubility stage args:

```json
{
  "pipeline_sequence": ["...", "solubility"],
  "stage_args": [
    ...,
    {"solvent_list": "my_solvents"}
  ]
}
```

The pipeline looks for `my_solvents.json` (or `my_solvents`) in `CONSTANT_FILES/solvent_lists/`. You can also supply an absolute path to a file anywhere on disk.

If `solvent_list` is omitted, the pipeline uses `default_list.json` (pure water). If that file is also missing, an emergency water-only fallback is used and a warning is emitted.

---

## Fallback behaviour

| Condition | Behaviour |
|-----------|-----------|
| `solvent_list` not specified | Uses `CONSTANT_FILES/solvent_lists/default_list.json` |
| `default_list.json` not found | Uses built-in water-only emergency combo; emits `[WARNING]` |
| Named solvent not found in `CONSTANT_FILES/solvents/` | That combo is skipped; emits `[WARNING]` |
| Combo `mol_frac` values do not sum to 1.0 | Warning emitted; calculation proceeds (values are used as-is) |
