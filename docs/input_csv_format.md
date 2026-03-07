# Input CSV Format

openCOSMO-RS Pipeline

---

## Overview

The pipeline accepts one or more CSV files as input. These are passed to the cleaning stage, which validates, standardises, and merges them into a unified `cleaned.csv` for downstream stages.

Multiple files are accepted — they are concatenated before processing. Duplicate molecules (same canonical SMILES) are collapsed to one entry per unique molecule.

---

## Providing input files

In `request.json`, use any of these fields:

| Field | Type | Description |
|-------|------|-------------|
| `input_csv` | string or list | Single CSV path or list of CSV paths |
| `input_csvs` | string or list | Alias for `input_csv` |
| `input_folder` | string | Directory path — all `.csv` files inside are loaded |
| `input_folders` | string or list | List of directories |

These can be combined. All located CSV files are merged.

---

## Column reference

### Required

| Column name | Type | Description |
|-------------|------|-------------|
| `smiles` | string | SMILES representation of the molecule. Must be parseable by RDKit. The pipeline canonicalises SMILES and generates an InChIKey from it. |

### Recommended

| Column name | Type | Description |
|-------------|------|-------------|
| `mol_name` | string | Human-readable molecule name. Used in logs, summaries, and output tables. If absent, the InChIKey is used as the identifier. |

### Optional — charge and spin

| Column name | Type | Description |
|-------------|------|-------------|
| `charge` | integer | Formal charge of the molecule. Used for quantum chemistry calculations. |
| `multiplicity` | integer | Spin multiplicity (1 = singlet, 2 = doublet, 3 = triplet, ...). |

**What happens when these columns are absent:**

If `charge` is not present in the CSV:
- RDKit is used to compute the formal charge from the SMILES.
- If RDKit returns a non-zero charge, it is used (`charge_source: "rdkit_computed"`).
- If RDKit returns zero, the default of 0 is applied (`charge_source: "default"`).
- A `[WARNING]` is emitted in the stage log showing how many molecules fell into each category.

If `multiplicity` is not present:
- RDKit is used to check for radical electrons (open-shell character).
- Closed-shell molecules get `multiplicity = 1` (`multiplicity_source: "rdkit_computed"`).
- If RDKit detects radicals but cannot determine multiplicity, `multiplicity = 1` is used as a default (`multiplicity_source: "default"`), and a warning is emitted.

> If your molecules have unusual charges or spin states, always include `charge` and `multiplicity` columns explicitly. The pipeline will use them directly (`charge_source: "user_input"`) and no warning is emitted.

### Optional — melting point

| Column name | Type | Description |
|-------------|------|-------------|
| `melting_temp_c` | float | Melting point in **degrees Celsius**. Used for solubility saturation corrections via the Myrdal-Yalkowsky model. |
| `melting_temp_source` | string | Source label for the melting point (e.g. `"ChemBook"`, `"DSC"`, `"literature"`). Stored in metadata for traceability. |

**What happens when melting point is absent:**
- The saturation correction in the solubility stage uses a fallback mode (or is disabled, depending on `saturation.enabled` in stage args).
- A `[WARNING]` is emitted if `melting_temp_c` is missing.
- If a value is present but outside the physically plausible range (roughly −273 °C to 2000 °C), it is rejected with a warning and treated as missing.

### Optional — identifiers

These columns are carried through the pipeline for traceability but do not affect calculations.

| Column name | Type | Description |
|-------------|------|-------------|
| `inchi_key` | string | If provided, used as a consistency check against the RDKit-computed InChIKey. Mismatches are flagged. |
| `cas` | string | CAS registry number. Stored in metadata. |
| `source` | string | Data provenance label. Stored in metadata. |
| `notes` | string | Free-text notes. Stored in metadata. |

---

## Column name normalisation

The cleaning stage normalises column names before processing:
- Strips leading/trailing whitespace
- Converts to lowercase
- Replaces spaces and hyphens with underscores

So `Melting Temp (C)`, `melting-temp-c`, and `melting_temp_c` all map to `melting_temp_c`.

---

## Handling multiple files / duplicate molecules

When multiple CSV files are provided:
1. All files are concatenated into `raw_combined.csv`.
2. Each row is cleaned and standardised independently.
3. Duplicate InChIKeys (same canonical molecule) are collapsed — the first occurrence wins for conflicting fields.
4. A warning is emitted for each duplicate.

---

## Example minimal CSV

```csv
smiles,mol_name
c1ccccc1,benzene
CC(=O)O,acetic acid
c1ccc(cc1)N,aniline
```

## Example full CSV

```csv
smiles,mol_name,charge,multiplicity,melting_temp_c,melting_temp_source
c1ccccc1,benzene,0,1,5.5,literature
CC(=O)O,acetic acid,0,1,16.6,ChemBook
[NH4+],ammonium,1,1,,
CC([O-])=O,acetate,-1,1,,
```

---

## Common issues

**RDKit cannot parse the SMILES**
The row is skipped with a `[SKIP]` log entry. Check the SMILES string for validity using RDKit or a SMILES validator.

**Melting point out of range**
Values outside the plausible physical range are discarded and a `[WARNING]` is emitted. Verify the value is in Celsius (not Kelvin or Fahrenheit).

**Charge/multiplicity mismatch**
If a user-supplied charge or multiplicity is inconsistent with the SMILES (e.g. a clearly neutral molecule with `charge = 1`), a warning is emitted. The user-supplied value is still used.
