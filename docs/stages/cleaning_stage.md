# Cleaning Stage  
openCOSMO‑RS Pipeline

The Cleaning Stage is the first scientific stage in the pipeline.  
It ingests raw CSV files, standardises them, validates required fields, canonicalises SMILES, generates InChIKeys, computes extensive molecule metadata, and produces a unified dataset suitable for downstream conformer generation and optimisation.

This document describes the Cleaning Stage in detail for both users and developers.

---

# 1. Purpose

The Cleaning Stage:

- Accepts one or more raw CSV files  
- Normalises column names  
- Validates required fields  
- Canonicalises SMILES  
- Generates InChIKeys  
- Computes physchem descriptors  
- Detects functional groups  
- Writes molecule‑level metadata (local + global)  
- Produces a unified `cleaned.csv` dataset  

It is the foundation for all downstream stages.

---

# 2. Canonical Input & Output

### **Input**
The stage accepts raw CSV files via:

- `input_csv`
- `input_csvs`
- `input_folder`
- `input_folders`

At least one must be provided.

### **Output**
The stage declares its canonical output:

```
cleaned.csv
```

and also writes:

```
molecule_metadata/<inchi_key>.json
```

These appear under:

```
jobs/J-.../outputs/
```

---

# 3. Parameters

The Cleaning Stage reads the following parameters from `pipeline_spec`:

| Parameter | Type | Required | Description |
|----------|------|----------|-------------|
| `input_csv` | str or list | optional | One or more CSV files |
| `input_folder` | str | optional | Folder containing CSV files |
| `overwrite_metadata` | bool | optional | Whether to overwrite global metadata |
| `stage_input` | unused | no | (Not used in this stage) |

At least one of `input_csv` or `input_folder` must be provided.

---

# 4. Strict Mode

Strict mode is enabled via:

```json
"cleaning": { "strict": true }
```

in `paths.json`.

Strict mode enforces:

- Missing SMILES → fail  
- Missing name column → fail  

Without strict mode, these become warnings and the file is skipped.

---

# 5. Execution Flow

The Cleaning Stage follows this sequence:

1. **Resolve input CSVs**  
   - Accepts folders or explicit file paths  
   - Validates existence  
   - Writes a merged `raw_combined.csv` to job inputs  

2. **Track items**  
   Each CSV file is treated as one “item” for resume logic.

3. **For each CSV file:**  
   - Load the file  
   - Standardise headers (fuzzy matching)  
   - Validate required fields  
   - Canonicalise SMILES  
   - Generate InChIKeys  
   - Compute metadata  
   - Write metadata JSON files  
   - Reorder columns  
   - Append to combined dataset  

4. **Write final cleaned dataset**  
   - Concatenate all processed frames  
   - Write `cleaned.csv`  

5. **Write warning summary**  

---

# 6. Header Standardisation

The stage uses fuzzy matching to map raw column names to standard names.

Example:

```
"Melting Point" → "melting_temp"
"Name" → "mol_name"
"experimental solubility /mol frac" → "experimental_solubility_mol_frac"
```

This is implemented via:

- `_normalise_header()`
- `_fuzzy_match()`
- `STANDARD_NAMES`

---

# 7. Required Fields

A CSV must contain:

- `smiles`
- `mol_name` or `mol_name_iupac`

If `charge` is missing:

- It is added with default value `0`
- A warning is recorded

---

# 8. SMILES Canonicalisation

SMILES are canonicalised using RDKit:

```python
Chem.MolToSmiles(mol, canonical=True)
```

Invalid SMILES are retained but flagged.

---

# 9. InChIKey Generation

For each SMILES:

```python
inchi.InchiToInchiKey(inchi.MolToInchi(mol))
```

Invalid SMILES → `inchi_key = None`.

These rows are retained but logged.

---

# 10. Metadata Generation

For each unique InChIKey, the stage computes:

### **Physchem descriptors**
- Rotatable bonds  
- Exact molecular weight  
- Heavy atom count  
- H‑bond donors/acceptors  
- TPSA  
- LogP  
- Aromatic rings  
- Fraction sp3  
- Bertz complexity  
- Molar refractivity  

### **Structural counts**
- N, O, S, halogens  
- Ring count  
- Heterocycle count  
- Double/triple bonds  

### **Functional groups**
Detected via SMARTS patterns:

- Amide  
- Ester  
- Alcohol  
- Amine  
- Carboxylic acid  
- Ketone  
- Aldehyde  
- Nitrile  
- Acrylate  
- Methacrylate  
- Vinyl ether  
- Vinyl ester  
- Styrenic  
- Ether  
- Carbonate  
- Urethane  
- Urea  
- Sulfonamide  
- Sulfonate  
- Phosphate  
- Halogenation flags  

Each group produces:

- A boolean flag (`is_amide`)  
- A count (`amide_count`)  

### **Provenance block**
Each metadata file includes:

```
source_file
cleaning_timestamp
pipeline_version
generated_by_request
generated_by_job
```

### **Local vs Global Metadata**
- Local metadata is always written  
- Global metadata is written only if:
  - `overwrite_metadata=True`, or  
  - the file does not already exist  

---

# 11. Warning System

The stage tracks:

- Missing melting temperature  
- Invalid SMILES  
- Defaulted charge  

At the end, a summary is printed:

```
Cleaning completed with warnings:
 - 3 molecule(s) missing melting_temp
 - 1 invalid SMILES encountered
 - 'charge' defaulted to 0 for 2 file(s)
```

---

# 12. Failure Modes

The stage fails if:

- No CSV files are provided  
- No CSV files can be read  
- Required fields missing (strict mode)  
- No valid rows processed  

Failures update:

```
pipeline_state.json
job_state.json
```

and halt the pipeline.

---

# 13. Developer Notes

### **Resume Support**
Each CSV file is treated as one “item”.  
If the stage is interrupted:

- Completed files are skipped  
- Pending files are processed  
- Metadata is not duplicated  

### **Canonical Output**
The stage declares:

```python
self.set_stage_output("cleaned.csv")
```

This is essential for:

- Stage chaining  
- Continuation  
- Reproducibility  

### **Global Metadata**
Global metadata lives in:

```
CONSTANT_FILES/molecule_metadata/
```

This allows cross‑request reuse.

---

# 14. Minimal Example

### pipeline_spec entry

```python
{
    "stage": "cleaning",
    "args": {
        "input_folder": "raw_data/",
        "overwrite_metadata": False
    }
}
```

### Running the stage

```
python3 main.py
```

---

# 15. Summary

The Cleaning Stage:

- Normalises raw data  
- Validates structure  
- Canonicalises SMILES  
- Generates InChIKeys  
- Computes rich metadata  
- Writes local + global metadata  
- Produces a unified cleaned dataset  

It is the foundation for all downstream scientific computation.

