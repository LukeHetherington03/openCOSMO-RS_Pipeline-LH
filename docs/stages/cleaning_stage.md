# Stage: Cleaning

## Purpose
Normalise raw molecule data from CSV into a clean, canonical dataset and per‑molecule metadata, ready for all downstream 3D/conformer stages.

The CleaningStage:
- standardises column names,
- validates and canonicalises SMILES,
- generates InChIKeys,
- computes molecule‑level descriptors,
- writes one metadata JSON per molecule,
- produces a unified clean_dataset.csv.

No conformers are generated here.

---

## Inputs

### Arguments
- `input_folder` or `input_folders` (optional)
- `input_csv` (optional; string or list of strings)

At least one CSV must be provided via folders and/or explicit paths.

### Each CSV must contain
- **SMILES** (or fuzzy‑matched equivalent)
- **At least one name column**: `Name` or `mol_name_iupac`

---

## Processing

### 1. Resolve input sources
- Collect all CSV files from:
  - `input_folder` / `input_folders`
  - `input_csv` / list of CSVs
- Fail if no CSVs found.

### 2. Standardise column names (fuzzy matching)
Examples:
- `"Name"` → `mol_name`
- `"mol_name"` → `mol_name_iupac`
- `"Melting point"` → `melting_temp`
- `"Melting point source"` → `melting_temp_source`
- `"Hfus"` → `Hfus`
- `"Gfus_mode"` → `Gfus_mode`
- etc.

### 3. Validate required fields
- Missing `smiles` → error (fail in strict mode).
- Missing any name field → error (fail in strict mode).
- Missing `charge` → default to `0` with info log.

### 4. Canonicalise SMILES
- RDKit canonicalisation applied to all SMILES.

### 5. Generate InChIKey
- From canonical SMILES.
- Invalid SMILES are logged and counted.

### 6. Compute molecule‑level descriptors (NEW)
For each unique InChIKey:
- `rotatable_bonds`
- `molecular_weight`
- `heavy_atom_count`
- `hbond_donors`
- `hbond_acceptors`
- `tpsa`
- `logp`
- `aromatic_rings`

### 7. Write molecule metadata JSON
One file per molecule:
`molecule_metadata/<inchi_key>.json`

Metadata includes:
- `metadata_version`
- `inchi_key`
- `smiles`
- `mol_name`
- `mol_name_iupac`
- `charge`
- `multiplicity`
- **new descriptors** (rotatable bonds, MW, TPSA, etc.)
- `melting_temp`
- `melting_temp_source`
- `Hfus`
- `Gfus_mode`
- `formula_calcd`
- `anionic`
- `provenance`:
  - `source_file`
  - `cleaning_timestamp`
  - `pipeline_version`

### 8. Log warnings
- Missing melting_temp (affects solubility)
- Invalid SMILES
- Defaulted charge

A warning summary is printed at the end.

---

## Outputs

### 1. `clean_dataset.csv`
A unified dataset containing:
- `key_inchi`
- `smiles` (canonical)
- `charge`
- `mol_name` (if present)
- `mol_name_iupac` (if present)
- all other normalised fields

**No** `lookup_id` or `conf_num` columns appear here.

### 2. `molecule_metadata/<inchi_key>.json`
One JSON per molecule containing:
- `metadata_version`
- `inchi_key`
- `smiles`
- `mol_name`
- `mol_name_iupac`
- `charge`
- `multiplicity`
- **rotatable_bonds**
- **molecular_weight**
- **heavy_atom_count**
- **hbond_donors**
- **hbond_acceptors**
- **tpsa**
- **logp**
- **aromatic_rings**
- `melting_temp`
- `melting_temp_source`
- `Hfus`
- `Gfus_mode`
- `formula_calcd`
- `anionic`
- `provenance`:
  - `source_file`
  - `cleaning_timestamp`
  - `pipeline_version`

---

## Guarantees
- Every molecule with a valid SMILES has:
  - a canonical SMILES,
  - a valid InChIKey,
  - a metadata JSON file.
- All optional metadata fields exist (defaults applied if missing).
- Missing melting_temp is explicitly logged and counted.
- No conformers are generated.
- No lookup_id or conf_num appears in metadata or clean_dataset.csv.

---

## Failure Modes

### No CSV files found
→ Stage fails immediately.

### Missing smiles or name fields
- **Non‑strict mode:** row/file skipped, stage continues.
- **Strict mode:** stage fails immediately.

