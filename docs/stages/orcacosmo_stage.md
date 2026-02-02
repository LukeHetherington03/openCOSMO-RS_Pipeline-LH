# Stage: ORCA COSMO

## Purpose
The ORCA COSMO Stage performs single‑point CPCM calculations using ORCA (TZVPD with optional TZVP fallback) to generate COSMO surface data and `.orcacosmo` files for each conformer.

It consumes:
- a conformer summary (energies‑like JSON) containing lookup_id, inchi_key, xyz_path
- molecule metadata (for SMILES lookup)

It produces:
- `.orcacosmo` files (one per conformer)
- `orcacosmo_results.json` (all results in one JSON)
- `orcacosmo_summary.csv` (human‑readable summary)
- `item_to_lookup_mapping.json`
- raw ORCA outputs and parser JSONs for debugging

This stage never modifies molecule‑level metadata or conformer geometry.

---

## Inputs

### Required
#### `summary_file`
A JSON list of conformer entries, each containing:

```json
{
  "lookup_id": "BQJCRH..._conf000",
  "inchi_key": "BQJCRHHNABKAKU-KBQPJGBKSA-N",
  "xyz_path": "xyz/BQJCRH..._conf000_opt.xyz",
  "energy": -13.456,
  "smiles": "CC(O)COC(=O)C=C"
}
Molecule metadata directory
Configured via:

Code
config["constant_files"]["metadata_dir"]
Each file: metadata/<inchi_key>.json

Used to retrieve canonical SMILES.

CPCM radii file
Code
config["constant_files"]["chemistry_dir"]/cpcm_radii.json
Contains:

element radii

cut_area

Parameters (args)
Parameter	Type	Default	Description
summary_file	string	required	Conformer summary JSON
enable_tzvp_fallback	bool	false	If TZVPD fails, try TZVP
strict (config)	bool	false	If true, any failure aborts the stage
Processing Logic
1. Load stage configuration
ORCA executable path

metadata directory

CPCM radii

fallback settings

ORCA version (best‑effort)

2. Load conformer XYZs
Read summary_file

Copy each conformer’s XYZ into inputs/ as <lookup_id>.xyz

3. Discover lookup_ids
All .xyz files in inputs/ become processing targets.

4. For each conformer
Assigned a short ORCA‑safe item number:

Code
item000, item001, item002, ...
Then:

a. Prepare workdir
Code
tmp_exec/item000/
b. Write ORCA input (TZVPD)
Short filenames:

Code
item000_tzvpd.inp
item000.log
item000.cpcm
item000.cpcm_corr
c. Run ORCA TZVPD
Check for valid .cpcm, .cpcm_corr, .log

If valid → continue

If invalid → fallback (if enabled)

d. Optional fallback: TZVP
Write item000_tzvp.inp

Run ORCA

Validate outputs

If still invalid → fail (strict) or skip (non‑strict)

e. Copy raw outputs
Copied to:

Code
raw_outputs/<lookup_id>.log
raw_outputs/<lookup_id>.cpcm
raw_outputs/<lookup_id>.cpcm_corr
f. Parse outputs
Using:

OrcaLogParser

OrcaCpcmParser

OrcaCpcmCorrParser

Parser JSONs are written to:

Code
parsed_outputs/<lookup_id>.log.json
parsed_outputs/<lookup_id>.cpcm.json
parsed_outputs/<lookup_id>.cpcm_corr.json
g. Build bundle for orchestrator
Bundle contains:

metadata (lookup_id, inchi_key, smiles, method_used, fallback)

paths to parser JSONs

XYZ path

h. Run OrcaCosmoOrchestrator
Produces:

Code
orcacosmo_outputs/<lookup_id>.orcacosmo
i. Record success
Stored in memory for final JSON + CSV.

Outputs
1. orcacosmo_outputs/<lookup_id>.orcacosmo
Final COSMO surface file reconstructed by the orchestrator.

2. orcacosmo_results.json
A single JSON containing all conformer results:

json
[
  {
    "lookup_id": "...",
    "inchi_key": "...",
    "smiles": "...",
    "item_number": 0,
    "method_used": "TZVPD",
    "fallback_triggered": false,
    "elapsed_seconds": 12.4,
    "orcacosmo_path": "orcacosmo_outputs/<lookup_id>.orcacosmo",
    "parsed": {
      "log": { ... },
      "cpcm": { ... },
      "cpcm_corr": { ... }
    },
    "provenance": {
      "orca_version": "6.1.1",
      "cpcm_radii_source": ".../cpcm_radii.json",
      "last_modified": "2026-02-01T23:40:00Z"
    }
  }
]
3. orcacosmo_summary.csv
Human‑readable summary:

| lookup_id | inchi_key | method_used | fallback_triggered | elapsed_seconds | orcacosmo_path | item_number |

4. item_to_lookup_mapping.json
Maps ORCA item numbers to lookup_ids:

json
{
  "0": "BQJCRH..._conf000",
  "1": "BQJCRH..._conf001"
}
5. Raw outputs (for debugging)
Code
raw_outputs/<lookup_id>.log
raw_outputs/<lookup_id>.cpcm
raw_outputs/<lookup_id>.cpcm_corr
6. Parser JSONs (for orchestrator)
Code
parsed_outputs/<lookup_id>.log.json
parsed_outputs/<lookup_id>.cpcm.json
parsed_outputs/<lookup_id>.cpcm_corr.json
Guarantees
ORCA filenames inside workdirs remain short (item000.*) to avoid filesystem limits.

External filenames use full lookup_id for clarity.

All results are aggregated into a single JSON file.

.orcacosmo files are always written if ORCA succeeds.

Fallback is only used if enabled.

No molecule‑level metadata is modified.

No geometry is changed.

Failure Modes
The stage fails if:

summary_file missing or invalid

CPCM radii file missing

ORCA executable missing

Both TZVPD and TZVP fail (strict mode)

XYZ missing for a conformer (strict mode)

In non‑strict mode:

failures are logged,

conformer is skipped,

stage continues.

Warning Summary (logged)
At the end of the stage:

Code
========== ORCA COSMO Summary ==========
Successful: X
Failed: Y
Fallback used: Z
Missing XYZ: N
========================================