# Stage: Solubility

## Purpose
The Solubility Stage performs COSMO‑RS solubility calculations for each molecule using
renumbered `.orcacosmo` conformers (required by the COSMO‑RS engine).  
It consumes:

- `orcacosmo_results.json` from the ORCA‑COSMO Stage  
- molecule metadata (melting point, SMILES, etc.)  
- solvent conformers from the global solvent library  

It produces:

- `solubility_results.json` (canonical, all molecules)
- `solubility_summary.csv` (human‑readable)
- per‑solute JSON files (`results/<inchi_key>.json`)
- per‑solute human summaries (`results/<inchi_key>_summary.txt`)
- renumbered conformer directories with mapping files

This stage never modifies conformer geometry or `.orcacosmo` content.

---

## Inputs

### Required
#### `summary_file`
A JSON file produced by ORCA‑COSMO Stage (`orcacosmo_results.json`), containing:

```json
{
  "lookup_id": "BQJCRH..._conf000",
  "inchi_key": "BQJCRHHNABKAKU-KBQPJGBKSA-N",
  "orcacosmo_path": "orcacosmo_outputs/BQJCRH..._conf000.orcacosmo"
}

Molecule metadata directory
Configured via:

Code
config["constant_files"]["metadata_dir"]
Each file: metadata/<inchi_key>.json

Used to retrieve:

SMILES

melting temperature

fusion enthalpy

Gfus model

Solvent conformer directory
Code
config["constant_files"]["solvent_dir"]/<solvent_name>/*.

orcacosmo
Parameters (args)
All parameters are optional unless stated.

Parameter	Type	Default	Description
summary_file	string	required	Path to ORCA‑COSMO results JSON
solvent_name	string	from defaults	Solvent to use (e.g., water)
temperature	float	from defaults	Temperature in Kelvin
SORcf	float	from defaults	COSMO‑RS SOR correction factor
calculations	string	"mixed_only"	COSMO‑RS calculation mode
parallel	bool	from defaults	Enable parallel execution
n_workers	int	from defaults	Number of parallel workers
Strict mode
Code
config["solubility"]["strict"] = true
If true:

any failure aborts the stage

missing metadata or conformers is fatal

Processing Logic
1. Load configuration
Read solubility defaults JSON

Determine solvent, temperature, SORcf

Determine parallel mode and worker count

Load metadata + solvent directories

2. Load ORCA‑COSMO summary
Group conformers by inchi_key:

Code
{
  "inchi_key": "...",
  "conformers": [
      ("lookup_id", "path/to/orcacosmo"),
      ...
  ]
}
3. For each solute (one per inchi_key)
a. Load metadata
From metadata/<inchi_key>.json or job‑local metadata.

b. Renumber solute conformers
COSMO‑RS requires sequential names:

Code
solute/<inchi_key>/<inchi_key>_c000.orcacosmo
solute/<inchi_key>/<inchi_key>_c001.orcacosmo
...
A mapping file is written:

Code
solute/<inchi_key>/mapping.json
Example:

json
{
  "c000": "BQJCRH..._conf000",
  "c001": "BQJCRH..._conf001"
}
c. Renumber solvent conformers
Copied from global solvent directory:

Code
solvent/<solvent_name>/<solvent_name>_c000.orcacosmo
...
d. Run COSMO‑RS wrapper
Inputs:

solute directory

solvent directory

SMILES

melting point

temperature

SORcf

Wrapper returns:

raw COSMO‑RS output

parsed result dict

number of solute/solvent conformers

e. Validate wrapper output
Ensure:

Code
result["x_solubility"]
exists.

f. Write outputs
Raw output → results/<inchi_key>_raw.txt

Machine‑readable JSON → results/<inchi_key>.json

Human summary → results/<inchi_key>_summary.txt

g. Record result
Appended to in‑memory results list for global summary.

Parallel Execution
Parallel execution is controlled by:

solubility_defaults.json

stage arguments (parallel, n_workers)

When enabled, each solute is processed in its own worker process using ProcessPoolExecutor.

Renumbering is isolated per solute → safe for parallelism.

Outputs
1. solubility_results.json
Canonical list of all solubility results:

json
[
  {
    "inchi_key": "...",
    "solubility_x": 1.23e-4,
    "json": "results/<inchi_key>.json",
    "summary": "results/<inchi_key>_summary.txt",
    "raw_output": "results/<inchi_key>_raw.txt",
    "n_confs": 12
  }
]
2. solubility_summary.csv
Human‑readable summary:

| inchi_key | solubility_x | n_confs | json | summary | raw_output |

3. Per‑solute JSON
Code
results/<inchi_key>.json
Contains:

solubility

metadata

renumbered mapping

provenance

4. Per‑solute human summary
Code
results/<inchi_key>_summary.txt
5. Renumbered conformer directories
Code
cosmo/solute/<inchi_key>/
cosmo/solvent/<solvent_name>/
Each solute directory contains:

<inchi_key>_c000.orcacosmo

<inchi_key>_c001.orcacosmo

mapping.json

Provenance
Each per‑solute JSON includes:

Code
"provenance": {
    "cosmors_version": "...",
    "parallel": true/false,
    "n_workers": 4,
    "elapsed_seconds": 12.3,
    "last_modified": "2026‑02‑01T23:40:00Z"
}
Guarantees
Renumbering is preserved exactly as required by COSMO‑RS.

Original conformer identity is preserved in mapping.json.

All results are aggregated into a single canonical JSON.

Strict mode enforces fail‑fast behaviour.

No geometry or .orcacosmo content is modified.

Solvent conformers are copied and renumbered consistently.

Failure Modes
The stage fails if:

summary_file missing or invalid

solvent conformers missing

wrapper output incomplete

strict mode enabled and any solute fails

In non‑strict mode:

failures are logged

solute is skipped

stage continues

Warning Summary (logged)
At the end of the stage:

Code
========== Solubility Summary ==========
Successful: X
Failed: Y
Missing metadata: Z
Missing conformers: N
========================================