# Stage: Optimisation

## Purpose
The OptimisationStage performs geometry optimisation on each conformer using a selected quantum or semi‑empirical backend.

It consumes a canonical `energies.json` (from Generation or Pruning) and produces:

- updated `energies.json` with optimisation provenance,
- optimised XYZ geometries,
- a human‑readable `optimisation_summary.csv`,
- a conformer‑level `summary.csv`,
- an updated `job_state.json`.

This stage never modifies molecule‑level metadata.

---

## Inputs

### Required

#### `inputs/energies.json`
A canonical conformer store produced by Generation or Pruning.

Each entry must match the `ConformerRecord` schema:

```json
{
  "lookup_id": "BQJCRH..._conf000",
  "inchi_key": "BQJCRHHNABKAKU-KBQPJGBKSA-N",
  "conf_num": 0,
  "xyz_path": "xyz/BQJCRH..._conf000.xyz",
  "energy": -12.345,
  "smiles": "CC(O)COC(=O)C=C",
  "provenance": { ... }
}
inputs/xyz/
Directory containing all XYZ files referenced in energies.json.

Parameters (args)
All parameters are optional unless stated.

Parameter	Type	Default	Description
engine	string	"gxtb"	Optimisation backend (xtb, gxtb, orca, orca_fast, orca_final, forcefield)
level	string	"normal"	Normalised optimisation level: loose, normal, tight, vtight
max_iter	int	250	Maximum optimisation iterations (XTB/gXTB)
global_fail_threshold	float	0.8	Abort if this fraction of conformers fail
strict (config)	bool	false	If true, any fatal issue aborts the stage
Normalised Optimisation Levels
The stage exposes a unified optimisation vocabulary:

Code
loose
normal
tight
vtight
gXTB mapping
Level	Flag
loose	--opt loose
normal	--opt
tight	--opt tight
vtight	--opt vtight
XTB mapping
Level	Iterations
loose	100
normal	250
tight	500
vtight	1000
ORCA mapping
Level	Keyword
loose	LooseOpt
normal	Opt
tight	TightOpt
vtight	VeryTightOpt
Forcefield
No effect — placeholder optimisation.

Processing Logic
1. Load conformers
Read energies.json into a ConformerSet.

Validate schema.

Group conformers by molecule.

2. Prepare XYZs
Copy all XYZs into inputs/xyz/.

Track missing XYZs.

If all XYZs missing → fail.

3. Optimise each conformer
For each conformer:

Run backend (XTB/gXTB/ORCA/forcefield).

Parse log file (energy, convergence, iterations).

Validate energy sanity.

Normalise XYZ filename → <lookup_id>_opt.xyz.

Compute:

convergence quality,

success score,

geometry version,

last_modified timestamp.

Record warnings.

Update ConformerRecord.provenance["optimisation"].

4. Global failure threshold
If failures exceed global_fail_threshold, abort optimisation early.

5. Usability filtering
A conformer is usable if:

XYZ exists,

energy is finite (if converged),

status is converged or partial.

If all conformers unusable → fail.

Outputs
All outputs are written to outputs/.

1. outputs/xyz/<lookup_id>_opt.xyz
Optimised geometries with normalised filenames.

2. outputs/summary.csv
Conformer‑level table:

| lookup_id | inchi_key | conf_num | energy | status | quality | success_score | xyz_path | log_path | elapsed_seconds |

3. outputs/optimisation_summary.csv
Molecule‑level summary:

| molecule_id | n_atoms | n_conformers_attempted | n_conformers_output | n_converged | n_failed | n_partial | avg_time_per_conf_s | total_time_s | engine | level |

4. outputs/energies.json
Canonical conformer store with updated optimisation provenance.

Example provenance block:

json
"optimisation": {
  "engine": "gxtb",
  "level": "tight",
  "status": "converged",
  "energy": -13.456,
  "xyz_path": "xyz/BQJCRH..._conf000_opt.xyz",
  "log_path": "log/BQJCRH..._gxtb.log",
  "elapsed_seconds": 12.3,
  "parser": {
    "iterations": 45,
    "converged": true,
    "elapsed_seconds": 11.9
  },
  "backend_meta": { ... },
  "engine_command": "...",
  "backend_version": "gXTB 6.5.1",
  "warnings": [],
  "convergence_quality": "very_good",
  "success_score": 0.8,
  "geometry_version": 1,
  "last_modified": "2026-02-01T23:40:00Z"
}
5. job_state.json
Updated with:

stage name,

engine,

level,

counts,

elapsed time.

Warning Summary (logged)
At the end of the stage, the log prints:

Code
========== Optimisation Summary ==========
Converged: X
Partial: Y
Failed: Z
Missing XYZ: N
==========================================
Failure Modes
The stage fails if:

energies.json missing or invalid,

all XYZs missing,

all conformers unusable,

strict mode enabled and any fatal error occurs,

global failure threshold exceeded.

Guarantees
Conformer identity is preserved (lookup_id, inchi_key, conf_num).

Only geometry and energy fields are updated.

Provenance is always appended, never overwritten.

Output energies.json is canonical and schema‑consistent.

XYZ filenames are normalised.

Molecule‑level metadata is never modified.

Versioning
The stage increments:

geometry_version per conformer,

last_modified timestamp.

This ensures reproducibility and traceability.