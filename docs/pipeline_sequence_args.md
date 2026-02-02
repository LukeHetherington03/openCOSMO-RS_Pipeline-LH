openCOSMO‑RS Pipeline  
Pipeline Sequence Specification & Stage Argument Guide
======================================================

This document defines all stages, arguments, defaults, and behaviours for
constructing a valid pipeline specification in the openCOSMO‑RS Pipeline.

A pipeline specification is a JSON file describing the ordered list of stages:

{
  "sequence": [
    { "stage": "cleaning", "args": {...} },
    { "stage": "generation", "args": {...} },
    { "stage": "pruning", "args": {...} },
    { "stage": "optimisation", "args": {...} },
    { "stage": "orcacosmo" },
    { "stage": "solubility", "args": {...} }
  ]
}

All heavy configuration (paths to binaries, COSMO constants, solubility defaults,
ORCA settings, etc.) is loaded from:

  config/paths.json
  config/solubility_defaults.json

Users only specify what they want to override.

----------------------------------------------------------------------
1. cleaning
----------------------------------------------------------------------

Purpose:
    Load molecules from CSV, sanitise, normalise, deduplicate.

Required arguments:
    input_csv : path to CSV containing SMILES or identifiers

Optional arguments:
    None.

Notes:
    - No conformer generation here.
    - No engine or seed needed.
    - Produces a cleaned molecule list.
    - Cleaning produces molecule metadata (InChI key, SMILES, charge,
      multiplicity, etc.).
    - Metadata is stored in molecule_metadata/<inchi_key>.json
    - If you are NOT using cleaning, you should use "continue_from" to maintain
      provenance. You do NOT need to specify the summary file manually; it is
      auto‑detected by the pipeline.

----------------------------------------------------------------------
2. generation
----------------------------------------------------------------------

Purpose:
    Generate conformers using one of several engines.

Supported engines:
    rdkit
    crest
    babel

Required arguments:
    None.

Optional arguments:
    num_confs     (default = 20)      Number of conformers
    engine        (default = "rdkit") Backend
    seed          (default = 42)      RNG seed (rdkit, crest)
    crest_flags   (default = none)    Additional CREST flags
    babel_method  (default = "confab") Babel conformer method

Notes:
    - Generation produces molecule metadata (InChI key, SMILES, charge,
      multiplicity, etc.).
    - All binaries (CREST, Babel) are read from config.

----------------------------------------------------------------------
3. pruning
----------------------------------------------------------------------

Purpose:
    Reduce conformer set based on energy or structural similarity.

Methods:
    topN
    energy_window
    boltzmann (future)

Required arguments:
    None.

Optional arguments (method‑specific):

  topN:
      max_confs (default = 1)

  energy_window:
      energy_cutoff   (default = 5.0 kcal/mol)
      rmsd_threshold  (default = none)

  boltzmann (future):
      boltzmann_temp       (default = 298.15 K)
      energy_gap_threshold (default = 10.0 kcal/mol)

Notes:
    - Default pruning behaviour: topN with max_confs = 1.
    - Boltzmann pruning will require both temperature and an energy‑gap cutoff
      to avoid pathological weighting.

----------------------------------------------------------------------
4. optimisation
----------------------------------------------------------------------

Purpose:
    Optimise geometries using a chosen backend.

Supported engines:
    forcefield
    xtb
    gxtb
    dft

Required arguments:
    engine

Optional arguments:
    max_iter (default = 250)
    level    (default = "normal")
             Allowed levels:
                 normal
                 loose
                 tight
                 vtight
                 vvtight

Backend behaviour:

  forcefield:
      RDKit MMFF/UFF

  xtb:
      Uses xtb binary from config
      Levels: normal, tight, loose

  gxtb:
      Uses xtb driver mode with gXTB binary
      Levels: normal, tight, loose, vtight, vvtight

  dft:
      Uses ORCA
      Levels: low, medium, high (mapped to ORCA keywords)

Notes:
    - All binaries (xtb, gxtb, orca) are read from config.
    - gXTB uses the validated driver command:
        xtb molecule.xyz --driver "<gxtb> -grad -c xtbdriver.xyz" --opt

----------------------------------------------------------------------
5. orcacosmo
----------------------------------------------------------------------

Purpose:
    Compute COSMO surfaces + polarizabilities using ORCA.

Required arguments:
    None.

Optional arguments:
    None.

Behaviour:
    - ORCA executable from config.
    - Always compute polarizabilities.
    - Never optimise geometry.
    - COSMO settings from internal defaults.

Notes:
    - This stage is fully automatic.
    - Users should not override anything.

----------------------------------------------------------------------
6. solubility
----------------------------------------------------------------------

Purpose:
    Compute solubility using openCOSMO‑RS.

Defaults:
    All defaults come from:
        config/solubility_defaults.json

    Example defaults file:

    {
      "default_solvent": "water",
      "calculation_type": "all",
      "temperature": 298.15,

      "saturation": {
        "enabled": true,
        "Gfus_mode": "MyrdalYalkowsky",
        "Hfus": "N/A",
        "SORcf": 1
      },

      "mixture": {
        "default_multiplicity": 1,
        "force_absolute_paths": true
      }
    }

Required arguments:
    None.

Optional arguments:
    solvent_name   (default = "water")
    temperature    (default = 298.15)

Advanced overrides (not recommended):
    calculation_type
    SORcf
    Hfus
    Gfus_mode
    default_multiplicity
    force_absolute_paths

Notes:
    - All defaults come from solubility_defaults.json.
    - Code does NOT hard‑code defaults.
    - Only solvent_name and temperature are expected to be overridden.
    - All binaries (openCOSMORS) are read from config.

----------------------------------------------------------------------
7. Example Pipeline Specs
----------------------------------------------------------------------

Minimal gXTB optimisation:

{
  "sequence": [
    { "stage": "cleaning", "args": { "input_csv": "mols.csv" } },
    { "stage": "generation" },
    { "stage": "pruning" },
    { "stage": "optimisation", "args": { "engine": "gxtb" } }
  ]
}

Full workflow:

pipeline_spec = 
[
    { "stage": "cleaning", "args": { "input_csv": "mols.csv" } },
    { "stage": "generation", "args": { "num_confs": 50, "engine": "rdkit" } },
    { "stage": "pruning", "args": { "method": "energy_window", "energy_cutoff": 3.0 } },
    { "stage": "optimisation", "args": { "engine": "gxtb", "level": "tight" } },
    { "stage": "orcacosmo" },
    { "stage": "solubility", "args": { "solvent_name": "methanol" } }
]


----------------------------------------------------------------------
8. Provenance Note
----------------------------------------------------------------------

If you are NOT using cleaning + input_csv, you should use:

    Request.continue_from(...)

This preserves provenance and ensures the correct summary file is
auto‑detected. You do NOT need to manually specify the summary file in
your pipeline spec.

----------------------------------------------------------------------
END OF DOCUMENT
----------------------------------------------------------------------
