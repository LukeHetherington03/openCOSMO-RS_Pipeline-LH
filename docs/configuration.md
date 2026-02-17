# Configuration Guide  
openCOSMO‑RS Pipeline

The openCOSMO‑RS Pipeline is fully **config‑driven**.  
This document explains how configuration works, how to modify it safely, and how the pipeline resolves settings at runtime.

---

# 1. Overview

The pipeline uses two layers of configuration:

1. **Global configuration**  
   Defined in `config/paths.json` and injected into every Request.

2. **Per‑stage configuration**  
   Defined in the pipeline specification (`pipeline_spec`) and passed to each stage as `args`.

The combination of these two layers determines:

- Which executables are used  
- Where data is stored  
- How each stage behaves  
- How provenance is recorded  
- How solubility calculations are performed  

---

# 2. Global Configuration (`paths.json`)

This file defines:

- Locations of external executables  
- Locations of COSMO‑RS bindings  
- Locations of CONSTANT_FILES  
- The base directory for all pipeline data  
- The Git version (injected automatically at runtime)

Example:

```json
{
  "base_dir": "/home/user/openCOSMO-RS_Pipeline/pipeline_data",

  "crest": {
    "home": "/home/user/software/crest",
    "executable": "/home/user/software/crest/crest"
  },

  "gxtb": {
    "home": "/home/user/software/g-xtb-main",
    "executable": "/home/user/software/g-xtb-main/binary/gxtb"
  },

  "xtb": {
    "home": "/home/user/software/xtb",
    "executable": "/home/user/software/xtb/bin/xtb"
  },

  "orca": {
    "home": "/home/user/software/orca_6_1_1",
    "executable": "/home/user/software/orca_6_1_1/orca"
  },

  "opencosmo": {
    "python_src": "/home/user/resources/openCOSMO-RS_py/src",
    "cpp_bindings": "/home/user/resources/openCOSMO-RS_cpp/bindings",
    "python_driver": "/home/user/openCOSMO-RS_Pipeline/modules/solubility_engine/runCOSMO_RS_cpp.py"
  },

  "constant_files": {
    "root": "/home/user/openCOSMO-RS_Pipeline/CONSTANT_FILES",
    "metadata_dir": "/home/user/openCOSMO-RS_Pipeline/CONSTANT_FILES/molecule_metadata",
    "solvent_dir": "/home/user/openCOSMO-RS_Pipeline/CONSTANT_FILES/solvents",
    "chemistry_dir": "/home/user/openCOSMO-RS_Pipeline/CONSTANT_FILES/chemistry"
  }
}
```

---

# 3. Git Version Injection

The pipeline automatically injects the Git version into the configuration at runtime:

```python
from modules.utils.git_version import get_git_version
config["pipeline_version"] = get_git_version()
```

This value appears in:

- Molecule metadata  
- Provenance blocks  
- Request metadata  

Example:

```
"pipeline_version": "main@abc1234"
```

This ensures full reproducibility.

---

# 4. Per‑Stage Configuration (`pipeline_spec`)

Each stage receives its own configuration block via `pipeline_spec`.

Example:

```python
pipeline_spec = [
    {"stage": "cleaning", "args": {"input_csv": input_csv, "overwrite_metadata": False}},
    {"stage": "generation", "args": {"engine": "rdkit", "num_confs": 10}},
    {"stage": "pruning", "args": {"n": 1}},
    {"stage": "optimisation", "args": {"engine": "xtb_opt_normal"}},
    {"stage": "orcacosmo", "args": {}},
    {"stage": "solubility", "args": {}},
]
```

Each `args` dictionary is merged with the global config and passed to the stage as:

```
self.parameters
self.config
```

---

# 5. How Configuration is Resolved

Each stage receives:

### `self.parameters`
- The per‑stage arguments from `pipeline_spec`
- Example: `{"engine": "xtb_opt_normal"}`

### `self.config`
- The entire global configuration from `paths.json`
- Plus the injected Git version

### `self.job`
- The job object, which includes:
  - request_id  
  - job_id  
  - stage name  
  - canonical input/output paths  

This separation ensures:

- Global config is stable  
- Stage config is explicit  
- No hidden defaults  
- Full reproducibility  

---

# 6. Canonical Inputs and Outputs

Each stage declares a **canonical output** using:

```python
self.set_stage_output("cleaned.csv")
```

The pipeline uses these canonical filenames to:

- Chain stages deterministically  
- Support resume  
- Support continuation  
- Support reproducibility  

The canonical filenames are defined in:

```
Job.STAGE_OUTPUTS
```

---

# 7. CONSTANT_FILES

The pipeline expects a directory containing:

```
CONSTANT_FILES/
    chemistry/
    solvents/
    molecule_metadata/
```

These files are referenced by:

- CleaningStage (metadata)
- ORCACOSMOStage (chemistry files)
- SolubilityStage (solvent files)

The location is defined in:

```
config["constant_files"]
```

---

# 8. Modifying Configuration Safely

### ✔ You may modify:
- Executable paths  
- CONSTANT_FILES paths  
- Base directory  
- Stage arguments  

### ✔ You should NOT modify:
- Canonical filenames  
- Stage input/output conventions  
- Internal pipeline directories  

### ✔ You must NOT modify:
- request.json (immutable record of user intent)
- job_state.json (managed by the pipeline)
- pipeline_state.json (managed by the pipeline)

---

# 9. Example: Full Configuration Flow

1. User runs `pl run`  
2. `paths.json` is loaded  
3. Git version is injected  
4. Request is created with:
   - global config  
   - pipeline_spec  
5. Each stage receives:
   - `self.parameters` (stage args)
   - `self.config` (global config)
6. Stages write canonical outputs  
7. Metadata includes:
   - Git version  
   - Request ID  
   - Job ID  
   - Timestamp  

---

# 10. Next Steps

Continue with:

- `docs/pipeline_architecture.md`  
- `docs/stages/cleaning.md`  
- `docs/resume_and_continuation.md`  

