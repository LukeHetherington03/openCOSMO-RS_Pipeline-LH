# Installation Guide  
openCOSMO‑RS Pipeline

This document explains how to install the external scientific tools and Python dependencies required to run the openCOSMO‑RS Pipeline.  
The pipeline is designed for Linux HPC environments and assumes a standard Python 3 installation with manually installed chemistry engines.

---

# 1. System Requirements

## Operating System
- Linux (recommended)
- Tested on HPC clusters and standard Linux workstations

## Hardware
- Multi‑core CPU (ORCA, XTB, and gXTB scale with threads)
- ≥ 16 GB RAM recommended
- ≥ 50 GB disk space for ORCA scratch + pipeline outputs

---

# 2. Add the Pipeline CLI Alias (Recommended)

The pipeline includes a full command‑line interface.  
To enable the `pl` command, add this to your `.bashrc`:

```bash
# openCOSMO‑RS Pipeline CLI
alias pl='python3 -m modules.cli.cli'
```

Reload your shell:

```bash
source ~/.bashrc
```

You can now run commands like:

```bash
pl run
pl env pip
pl env xtb
pl queue status
pl help
```

No other shell configuration is required — all executables are referenced directly via `config/paths.json`.

---

# 3. Python Dependencies

The pipeline uses a small set of Python libraries.  
Install them with:

```bash
pip install numpy scipy pandas rdkit-pypi openbabel
```

To verify your environment, run:

```bash
pl env pip
```

Example output:

```
======================================================================
Pip Package Validation
======================================================================
✓ numpy                     version 1.26.4
✓ scipy                     version 1.14.1
✓ pandas                    version 2.3.3
✓ rdkit                     version 2025.09.1
✓ openbabel                 version 3.1.0
```

If all checks pass, your Python environment is ready.

---

# 4. Required External Tools

The pipeline depends on the following external executables and libraries:

| Component | Purpose | Required |
|----------|----------|----------|
| **ORCA 6.x** | DFT + COSMO calculations | ✔ |
| **XTB** | Fast semi‑empirical optimisation | ✔ |
| **gXTB** | Gradient driver for XTB | ✔ |
| **CREST** | Conformer sampling | Optional |
| **COSMO‑RS Python bindings** | Solubility prediction | ✔ |
| **COSMO‑RS C++ bindings** | High‑performance solubility engine | ✔ |
| **Git** | Version provenance | ✔ |

All tools are referenced directly in `config/paths.json`.  
No environment modules or shell exports are required.

---

# 5. ORCA Installation

ORCA is required for:
- TZVPD / TZVP COSMO calculations  
- ORCA‑based optimisation (optional)

### Download ORCA
Official download (registration required):  
https://orcaforum.kofo.mpg.de

### Install
Extract ORCA to a directory such as:

```
/home/<user>/software/orca_6_1_1/
```

Ensure the executable exists:

```
/home/<user>/software/orca_6_1_1/orca
```

---

# 6. XTB Installation

XTB is used for:
- Fast geometry optimisation  
- gXTB driver integration  

### Repository
https://github.com/grimme-lab/xtb

### Installation
Follow the instructions in the repo README.  
You will obtain a compiled binary:

```
xtb
```

Place it somewhere like:

```
/home/<user>/software/xtb/bin/xtb
```

---

# 7. gXTB Installation

gXTB is required for:
- Gradient‑based optimisation using XTB  
- The `--driver` interface used by the pipeline  

### Repository
https://github.com/grimme-lab/gxtb

### Installation
Follow the README instructions to compile gXTB.  
You will obtain:

```
gxtb
```

Place it somewhere like:

```
/home/<user>/software/g-xtb-main/binary/gxtb
```

---

# 8. CREST Installation (Optional)

CREST is used only for conformer generation if selected.

### Repository
https://github.com/grimme-lab/crest

### Installation
Follow the README instructions.  
You will obtain:

```
crest
```

Place it somewhere like:

```
/home/<user>/software/crest/crest
```

---

# 9. COSMO‑RS Bindings (Python + C++)

The solubility stage requires:
- Python bindings for COSMO‑RS  
- C++ shared libraries  
- A Python driver script that calls the C++ engine  

### Repositories (your internal versions)
You must provide:

```
python_src      → COSMO‑RS Python bindings
cpp_bindings    → COSMO‑RS C++ shared libraries (.so)
python_driver   → runCOSMO_RS_cpp.py
```

Example layout:

```
/home/<user>/resources/openCOSMO-RS_py/src
/home/<user>/resources/openCOSMO-RS_cpp/bindings
/home/<user>/openCOSMO-RS_Pipeline/modules/solubility_engine/runCOSMO_RS_cpp.py
```

---

# 10. CONSTANT_FILES Directory

The pipeline requires a directory containing:

- `chemistry/` (CPCM radii, etc.)
- `solvents/` (COSMO‑RS solvent files)
- `molecule_metadata/` (global metadata store)

Example:

```
/home/<user>/openCOSMO-RS_Pipeline/CONSTANT_FILES/
```

This directory is referenced in `paths.json`.

---

# 11. Configure paths.json

The pipeline uses a single configuration file:

```
config/paths.json
```

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

# 12. Verification

Test ORCA:

```bash
/home/user/software/orca_6_1_1/orca -h
```

Test XTB:

```bash
xtb --version
```

Test gXTB:

```bash
gxtb -h
```

Test CREST:

```bash
crest --help
```

Test COSMO‑RS Python bindings:

```bash
python3 -c "import cosmo_rs; print('OK')"
```

---

# 13. Next Steps

Once installation is complete, continue with:

- `docs/configuration.md`  
- `docs/pipeline_architecture.md`  
- `docs/stages/cleaning.md`  
- `docs/resume_and_continuation.md`  

