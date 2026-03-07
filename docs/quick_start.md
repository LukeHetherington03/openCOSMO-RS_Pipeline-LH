# Quick Start Guide

openCOSMO-RS Pipeline

---

## Prerequisites

- Python 3.10+
- Access to the software binaries (ORCA, XTB, gXTB, CREST)
- Access to the openCOSMO-RS library files

---

## 1. Install software and resources

The pipeline requires external binaries and libraries that are **not** included in the repository. These must be placed in your **home directory**, not inside the project folder.

> **Why?** XTB and gXTB use Fortran fixed-length character buffers internally. If the path to their executables exceeds this fixed length, the path is silently truncated and the tool fails with a cryptic error. Keeping software in a short path like `~/software/` prevents this. Installing inside a deep project directory (e.g. `/home/user/projects/my-research/openCOSMO-RS_Pipeline-LH/software/`) will cause silent failures.

### Expected directory layout

```
~/software/
    crest/
        crest                    # CREST executable
    xtb/
        bin/
            xtb                  # XTB executable
    g-xtb-main/
        binary/
            gxtb                 # gXTB executable
    orca_6_1_1/
        orca                     # ORCA executable

~/resources/
    openCOSMO-RS_py/
        src/                     # Python openCOSMO-RS source
    openCOSMO-RS_cpp/
        bindings/                # C++ bindings
```

Obtain these from your institution or the respective project repositories and place them at the paths above.

---

## 2. Clone the pipeline repository

```bash
git clone <repository-url> ~/openCOSMO-RS_Pipeline-LH
cd ~/openCOSMO-RS_Pipeline-LH
```

---

## 3. Install Python dependencies

```bash
pip install -r requirements.txt
```

---

## 4. Add the CLI alias

Add the following line to your `~/.bashrc` (or `~/.zshrc`):

```bash
# openCOSMO-RS Pipeline
alias pl='python3 -m modules.cli.cli'
```

Then reload your shell:

```bash
source ~/.bashrc
```

The `pl` command must be run from the project root directory. It reads configuration from `config/paths.json` relative to the working directory.

---

## 5. Update the configuration

Edit `config/paths.json` to point to your installed software and resources.

```json
{
  "base_dir": "/home/<you>/openCOSMO-RS_Pipeline-LH/pipeline_data",

  "crest": {
    "home":       "/home/<you>/software/crest",
    "executable": "/home/<you>/software/crest/crest"
  },
  "gxtb": {
    "home":       "/home/<you>/software/g-xtb-main",
    "executable": "/home/<you>/software/g-xtb-main/binary/gxtb"
  },
  "xtb": {
    "home":       "/home/<you>/software/xtb",
    "executable": "/home/<you>/software/xtb/bin/xtb"
  },
  "orca": {
    "home":       "/home/<you>/software/orca_6_1_1",
    "executable": "/home/<you>/software/orca_6_1_1/orca"
  },
  "opencosmo": {
    "python_src":    "/home/<you>/resources/openCOSMO-RS_py/src",
    "cpp_bindings":  "/home/<you>/resources/openCOSMO-RS_cpp/bindings",
    "python_driver": "/home/<you>/openCOSMO-RS_Pipeline-LH/modules/solubility_engine/runCOSMO_RS_cpp.py"
  }
}
```

> **Reminder:** Keep all paths short. The `xtb` and `gxtb` entries are especially sensitive to long paths due to Fortran buffer truncation. Use `~/software/` or similar short paths, not nested project subdirectories.

Edit `config/resource_allocation.json` to match your machine:

```json
{
  "total_cores": 24,
  "max_cores": 20
}
```

Set `max_cores` to approximately 80% of physical cores to leave headroom for the OS.

---

## 6. Validate your environment

```bash
pl env check
```

This validates all configured executables and resource paths. Fix any errors before running calculations.

---

## 7. Start the queue worker

The pipeline uses a persistent background worker to execute requests.

```bash
pl q start
```

Check it is running:

```bash
pl q status
```

---

## 8. Prepare your input CSV

Your input CSV must contain at minimum:

| Column | Required | Description |
|--------|----------|-------------|
| `smiles` | Yes | SMILES string |
| `mol_name` | Recommended | Human-readable name |
| `charge` | Optional | Formal charge (integer). If absent, RDKit or default (0) is used. |
| `multiplicity` | Optional | Spin multiplicity. If absent, RDKit (closed-shell) or default (1) is used. |
| `melting_temp_c` | Optional | Melting point in Celsius (for solubility saturation corrections) |
| `melting_temp_source` | Optional | Source label for the melting point value |

If `charge` or `multiplicity` are absent, a pipeline warning is emitted describing how many molecules used RDKit-computed or default values.

---

## 9. Submit a request

Create a `request.json` file:

```json
{
  "request_name": "my_first_run",
  "input_csv": "/path/to/your/molecules.csv",
  "pipeline_sequence": ["cleaning", "generation", "optimisation", "orcacosmo", "solubility"],
  "stage_args": [
    {},
    {"engine": "rdkit", "num_confs": 20},
    {"engine": "gxtb_opt_normal"},
    {"default_basis": "TZVP"},
    {}
  ]
}
```

Submit it:

```bash
pl r submit /path/to/request.json
```

Monitor progress:

```bash
pl r list
pl r status <request-id>
pl r logs <request-id>
```

---

## 10. Retrieve results

When the pipeline completes, outputs are collected in:

```
pipeline_data/requests/<request-id>/final_outputs/
    solubility_results.json          # full machine-readable results
    solubility_results.csv           # flat per-molecule/solvent table
    solubility_human_summary.txt     # aligned human-readable table
    molecule_metadata/               # per-molecule JSON descriptors
    pipeline_warnings.txt            # all warnings from the run
```

---

## Common commands

```bash
pl q start                   # Start the worker
pl q stop                    # Graceful stop (running items finish)
pl q kill                    # Immediate stop
pl q status                  # Worker state and queue counts
pl q list                    # Pending and running requests

pl r list                    # All requests
pl r status <id>             # Pipeline state for a request
pl r logs <id>               # Stage logs
pl r submit <request.json>   # Submit a new request

pl env check                 # Full environment validation
pl env software              # Check executables only
pl env resources             # Check openCOSMO paths only
```

---

## Troubleshooting

**XTB or gXTB fails with a path error / truncated filename**
The executable or working directory path is too long. Move software to `~/software/` and keep the project in a short path.

**`pl` command not found**
You are either in the wrong directory or the alias is not loaded. Run from the project root and check `~/.bashrc`.

**Worker not processing requests**
Run `pl q status` to check the worker state. If it is not running, start it with `pl q start`.

**Stage fails immediately**
Run `pl env check` to verify all configured executables are accessible. Check `pl r logs <id>` for the specific error.
