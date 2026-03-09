# openCOSMO-RS Pipeline

A fully reproducible, config-driven workflow for predicting molecular solubility using conformer generation, quantum chemistry, and COSMO-RS modelling.

---

## What it does

The pipeline automates the full solubility prediction workflow:

1. **Cleaning** — validate SMILES, generate InChIKeys, resolve charge/multiplicity/melting point, compute physchem descriptors
2. **Generation** — generate 3D conformers (RDKit, CREST, or OpenBabel)
3. **Pruning** — filter conformers by energy window and RMSD (optional)
4. **Optimisation** — geometry optimise conformers (gXTB, XTB, ORCA DFT, or force fields)
5. **ORCA COSMO** — run ORCA CPCM single-point calculations, produce `.orcacosmo` surface files
6. **Solubility** — run openCOSMO-RS against a configurable solvent list, apply saturation correction

All stages are checkpointed and resumable. A background queue worker manages execution. A CLI handles submission and monitoring.

---

## Quick start

See [docs/quick_start.md](docs/quick_start.md) for the full setup guide. In brief:

**1. Place software and resources in your home directory** (not the project folder — long paths cause silent truncation in XTB/gXTB):

```
~/software/orca_6_1_1/orca
~/software/xtb/bin/xtb
~/software/g-xtb-main/binary/gxtb
~/software/crest/crest
~/resources/openCOSMO-RS_py/src/
~/resources/openCOSMO-RS_cpp/bindings/
```

**2. Install Python dependencies:**

```bash
pip install -r requirements.txt
```

**3. Add the CLI alias to `~/.bashrc`:**

```bash
alias pl='python3 -m modules.cli.cli'
```

**4. Update `config/paths.json`** with your software paths.

**5. Validate your environment:**

```bash
pl env check
```

**6. Start the worker:**

```bash
pl q start
```

---

Edit `modules/main.py` to configure your pipeline and input data, then run:

```bash
python3 -m modules.main
```

The pipeline spec is a list of stage dicts — each entry has a `stage` name and an `args` dict of parameters for that stage:

```python
input_csv = "/path/to/molecules.csv"

pipeline_spec = [
    {"stage": "cleaning",     "args": {"input_csv": input_csv, "overwrite_metadata": True}},
    {"stage": "generation",   "args": {"engine": "rdkit", "n": 5}},
    {"stage": "pruning",      "args": {"n": 1}},
    {"stage": "optimisation", "args": {"engine": "gxtb_opt_normal"}},
    {"stage": "orcacosmo",    "args": {}},
    {"stage": "solubility",   "args": {"solvent_list":"default_list"}},
]

parameters = {
    "title": "my_run",
    "resources": {"cpus": 20, "memory_gb": 64},
    "config": config,
}
```

You can repeat stages (e.g. multiple sequential optimisation passes with different engines). Set `USE_QUEUE = True` in `main.py` to enqueue the run for the background worker, or `False` to run directly in the foreground.

---

## Common commands

```bash
pl q start                  # Start the queue worker
pl q stop                   # Graceful stop
pl q status                 # Worker state and queue counts

pl r list                   # List all requests
pl r status <id>            # Pipeline state for a request
pl r logs <id>              # Stage logs
pl r info <id>              # Detailed request info

pl env check                # Full environment validation
pl env software             # Validate executables only
pl env resources            # Validate openCOSMO-RS paths
```

---

## Outputs

When a full pipeline completes, results are in:

```
pipeline_data/requests/<id>/final_outputs/
    solubility_results.json          machine-readable results
    solubility_results.csv           flat per-molecule/solvent table
    solubility_human_summary.txt     aligned human-readable table
    molecule_metadata/               per-molecule JSON descriptors
    pipeline_warnings.txt            all warnings from the run
```

---

## Data directory layout

```
pipeline_data/requests/<request_id>/
    request.json              submitted pipeline spec
    request_state.json        pipeline state (stage, status)
    request.log               high-level request log
    jobs/
        J-<timestamp>-<stage>/
            inputs/           canonical input from previous stage
            outputs/          canonical output + checkpoints
            stage_logs/
                stage.log     full stage execution log
                stage_context.log
            job_state.json    item-level progress
    final_outputs/            collected results and warnings
```

---

## Documentation

| Document | Contents |
|----------|----------|
| [docs/quick_start.md](docs/quick_start.md) | Setup, installation, first run |
| [docs/pipeline_overview.md](docs/pipeline_overview.md) | Scientific workflow, execution model, design principles |
| [docs/arguments_guide.md](docs/arguments_guide.md) | All pipeline parameters and stage arguments |
| [docs/input_csv_format.md](docs/input_csv_format.md) | Input CSV column reference |
| [docs/output_format.md](docs/output_format.md) | Output file schemas |
| [docs/solvent_list_format.md](docs/solvent_list_format.md) | Solvent list format and available solvents |
| [docs/developer_guide.md](docs/developer_guide.md) | Architecture, adding backends/stages, BaseStage contract |
| [docs/troubleshooting.md](docs/troubleshooting.md) | Common failures and diagnostics |

All documents are available as `.md` and `.docx`.

---

## Authors

Luke Hetherington
