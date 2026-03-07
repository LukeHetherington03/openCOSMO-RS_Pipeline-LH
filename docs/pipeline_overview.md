# Pipeline Overview

openCOSMO-RS Pipeline

---

## Purpose

The openCOSMO-RS Pipeline is a fully automated, reproducible scientific workflow for predicting molecular solubility in arbitrary solvents. It chains together conformer generation, geometry optimisation, quantum-chemical COSMO surface calculations (ORCA), and COSMO-RS thermodynamic modelling (openCOSMO-RS) into a single managed system with checkpointing, resumability, provenance tracking, and a CLI for job management.

---

## Scientific workflow

The pipeline implements the following sequence of scientific steps:

```
Raw molecule CSV(s)
      |
      v
1. Cleaning
   - Validate and canonicalise SMILES
   - Generate InChIKeys
   - Resolve charge, multiplicity, melting point
   - Compute physchem descriptors
   - Write per-molecule metadata
      |
      v  cleaned.csv
2. Generation
   - Generate 3D conformers (RDKit ETKDGv3, CREST, or OpenBabel)
   - One set of conformers per unique molecule
      |
      v  energies.json
3. Pruning  (optional)
   - Filter by energy window and/or RMSD threshold
   - Reduce conformer count before expensive QM steps
      |
      v  energies.json
4. Optimisation
   - Geometry optimise each conformer
   - Engines: gXTB, XTB, ORCA DFT, or force fields (MMFF94/UFF)
   - Each conformer optimised independently in parallel
      |
      v  energies.json
5. ORCA COSMO
   - Run ORCA CPCM single-point calculation per conformer
   - Produces .orcacosmo surface charge density files
   - Fallback basis set if primary fails
      |
      v  orcacosmo_summary.json
6. Solubility
   - Run openCOSMO-RS for each molecule against each solvent combo
   - Apply saturation correction using melting point (Myrdal-Yalkowsky)
   - Output log10 solubility in mole fraction
      |
      v  solubility_results.json / .csv / _human_summary.txt
```

Not all stages are required in every run. The pipeline supports partial sequences — for example running only cleaning and generation, or starting from an existing `orcacosmo_summary.json` and running only solubility.

---

## Execution model

The pipeline is driven by a **queue worker** — a background process that picks up submitted requests and executes them one at a time.

```
User submits request.json  →  pl r submit request.json
                                    |
                                    v
                           Queue (SQLite)
                                    |
                            QueueWorker polls
                                    |
                                    v
                           PipelineRunner.run()
                                    |
                     ┌──────────────┼──────────────┐
                     v              v              v
               Job (cleaning)  Job (gen)  ...  Job (solubility)
                     |
               Stage.execute()
                     |
               ProcessPoolExecutor
                     |
           [worker] [worker] [worker] ...
```

Each stage runs its items (molecules or conformers) in parallel using a `ProcessPoolExecutor`. The number of workers and cores per item are allocated automatically based on the node's CPU budget, with overrides available in the request.

---

## Request and job structure

A **Request** is the top-level container for a pipeline run. It holds:
- The submitted `request.json` (pipeline spec, stage arguments)
- A state file tracking which stages have completed
- A high-level `request.log`
- One **Job** per stage

A **Job** owns a single stage execution:
- `inputs/` — the canonical output from the previous stage, plus resolved parameters
- `outputs/` — the canonical output of this stage, plus checkpoints and auxiliary files
- `stage_logs/stage.log` — full execution log for this stage
- `job_state.json` — item-level progress (pending / completed / failed counts)

```
pipeline_data/requests/<request_id>/
    request.json
    request_state.json
    request.log
    jobs/
        J-<timestamp>-cleaning/
            inputs/
            outputs/
                cleaned.csv
                molecule_metadata/
                checkpoints/
            stage_logs/
                stage.log
                stage_context.log
        J-<timestamp>-generation/
            ...
        ...
    final_outputs/
        solubility_results.json
        solubility_results.csv
        solubility_human_summary.txt
        molecule_metadata/
        pipeline_warnings.txt
```

---

## Key design principles

### Canonical inputs and outputs

Every stage declares exactly one **canonical output file**. This file is copied into the next job's `inputs/` directory as `stage_input`. There is no path guessing — the handoff is deterministic and recorded in `job_state.json`.

| Stage | Canonical output |
|-------|----------------|
| cleaning | `cleaned.csv` |
| generation | `energies.json` |
| pruning | `energies.json` |
| optimisation | `energies.json` |
| orcacosmo | `orcacosmo_summary.json` |
| solubility | `solubility_results.json` |

### Checkpointing and resumability

Every parallel worker writes an atomic checkpoint after completing its item. If the pipeline is interrupted (worker stopped, node crash, HPC preemption), resubmitting the request resumes from the last completed stage. Within that stage, items with existing checkpoints are skipped and logged as `[RESUME]`.

### Provenance tracking

The `energies.json` ConformerSet format accumulates an `optimisation_history` list as conformers pass through generation, pruning, and optimisation stages. Each entry records the engine used, level, energy, geometry path, timing, tool version, and timestamp. The final `solubility_results.json` embeds this full history per conformer.

Molecule metadata (charge, multiplicity, melting point, physchem descriptors, functional groups) is written per-molecule at the cleaning stage and flows through the pipeline unchanged.

### Warning transparency

All non-fatal decisions that could affect results are emitted as `[WARNING]` log entries:
- Default charge or multiplicity applied (with per-file breakdown of user_input / rdkit_computed / default counts)
- Missing melting point (affects saturation correction)
- Basis set fallback triggered in ORCA COSMO
- Missing solvent list (using fallback)
- Resource budget reduced

At the end of every run, all warnings are collected, deduplicated, and written to `final_outputs/pipeline_warnings.txt` for review.

### Config-driven, no hard-coded paths

All software paths, resource limits, and stage defaults are in `config/`. No path is hard-coded in Python. This allows the pipeline to run on different machines by updating `config/paths.json` only.

### Strict mode

Each stage supports a `strict` flag. When enabled, any per-item failure aborts the entire pool immediately rather than continuing to the next item. Disabled by default — the pipeline is designed to produce partial results and report failures rather than halt on any single molecule.

---

## Parallelism and resource allocation

Before each stage, the **ResourceAllocator** computes the optimal worker configuration:

```
max_cores     (from config/resource_allocation.json or request override)
cores_per_item (from stage defaults or request override)
n_workers     = floor(max_cores / cores_per_item), capped at n_items
```

For ORCA and XTB/gXTB, `cores_per_item` is passed directly to the tool as `%pal nprocs` or `OMP_NUM_THREADS` respectively. For RDKit and force-field backends, parallelism is at the molecule level only (serial per molecule, parallel across molecules).

---

## Solvent model

Solvents are defined by pre-computed COSMO files stored in `CONSTANT_FILES/solvents/<name>/`. A **solvent combo** is a pure solvent or mixture of solvents at specified mole fractions, defined in a solvent list JSON file.

The solubility calculation uses openCOSMO-RS (Python + C++ bindings) to compute the activity coefficient of the solute in the solvent mixture, then applies a saturation correction using the Myrdal-Yalkowsky fusion model if a melting point is available.

---

## CLI reference summary

```bash
# Start/stop the worker
pl q start
pl q stop
pl q status

# Submit and monitor requests
pl r submit request.json
pl r list
pl r status <id>
pl r logs <id>

# Environment validation
pl env check
```

Full documentation: see [quick_start.md](quick_start.md) for setup, [arguments_guide.md](arguments_guide.md) for all parameters, [developer_guide.md](developer_guide.md) for extending the pipeline.
