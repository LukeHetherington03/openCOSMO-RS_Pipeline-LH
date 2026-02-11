openCOSMO‑RS Pipeline
=====================

A modular, deterministic workflow for generating COSMO‑RS solubility predictions
from SMILES through conformer generation, optimisation, ORCA‑COSMO calculations,
and final COSMO‑RS evaluation.

The pipeline is designed for:
• reproducible computational chemistry workflows  
• large‑scale dataset generation  
• ML‑correctable solubility prediction pipelines  
• transparent provenance and auditability  
• clean stage boundaries with explicit input/output contracts  

-------------------------------------------------------------------------------
Purpose ✨
-------------------------------------------------------------------------------

The openCOSMO‑RS Pipeline provides an end‑to‑end, fully auditable workflow for
turning molecular structures into high‑quality COSMO‑RS property predictions.
It emphasises determinism, reproducibility, and scalability, making it suitable
for both scientific studies and machine‑learning dataset generation.

-------------------------------------------------------------------------------
Features 🚀
-------------------------------------------------------------------------------

• Fully modular stage architecture  
  Each stage is isolated, deterministic, and writes its own canonical inputs
  and outputs. No hidden logic or cross‑stage side effects.

• New execution engine  
  Queue + runner + worker model replaces the old monolithic controller,
  enabling clean job boundaries and robust resumption.

• Atomic writes and strict provenance tracking  
  Every stage records its state, inputs, outputs, and logs in a structured,
  machine‑readable format.

• Clean separation of chemistry constants and configuration  
  CPCM radii, solvent files, and metadata live in CONSTANT_FILES/, while
  executable paths and defaults live in config/.

• ORCA‑COSMO stage rebuilt  
  - deterministic CPCM radii loading  
  - safe ORCA execution environment  
  - TZVPD → TZVP fallback logic  
  - canonical output: orcacosmo_summary.json  

• Solubility stage rebuilt  
  - Python/C++ COSMO‑RS wrapper  
  - per‑molecule directory structure  
  - mixture_inputs.txt generation  
  - raw COSMO‑RS output preserved  
  - solubility_results.json + human‑readable summary  

• High‑throughput quantum chemistry  
  Integrated support for gXTB, ORCA, and openCOSMO‑RS enables workflows from
  fast screening to high‑accuracy DFT‑corrected datasets.

• Robust job management and resumability  
  Each job tracks pending, completed, and failed items. Interrupted runs can
  resume safely without recomputing finished work.

• ML‑ready outputs  
  Structured JSON, canonical identifiers, and full provenance make the results
  ideal for downstream machine‑learning pipelines.

-------------------------------------------------------------------------------
Pipeline Stages (Overview Only) 🔧
-------------------------------------------------------------------------------

1. Cleaning Stage  
   Validates SMILES, canonicalises structures, and prepares identifiers.

2. Generation Stage  
   Produces conformers using gXTB/xtb.

3. Pruning Stage  
   Selects representative conformers and prepares XYZs for ORCA.

4. Optimisation Stage  
   Runs gXTB or ORCA geometry optimisation and records energies.

5. ORCA‑COSMO Stage  
   Runs ORCA CPCM calculations and reconstructs .orcacosmo files.  
   Outputs orcacosmo_summary.json.

6. Solubility Stage  
   Runs COSMO‑RS using Python/C++ bindings and writes solubility_results.json.

7. Cleaning Stage  
   Removes temporary execution directories.

Each stage has its own dedicated .md file with full logic, contracts, and examples.

-------------------------------------------------------------------------------
Running the Pipeline ▶️
-------------------------------------------------------------------------------

From the repository root:

    python3 -m modules.main

A new request directory will be created under:

    pipeline_data/requests/R-<timestamp>/

All stage outputs, logs, and provenance files are stored inside this folder.

A nohup‑safe launcher is also available:

    python3 -m modules.main_nohup

This automatically re‑executes itself inside a nohup environment for unattended
execution.

-------------------------------------------------------------------------------
Configuration ⚙️
-------------------------------------------------------------------------------

Configuration lives in the config/ directory:
• paths.json — ORCA paths, chemistry directories  
• *_defaults.json — stage‑specific defaults  
• resource_allocation.json — CPU allocation per stage  

Constant files (CPCM radii, solvent files, molecule metadata) live in
CONSTANT_FILES/.
Project Structure (High‑Level Overview) 📁
=========================================

The openCOSMO‑RS Pipeline is organised into clear, purpose‑driven modules.
Each folder has a single responsibility, supporting a clean, maintainable,
and auditable workflow from SMILES → COSMO‑RS solubility predictions.

-------------------------------------------------------------------------------
Top‑Level Directories
-------------------------------------------------------------------------------

config/  
    Contains all configuration files used by the pipeline:
    • paths.json — absolute paths to ORCA, gXTB, xtb, COSMO‑RS binaries  
    • *_defaults.json — stage‑specific default parameters  
    • resource_allocation.json — CPU/memory allocation per stage  
    These files define deterministic behaviour for every stage.

CONSTANT_FILES/  
    Chemistry constants and metadata used across the pipeline:
    • chemistry/ — CPCM radii, element parameters  
    • solvents/ — COSMO‑RS solvent files  
    • molecule_metadata/ — curated metadata for molecules  
    These files never change during execution and ensure reproducibility.

docs/  
    Developer and user documentation, including stage contracts, examples,
    and architectural notes.

modules/  
    Core implementation of the pipeline. Each subfolder has a dedicated role:

    build/  
        Request and Job management:
        • request_manager.py — creates and tracks pipeline requests  
        • job_manager.py — handles job creation, state, and resumption  
        • log_helper.py — consistent logging across all stages  
        This layer defines the execution model (items, pending, completed).

    cli/  
        Command‑line entrypoints and helper scripts for running the pipeline
        interactively or in batch mode.

    execution/  
        The execution engine:
        • runner.py — sequential job runner  
        • queue_worker.py — background worker for queued jobs  
        • scheduler tools — optional queue‑based orchestration  
        This layer replaces the old monolithic controller.

    parsers/  
        ORCA and COSMO‑RS output parsers:
        • extract energies, CPCM surfaces, COSMO files  
        • produce structured JSON summaries  
        Ensures clean, machine‑readable outputs for downstream stages.

    post_analysis/  
        Tools for visualisation and benchmarking:
        • plotting utilities  
        • dataset comparison tools  
        • error analysis and ML‑readiness checks  
        Outputs written to post_analysis_results/.

    request_tools/  
        Helpers for inspecting, resuming, or modifying existing requests.

    solubility_engine/  
        Python/C++ bindings for openCOSMO‑RS:
        • mixture_inputs.txt generation  
        • COSMO‑RS execution wrapper  
        • solubility result parsing  
        This is the core of the solubility stage.

    stages/  
        Implementation of each pipeline stage:
        • cleaning_stage.py  
        • generation_stage.py  
        • pruning_stage.py  
        • optimisation_stage.py  
        • orcacosmo_stage.py  
        • solubility_stage.py  
        Each stage is isolated, deterministic, and writes a canonical output.

    utils/  
        Shared utilities:
        • file helpers  
        • chemistry helpers  
        • environment setup  
        • safe directory creation  
        These functions support all stages without introducing hidden logic.

pipeline_data/  
    Automatically generated during execution:
    • requests/ — each request gets its own folder  
    • jobs/ — each stage execution is a job with its own state  
    • logs/ — request‑level logs  
    • raw_outputs/ — ORCA logs, CPCM files, COSMO‑RS raw output  
    • parsed_outputs/ — JSON summaries  
    This directory contains the full provenance of every run.

tests/  
    Unit tests and integration tests for pipeline components.

-------------------------------------------------------------------------------
Overall Module Responsibilities 🧩
-------------------------------------------------------------------------------

• build/ — defines the execution model (Request, Job, state tracking)  
• execution/ — runs jobs in sequence or via a queue worker  
• stages/ — chemistry logic for each pipeline step  
• solubility_engine/ — COSMO‑RS backend (Python/C++ bindings)  
• parsers/ — converts ORCA/COSMO outputs into structured JSON  
• post_analysis/ — visualisation and benchmarking tools  
• utils/ — shared helpers with no side effects  
• config/ + CONSTANT_FILES/ — deterministic configuration and chemistry data  
• pipeline_data/ — all generated outputs, logs, and provenance  

-------------------------------------------------------------------------------
