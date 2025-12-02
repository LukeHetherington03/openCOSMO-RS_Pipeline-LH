# openCOSMO-RS_Pipeline-LH

A modular workflow manager for COSMO-RS property prediction.  
This pipeline automates conformer generation, optimisation, pruning, and COSMO file creation, with a reproducible folder structure for clean data handling and downstream analysis.

---

## 🚀 Features
- **Workflow orchestration** via `main.py` and `work_flow_manager.py`
- **Modular design**: each stage (generation, optimisation, pruning, COSMO file creation) is a separate Python module
- **Data provenance**: structured folder hierarchy (`pipeline_data/`) for raw, cleaned, and processed outputs
- **Resource management**: runtime benchmarking and error handling for scalable workflows
- **Extensible**: modules can be swapped or extended for new computational chemistry tasks

---

## 📂 Repository Structure
openCOSMO-RS_Pipeline-LH/ 
├── modules/ # Modular pipeline components 
│ ├── main.py # Entry point, runs workflow manager 
│ ├── conformer_generation.py 
│ ├── conformer_optimisation.py 
│ ├── conformer_pruning.py 
│ ├── cosmo_file_generation.py 
│ ├── data_cleaning.py 
│ ├── molecule_utils.py 
│ ├── resource_management.py 
│ └── work_flow_manager.py 
├── pipeline_data/ # Staged data folders 
│ ├── 1_raw_data/ 
│ ├── 2_clean_data/ 
│ ├── 3_conformer_xyz/ 
│ ├── 4_pruned_conformers/ 
│ ├── 5_conformer_xyz_optimised/ 
│ ├── 6_cosmo_files/ 
│ └── 7_reports/ 
└── README.md

## 📦 Dependencies
This project requires Python 3.9+ and the following packages:

- **Core scientific stack**
  - `numpy`
  - `scipy`
  - `pandas`
- **Chemistry toolkits**
  - `rdkit` (for molecule handling and conformer generation)
  - `xtb` (via command line, for semiempirical optimisation)
  - `orca` (external quantum chemistry engine, for CPCM/COSMO runs)
- **Workflow utilities**
  - `tqdm` (progress bars)
  - `logging` (Python standard library, structured logs)
  - `argparse` (command-line interface)

> ⚠️ Note: ORCA and xTB must be installed separately and available in your `$PATH`. COSMO-RS analysis requires the openCOSMO package.

---

## 🖥️ Usage
1. Place input molecules in `pipeline_data/1_raw_data/`.
2. Run the workflow manager:
   ```bash
   python main.py

Code