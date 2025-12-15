import os
import pandas as pd
import csv

from .data_cleaning import DataCleaner
from .conformer_generation import ConformerGenerator
from .conformer_pruning import ConformerPruner
from .conformer_optimisation import ConformerOptimiser
from .cosmo_file_generation import CosmoFileGenerator

from .resource_management import ResourceManager, Benchmark
from .molecule_utils import MoleculeUtils

OPTIMISATION_METHODS = {
    "xtb": {"runner": "_run_xtb"},
    "gxtb": {"runner": "_run_gxtb"},
    "forcefield": {"runner": "_run_forcefield"},
    "dft": {"runner": "_run_dft", "mode": "standard"},
    "dft_fast": {"runner": "_run_dft", "mode": "fast"},
    "dft_final": {"runner": "_run_dft", "mode": "final"},
    "dft_cpcm": {"runner": "_run_dft", "mode": "cpcm"},
}

class WorkflowManager:
    def __init__(self, base_dir="pipeline_data"):
        self.base_dir = base_dir
        self.folders = {
            "raw_data": os.path.join(base_dir, "1_raw_data"),
            "clean_data": os.path.join(base_dir, "2_clean_data"),
            "conformer_xyz": os.path.join(base_dir, "3_conformer_xyz"),
            "pruned_conformers": os.path.join(base_dir, "4_pruned_conformers"),
            "optimised_conformers": os.path.join(base_dir, "5_conformer_xyz_optimised"),
            "cosmo_files": os.path.join(base_dir, "6_cosmo_files"),
            "reports": os.path.join(base_dir, "7_reports"),
            "temp_exec": os.path.join(base_dir, "tmp_exec"),
        }
        for f in self.folders.values():
            os.makedirs(f, exist_ok=True)

        self.last_outputs = {}

    # ------------------ All-in-one Pipeline ------------------
    def run_pipeline(self, stages):
        """
        Run pipeline stages in the order provided.
        Each element in `stages` is a dict with keys:
        - "stage": string name of the stage
        - "args": dictionary of arguments for that stage
        """
        results = {}
        for step in stages:
            stage = step["stage"]
            args = step.get("args")

            if stage == "cleaning" and args:
                results["clean_csv"] = self.run_data_cleaning(args)

            elif stage == "generation" and args:
                results["conformer_xyz"] = self.run_conformer_generation(args)

            elif stage == "pruning" and args:
                results["pruned_conformers"] = self.run_conformer_pruning(args)

            elif stage == "optimisation" and args:
                results["optimised_conformers"] = self.run_conformer_optimisation(args)

            elif stage == "cosmo" and args:
                results["cosmo_files"] = self.run_orca_cosmo_step(args)

            else:
                raise ValueError(f"Unknown stage: {stage}")

        return results


    # ------------------ Cleaning ------------------
    def run_data_cleaning(self, cleaning_args=None):
        if cleaning_args is None or "dataset" not in cleaning_args:
            raise ValueError("cleaning_args must include 'dataset' (dataset name).")

        dataset = cleaning_args["dataset"]
        input_subdir = cleaning_args.get("input_subdir")

        # If input_subdir not provided, auto-detect CSV in raw_data/<dataset>
        if input_subdir is None:
            dataset_dir = os.path.join(self.folders["raw_data"], dataset)
            if not os.path.isdir(dataset_dir):
                raise FileNotFoundError(f"Raw dataset folder not found: {dataset_dir}")

            csv_files = [f for f in os.listdir(dataset_dir) if f.endswith(".csv")]
            if len(csv_files) == 0:
                raise FileNotFoundError(f"No CSV files found in {dataset_dir}")
            elif len(csv_files) > 1:
                raise ValueError(f"Multiple CSV files found in {dataset_dir}, please specify 'input_subdir'")
            else:
                input_subdir = os.path.join(dataset, csv_files[0])

        raw_path = os.path.join(self.folders["raw_data"], input_subdir)
        if not os.path.exists(raw_path):
            raise FileNotFoundError(f"Raw input file not found: {raw_path}")

        # Build cleaned filename
        cleaned_filename = cleaning_args.get("cleaned_filename")
        if cleaned_filename is None:
            base, ext = os.path.splitext(os.path.basename(raw_path))
            cleaned_filename = f"{dataset}/{base}_clean.csv"

        clean_path = os.path.join(self.folders["clean_data"], cleaned_filename)
        os.makedirs(os.path.dirname(clean_path), exist_ok=True)

        cleaner = DataCleaner(raw_path, clean_path)
        clean_csv = cleaner.clean()  # returns path to cleaned CSV

        return {
            "dataset": dataset,
            "raw_csv": raw_path,
            "clean_csv": clean_csv,
        }



    # ------------------ Generation ------------------
    def run_conformer_generation(self, generation_args: dict):
        """
        Run conformer generation stage.

        Required arguments:
        - dataset (str): name of dataset
        - engine (str): conformer generation engine ('rdkit', 'crest', 'openbabel')

        Optional arguments:
        - clean_csv (str): relative path to cleaned CSV; if not provided, inferred from dataset
        - n_conformers (int): number of conformers to generate; defaults to 100
        - seed (int): random seed; defaults to 42

        Behaviour:
        - Generates conformers into a dataset/engine folder.
        - Creates a master CSV recording provenance.
        - Stores master CSV path in self.last_outputs for downstream auto‑detection.
        """

        dataset = generation_args.get("dataset")
        if not dataset:
            raise ValueError("Missing required argument: 'dataset'")

        engine = generation_args.get("engine")
        if not engine:
            raise ValueError("Missing required argument: 'engine'")
        method = engine.lower()
        valid_methods = {"rdkit", "crest", "openbabel"}
        if method not in valid_methods:
            raise ValueError(f"Invalid engine '{engine}'. Must be one of {', '.join(valid_methods)}.")

        clean_csv_rel = generation_args.get("clean_csv")
        if clean_csv_rel:
            clean_csv_path = os.path.join(self.folders["clean_data"], clean_csv_rel)
        else:
            clean_csv_path = os.path.join(self.folders["clean_data"], f"{dataset}/{dataset}_clean.csv")
        if not os.path.exists(clean_csv_path):
            raise FileNotFoundError(f"Cleaned CSV not found: {clean_csv_path}")

        out_dir_rel = os.path.join(dataset, method)
        out_dir = os.path.join(self.folders["conformer_xyz"], out_dir_rel)
        os.makedirs(out_dir, exist_ok=True)

        generator = ConformerGenerator(dataset=dataset, method=method, output_dir=out_dir)

        n_confs = generation_args.get("n_conformers", 100)
        seed = generation_args.get("seed", 42)

        df = pd.read_csv(clean_csv_path)
        for _, row in df.iterrows():
            generator.generate(
                smiles=row["smiles"],
                charge=int(row.get("charge", 0)),
                num_confs=n_confs,
                seed=seed
            )

        master_csv_path = os.path.join(out_dir, f"{dataset}_{method}_{n_confs}_master.csv")

        # --- Remember master CSV for downstream stages ---
        self.last_outputs["master_csv"] = master_csv_path

        return {
            "dataset": dataset,
            "engine": method,
            "n_conformers": n_confs,
            "conformer_xyz_folder": out_dir,
            "master_csv": master_csv_path
        }


    # ------------------ Pruning ------------------
    def run_conformer_pruning(self, pruning_args: dict):
        """
        Run conformer pruning stage.

        Required arguments:
        - dataset (str): name of dataset (e.g. 'acr_t5')
        - engine (str): conformer generation engine ('rdkit', 'crest', 'openbabel')

        Optional arguments:
        - method (str): pruning method; defaults to 'topN'
        - master_csv (str): explicit path to master CSV; if not provided, auto-located
        - lookup_csv (str): optional existing lookup CSV to restrict pruning subset
        - top_n (int): number of conformers to keep if using 'topN'
        - energy_window (float): kcal/mol window if using 'energy_window'
        - lower_pct / upper_pct (float): bounds if using 'percentile'
        - any other method-specific parameters

        Behaviour:
        - Resolves master CSV using find_master_csv (explicit arg → remembered → auto-search).
        - Creates a new pruning output folder under
            `pipeline_data/4_pruned_conformers/<dataset>/<engine>/`.
        - Writes a lookup CSV recording pruned conformers.
        - Returns paths and provenance in a dictionary.
        """
        dataset = pruning_args["dataset"]
        engine = pruning_args["engine"]
        method = pruning_args.get("method", "topN")

        # --- Master CSV resolution ---
        master_csv = self.find_master_csv(dataset, engine, pruning_args.get("master_csv"))

        # --- Ensure base pruning folder exists ---
        pruned_base = os.path.join(self.folders["pruned_conformers"], dataset, engine)
        os.makedirs(pruned_base, exist_ok=True)

        # --- Create output folder for pruned conformers ---
        step_index = self._next_step_index(pruned_base, "lookup")

        # Build a suffix that includes method + key parameter
        if method == "topN":
            param_str = f"topN_{pruning_args.get('top_n', 10)}"
        elif method == "energy_window":
            param_str = f"energyWindow_{pruning_args.get('energy_window', 5.0)}"
        elif method == "percentile":
            lower = pruning_args.get("lower_pct", 0)
            upper = pruning_args.get("upper_pct", 100)
            param_str = f"percentile_{lower}-{upper}"
        else:
            param_str = method  # fallback

        suffix = f"{param_str}_step{step_index}"
        out_dir = os.path.join(pruned_base, suffix)
        os.makedirs(out_dir, exist_ok=True)

        # --- Build lookup filename ---
        lookup_filename = f"lookup_{dataset}_{engine}_{suffix}.csv"
        lookup_path = os.path.join(out_dir, lookup_filename)

        # --- Run pruning ---
        params = {k: v for k, v in pruning_args.items()
                if k not in ("dataset", "engine", "method", "master_csv")}
        
        pruner = ConformerPruner(
            master_csv=master_csv,
            output_dir=out_dir,
            lookup_csv=pruning_args.get("lookup_csv"),
            energy_source=pruning_args.get("energy_source")
        )

        lookup_path, n_selected, energy_source = pruner.prune(
            method=method,
            output_file=lookup_path,
            **params
        )

        # --- Record outputs ---
        self.last_outputs["pruned_conformers"] = os.path.relpath(out_dir, self.folders["pruned_conformers"])
        self.last_outputs["lookup_csv"] = lookup_path
        self.last_outputs["master_csv"] = master_csv

        return {
            "dataset": dataset,
            "engine": engine,
            "method": method,
            "lookup_csv": lookup_path,
            "pruned_conformers_folder": out_dir,
            "master_csv": master_csv,
            "energy_source": energy_source,
            "n_selected": n_selected
        }



    def run_conformer_optimisation(self, opt_args: dict):
        dataset = opt_args["dataset"]
        engine = opt_args["engine"]

        # --- Resolve lookup CSV ---
        lookup_csv = opt_args.get("lookup_csv") or self.last_outputs.get("lookup_csv")

        # --- Resolve input XYZ directory ---
        xyz_dir = opt_args.get("xyz_dir")
        if not xyz_dir:
            if self.last_outputs.get("optimised_conformers"):
                xyz_dir = os.path.join(self.folders["optimised_conformers"],
                                    self.last_outputs["optimised_conformers"])
            else:
                genmethod = opt_args.get("conf_prod", engine)
                xyz_dir = os.path.join(self.folders["conformer_xyz"], dataset, genmethod)

        # --- Output folder ---
        opt_base = os.path.join(self.folders["optimised_conformers"], dataset, engine)
        step_index = self._next_step_index(opt_base, "optimised")
        out_dir = os.path.join(opt_base, f"step{step_index}")
        os.makedirs(out_dir, exist_ok=True)

        # --- Instantiate optimiser ---
        optimiser = ConformerOptimiser(
            input_xyz_dir=xyz_dir,
            lookup_csv=lookup_csv,
            output_dir=out_dir,
            tmp_exec=self.folders["temp_exec"],
            engine=engine,
            master_csv=self.last_outputs.get("master_csv"),
            **{k: v for k, v in opt_args.items()
            if k not in ("dataset", "engine", "xyz_dir", "lookup_csv", "tmp_exec")}
        )

        # --- Run optimisation ---
        summary_path = optimiser.optimise()

        # --- Record outputs ---
        self.last_outputs.update({
            "optimised_conformers": os.path.relpath(out_dir, self.folders["optimised_conformers"]),
            "lookup_csv": lookup_csv,
            "xyz_dir": xyz_dir,
            "optimisation_summary": summary_path
        })

        return {
            "dataset": dataset,
            "engine": engine,
            "lookup_csv": lookup_csv,
            "xyz_dir": xyz_dir,
            "optimised_conformers_folder": out_dir,
            "summary_csv": summary_path
        }



    # ------------------ COSMO ------------------
    def run_orca_cosmo_step(self, cosmo_args=None):
        """
        Generate COSMO files using ORCA.

        Parameters
        ----------
        cosmo_args : dict
            Dictionary of arguments for COSMO stage.
            Keys:
                - "input_subdir" : relative path to folder of optimised conformers
                - "lookup_csv"   : optional pruning lookup CSV
                - "method"       : functional (default "B3LYP")
                - "basis"        : basis set (default "def2-SVP")
                - "solvent"      : solvent model (default "Water")
                - "charge"       : molecular charge (default 0)
                - "multiplicity" : spin multiplicity (default 1)

        Returns
        -------
        results : dict
            Dictionary keyed by InChIKey with values:
                { "dir": output directory, "n_confs": number of conformers processed }
        """
        if cosmo_args is None:
            raise ValueError("cosmo_args must be provided.")

        # Resolve input folder
        input_subdir = cosmo_args.get("input_subdir", self.last_outputs.get("optimised_conformers"))
        if input_subdir is None:
            raise ValueError("No input_subdir provided for COSMO stage.")

        optimisation_dir = os.path.join(self.folders["optimised_conformers"], input_subdir)
        if not os.path.exists(optimisation_dir):
            raise FileNotFoundError(f"Optimisation directory not found: {optimisation_dir}")

        # Filter xyz files via MoleculeUtils
        lookup_csv = cosmo_args.get("lookup_csv")
        xyz_files = self.filter_xyz_with_lookup(optimisation_dir, lookup_csv)

        # Defaults
        method = cosmo_args.get("method", "B3LYP")
        basis = cosmo_args.get("basis", "def2-SVP")
        solvent = cosmo_args.get("solvent", "Water")
        charge = cosmo_args.get("charge", 0)
        multiplicity = cosmo_args.get("multiplicity", 1)

        generator = CosmoFileGenerator(cosmo_root=self.folders["cosmo_files"])
        results = {}

        # Loop handled here in WFM
        for xyz_file in xyz_files:
            try:
                orcacosmo_path = generator.run_cosmo_for_xyz(
                    xyz_file,
                    method=method,
                    basis=basis,
                    solvent=solvent,
                    charge=charge,
                    multiplicity=multiplicity
                )
                print(f"Generated {orcacosmo_path}")

                lookup_id = os.path.splitext(os.path.basename(xyz_file))[0]
                inchi = lookup_id.split("_conf")[0]
                mol_dir = os.path.join(self.folders["cosmo_files"], inchi)

                results.setdefault(inchi, {"dir": mol_dir, "n_confs": 0})
                results[inchi]["n_confs"] += 1

            except RuntimeError as e:
                print(f"Skipping {xyz_file}: {e}")

        if not results:
            raise RuntimeError("No COSMO files generated; all ORCA jobs failed.")

        self.last_outputs["cosmo_files"] = self.folders["cosmo_files"]
        return results

    def _next_step_index(self, folder: str, prefix: str) -> int:
        """
        Determine the next step index for files in a folder.
        If the folder does not exist yet, return 1.
        """
        if not os.path.exists(folder):
            return 1

        existing = [
            f for f in os.listdir(folder)
            if f.startswith(prefix) and f.endswith(".csv")
        ]
        if not existing:
            return 1

        # Extract step numbers from filenames like lookup_topN_step3.csv
        indices = []
        for f in existing:
            parts = f.replace(".csv", "").split("_")
            for p in parts:
                if p.startswith("step"):
                    try:
                        indices.append(int(p.replace("step", "")))
                    except ValueError:
                        pass
        return max(indices) + 1 if indices else 1


    def filter_xyz_with_lookup(self, optimisation_dir, lookup_csv=None):
        """
        Collect .xyz files from optimisation_dir, optionally filter by lookup CSV.

        Parameters
        ----------
        optimisation_dir : str
            Path to folder containing optimised conformers.
        lookup_csv : str, optional
            Path to pruning lookup CSV. If provided, only conformers listed here are processed.

        Returns
        -------
        xyz_files : list of str
            List of file paths to .xyz files.
        """
        xyz_files = []
        for root, _, files in os.walk(optimisation_dir):
            for f in files:
                if f.endswith(".xyz"):
                    xyz_files.append(os.path.join(root, f))

        if not xyz_files:
            raise RuntimeError(f"No XYZ files found in {optimisation_dir}")

        if lookup_csv:
            if not os.path.exists(lookup_csv):
                raise FileNotFoundError(f"Lookup CSV not found: {lookup_csv}")
            allowed_ids = set()
            with open(lookup_csv, newline="") as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    allowed_ids.add(row.get("lookup_id") or row.get("filename"))
            xyz_files = [f for f in xyz_files if os.path.splitext(os.path.basename(f))[0] in allowed_ids]

        if not xyz_files:
            raise RuntimeError(f"No XYZ files left to process after filtering with lookup {lookup_csv}")

        return xyz_files
    


    def find_master_csv(self, dataset: str, engine: str, explicit_path: str = None) -> str:
        """
        Locate the master CSV for a given dataset/engine.

        Resolution order:
        1. If explicit_path is provided, validate and return it.
        2. If self.last_outputs contains a master_csv, validate and return it.
        3. Otherwise, search in the conformer folder for files starting with '_master'.

        Args:
            dataset (str): dataset name (e.g. 'acr_t5')
            engine (str): conformer generation engine ('rdkit', 'crest', 'openbabel')
            explicit_path (str, optional): direct path to master CSV

        Returns:
            str: absolute path to master CSV

        Raises:
            FileNotFoundError: if no valid master CSV can be found
            ValueError: if dataset or engine are missing
        """

        if not dataset:
            raise ValueError("Missing required argument: 'dataset'")
        if not engine:
            raise ValueError("Missing required argument: 'engine'")

        # 1. Explicit path
        if explicit_path:
            if os.path.exists(explicit_path):
                return explicit_path
            raise FileNotFoundError(f"Specified master CSV not found: {explicit_path}")

        # 2. Remembered from last_outputs
        remembered = self.last_outputs.get("master_csv")
        if remembered and os.path.exists(remembered):
            return remembered

        # 3. Auto-search in conformer folder
        base_dir = os.path.join(self.folders["conformer_xyz"], dataset, engine)
        if not os.path.exists(base_dir):
            raise FileNotFoundError(f"Conformer folder not found: {base_dir}")

        master_csv = next(
            (os.path.join(base_dir, f) for f in os.listdir(base_dir) if f.startswith("_master")),
            None
        )
        if master_csv is None:
            raise FileNotFoundError(f"No master CSV found in {base_dir}")

        return master_csv



    def find_lookup_csv(self, dataset: str, conf_prod: str) -> str:
        """
        Locate the most recent pruning lookup CSV for a given dataset/engine.
        Strategy: search pruned_conformers/<dataset>/<conf_prod>/,
        pick the lookup file from the pruning step with the fewest conformers.
        """
        pruned_base = os.path.join(self.folders["pruned_conformers"], dataset, conf_prod)
        if not os.path.exists(pruned_base):
            raise FileNotFoundError(f"No pruned conformers found in {pruned_base}")

        candidate_dirs = [
            os.path.join(pruned_base, d)
            for d in os.listdir(pruned_base)
            if os.path.isdir(os.path.join(pruned_base, d))
        ]
        if not candidate_dirs:
            raise FileNotFoundError(f"No pruning steps found under {pruned_base}")

        lookup_csvs = []
        for d in candidate_dirs:
            for f in os.listdir(d):
                if f.startswith("lookup_") and f.endswith(".csv"):
                    path = os.path.join(d, f)
                    # count rows to decide "most recent" (fewest conformers)
                    nrows = sum(1 for _ in open(path)) - 1  # minus header
                    lookup_csvs.append((nrows, path))

        if not lookup_csvs:
            raise FileNotFoundError(f"No lookup CSVs found under {pruned_base}")

        # select the one with the least conformers
        _, lookup_csv = min(lookup_csvs, key=lambda x: x[0])
        return lookup_csv

    def run_roundrobin_optimisation(self,
                                    lookup_csv: str,
                                    input_xyz_dir: str,
                                    opt_engine: str = "dft",
                                    top_n: int = 5,
                                    **kwargs):
        """
        Run conformer optimisations in round-robin order across molecules.

        Args:
            lookup_csv (str): Path to pruned lookup CSV.
            input_xyz_dir (str): Directory containing conformer .xyz files.
            opt_engine (str): optimisation backend ('dft', 'dft_fast', etc.).
            top_n (int): number of conformers per molecule to optimise.
            kwargs: engine-specific parameters (functional, basis, etc.).

        Behaviour:
            - Loads lookup CSV.
            - Groups conformers by molecule ID.
            - Iterates round-robin: first conformer of each molecule, then second, etc.
            - Runs optimiser for each conformer in order.
            - Updates master CSV after each run.
        """
        import pandas as pd, os

        df_lookup = pd.read_csv(lookup_csv)

        # Group by molecule ID (inchi_key)
        grouped = df_lookup.groupby("inchi_key")

        # --- Construct output folder path ---
        dataset = os.path.basename(os.path.dirname(os.path.dirname(lookup_csv)))
        conf_prod = os.path.basename(os.path.dirname(lookup_csv))
        prune_step = os.path.basename(os.path.dirname(lookup_csv))
        out_dir_rel = os.path.join(dataset, conf_prod, prune_step, opt_engine)
        out_dir = os.path.join(self.folders["optimised_conformers"], out_dir_rel)
        os.makedirs(out_dir, exist_ok=True)

        results = []
        for rank in range(top_n):
            for inchi, group in grouped:
                confs = group.sort_values("energy_rank").reset_index(drop=True)
                if rank < len(confs):
                    lookup_id = confs.loc[rank, "lookup_id"]
                    xyz_file = os.path.join(input_xyz_dir, f"{lookup_id}.xyz")

                    optimiser = ConformerOptimiser(
                        input_xyz_dir=input_xyz_dir,
                        lookup_csv=lookup_csv,
                        output_dir=out_dir,
                        tmp_exec=self.folders["temp_exec"],
                        engine=opt_engine,
                        master_csv=self.last_outputs.get("master_csv"),
                        **kwargs
                    )
                    result = optimiser._run_engine(xyz_file)
                    results.append(result)

                    # Update master CSV with energy + method
                    updates = result["updates"]
                    updates["optimisation_method"] = opt_engine
                    self._update_master_csv(updates)

        return {
            "results": results,
            "output_dir": out_dir,
            "lookup_csv": lookup_csv,
            "engine": opt_engine
        }
