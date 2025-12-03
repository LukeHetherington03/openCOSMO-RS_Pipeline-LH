import os
import pandas as pd

from .data_cleaning import DataCleaner
from .conformer_generation import ConformerGenerator
from .conformer_pruning import ConformerPruner
from .conformer_optimisation import ConformerOptimiser
from .cosmo_file_generation import CosmoFileGenerator
from .resource_management import ResourceManager, Benchmark


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

        # Track last outputs for chaining
        self.last_outputs = {}


    # ------------------ All-in-one Pipeline ------------------
    def run_pipeline(self,
                    raw_relpath=None,              # str: raw input file under 1_raw_data/
                    generation_statement=None,     # dict: {engine, n_conformers, seed}
                    pruning_args=None,             # dict: {method_subdir?, method, params}
                    optimisation_args=None,        # dict: {pruned_subdir?, engine, params}
                    cosmo_args=None):              # dict: {optimisation_subdir?, method, basis, solvent, charge, multiplicity}
        ...

        """
        Run the full pipeline. Each stage uses provided arguments or defaults to last output.
        """
        # Step 1: Cleaning
        if raw_relpath:
            self.run_data_cleaning(raw_relpath)

        # Step 2: Generation
        if generation_statement:
            self.run_conformer_generation(generation_statement=generation_statement)

        # Step 3: Pruning
        if pruning_args:
            self.run_conformer_pruning(**pruning_args)

        # Step 4: Optimisation
        if optimisation_args:
            self.run_conformer_optimisation(**optimisation_args)

        # Step 5: COSMO
        if cosmo_args:
            return self.run_orca_cosmo_step(**cosmo_args)

        return self.last_outputs

    # ------------------ Cleaning ------------------
    def run_data_cleaning(self, raw_relpath, cleaned_filename=None):
        raw_path = os.path.join(self.folders["raw_data"], raw_relpath)
        if cleaned_filename is None:
            base, ext = os.path.splitext(raw_relpath)
            cleaned_filename = base + "_clean" + ext
        clean_path = os.path.join(self.folders["clean_data"], cleaned_filename)
        os.makedirs(os.path.dirname(clean_path), exist_ok=True)

        cleaner = DataCleaner(raw_path, clean_path)
        inp_file = cleaner.clean()
        self.last_outputs["clean_data"] = clean_path
        return inp_file

    # ------------------ Generation ------------------
    def run_conformer_generation(self, inp_relpath=None, generation_statement=None):
        if inp_relpath is None:
            inp_relpath = os.path.relpath(self.last_outputs.get("clean_data"), self.folders["clean_data"])
        inp_path = os.path.join(self.folders["clean_data"], inp_relpath)

        method = generation_statement["engine"].lower()
        out_dir = os.path.join(self.folders["conformer_xyz"], os.path.dirname(inp_relpath), method)
        os.makedirs(out_dir, exist_ok=True)

        generator = ConformerGenerator(method=method, output_dir=out_dir)
        with open(inp_path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                inchi_key, smiles, _, charge = parts
                generator.generate(
                    smiles=smiles,
                    charge=int(charge),
                    num_confs=generation_statement.get("n_conformers"),
                    seed=generation_statement.get("seed", 42)
                )

        self.last_outputs["conformer_xyz"] = out_dir
        return out_dir

    # ------------------ Pruning ------------------
    def run_conformer_pruning(self, method_subdir=None, method="energy_window", **params):
        if method_subdir is None:
            method_subdir = os.path.relpath(self.last_outputs.get("conformer_xyz"), self.folders["conformer_xyz"])
        dataset, engine = method_subdir.split("/")

        # suffix naming
        if method == "topN":
            suffix = f"{method}_{params.get('top_n', 10)}"
        elif method == "energy_window":
            suffix = f"{method}_{params.get('energy_window', 5.0)}"
        elif method == "percentile":
            suffix = f"{method}_{params.get('lower_pct', 0)}_{params.get('upper_pct', 100)}"
        elif method == "rot_bond":
            suffix = f"{method}_{params.get('max_rot_bonds', 8)}"
        elif method == "rmsd":
            suffix = f"{method}_{params.get('rmsd_threshold', 0.5)}"
        else:
            suffix = method

        energies_csv = os.path.join(self.folders["conformer_xyz"], method_subdir, "_conformer_energies.csv")
        out_dir = os.path.join(self.folders["pruned_conformers"], dataset, engine, suffix)
        os.makedirs(out_dir, exist_ok=True)

        lookup_filename = f"lookup_{dataset}_{engine}_{suffix}.csv"
        lookup_path = os.path.join(out_dir, lookup_filename)

        pruner = ConformerPruner(energies_csv, out_dir)
        result = pruner.prune(method=method, output_file=lookup_path, **params)

        self.last_outputs["pruned_conformers"] = out_dir
        return result

    # ------------------ Optimisation ------------------
    def run_conformer_optimisation(self, pruned_subdir=None, engine="gxtb", **params):
        if pruned_subdir is None:
            pruned_subdir = os.path.relpath(self.last_outputs.get("pruned_conformers"), self.folders["pruned_conformers"])
        dataset, eng, method_suffix = pruned_subdir.split("/")

        lookup_filename = f"lookup_{dataset}_{eng}_{method_suffix}.csv"
        lookup_csv = os.path.join(self.folders["pruned_conformers"], pruned_subdir, lookup_filename)

        if engine == "gxtb":
            suffix = f"{engine}_{params.get('level','GFN2')}"
        elif engine == "orca":
            suffix = f"{engine}_{params.get('functional','B3LYP')}_{params.get('basis','def2-SVP')}"
        else:
            suffix = engine

        out_dir = os.path.join(self.folders["optimised_conformers"], pruned_subdir, suffix)
        os.makedirs(out_dir, exist_ok=True)

        optimiser = ConformerOptimiser(
            lookup_csv,
            out_dir,
            engine=engine,
            dataset=dataset,
            generation_engine=eng,
            pruning_method=method_suffix,
            **params
        )
        result = optimiser.optimise()

        self.last_outputs["optimised_conformers"] = out_dir
        return result

    # ------------------ COSMO ------------------
    def run_orca_cosmo_step(self, optimisation_subdir=None,
                            method="B3LYP", basis="def2-SVP", solvent="Water",
                            charge=0, multiplicity=1):
        if optimisation_subdir is None:
            optimisation_subdir = os.path.relpath(self.last_outputs.get("optimised_conformers"), self.folders["optimised_conformers"])
        optimisation_dir = os.path.join(self.folders["optimised_conformers"], optimisation_subdir)

        generator = CosmoFileGenerator(cosmo_root=self.folders["cosmo_files"])
        results = generator.generate_orca_cosmo_from_folder(
            optimisation_dir,
            method=method,
            basis=basis,
            solvent=solvent,
            charge=charge,
            multiplicity=multiplicity
        )
        self.last_outputs["cosmo_files"] = self.folders["cosmo_files"]
        return results
