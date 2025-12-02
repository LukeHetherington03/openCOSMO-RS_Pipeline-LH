"""
conformer_optimisation.py

Module for geometry optimisation of pruned conformers.
Reads a pruned_conformers_lookup_<method>.csv file, locates the corresponding .xyz files,
and performs geometry optimisation using either xTB (GFN-xTB) or ORCA (DFT).

Available engines:
    - gxtb : GFN-xTB geometry optimisation
        Parameters:
            - level : str (default="GFN2-xTB")
            - max_iter : int (default=250)

    - orca : ORCA DFT geometry optimisation
        Parameters:
            - functional : str (default="B3LYP")
            - basis : str (default="def2-SVP")
            - max_iter : int (default=250)

Outputs:
    - Optimised geometries written to 5_conformer_xyz_optimised/<engine>/<method>_<param>/
    - Each file named <inchi_key>_confN_opt.xyz
    - A log CSV summarising optimisation results: energies, convergence status, runtime
"""

import os
import subprocess
import pandas as pd


class ConformerOptimiser:
    def __init__(self,
                 lookup_csv,
                 output_dir,
                 engine="gxtb",
                 dataset=None,
                 generation_engine=None,
                 pruning_method=None,
                 **params):
        """
        Initialise a conformer optimisation job.

        Parameters
        ----------
        lookup_csv : str
            Path to the lookup CSV produced by pruning stage.
            Naming convention: lookup_<dataset>_<engine>_<method_param>.csv
            This file contains the conformer IDs to be optimised.
        output_dir : str
            Directory where optimisation results will be written.
            Subfolders 'xyz/' and 'log/' will be created inside this directory.
        engine : str, default="gxtb"
            Optimisation backend to use. Supported values:
                - "gxtb" : GFN-xTB quantum chemistry engine
                - "orca" : ORCA quantum chemistry engine
        dataset : str, optional
            Name of the dataset (e.g. "acrylates").
            Used for provenance in optimisation summary.
        generation_engine : str, optional
            Name of the conformer generation engine (e.g. "rdkit", "crest").
            Used for provenance in optimisation summary.
        pruning_method : str, optional
            Name of the pruning method and parameters (e.g. "topN_5").
            Used for provenance in optimisation summary.
        params : dict
            Engine-specific parameters, e.g.:
                - level="GFN2-xTB" for gxtb
                - functional="B3LYP", basis="def2-SVP" for orca
                - max_iter=250 for iteration limits

        Attributes
        ----------
        self.lookup_csv : str
            Path to the lookup file.
        self.output_dir : str
            Root directory for optimisation outputs.
        self.engine : str
            Lowercased engine name ("gxtb" or "orca").
        self.params : dict
            Engine-specific parameters.
        self.dataset : str
            Dataset name for provenance.
        self.generation_engine : str
            Generation engine name for provenance.
        self.pruning_method : str
            Pruning method string for provenance.
        """
        self.lookup_csv = lookup_csv
        self.output_dir = output_dir
        self.engine = engine.lower()
        self.params = params

        # Provenance fields for summary
        self.dataset = dataset
        self.generation_engine = generation_engine
        self.pruning_method = pruning_method

        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)


    def optimise(self):
        """
        Run optimisation for all conformers in the lookup file.
        Writes _optimisation_summary.csv at the top of the engine folder.
        """
        df = pd.read_csv(self.lookup_csv)
        results = []
        total_jobs = len(df)

        for idx, row in enumerate(df.itertuples(index=False), start=1):
            lookup_id = row.lookup_id
            xyz_file = self._find_xyz_file(lookup_id)
            if not xyz_file:
                print(f"[{idx}/{total_jobs}] Missing XYZ file for {lookup_id}, skipping.")
                continue

            if self.engine == "gxtb":
                energy, status, out_file, log_file = self._run_gxtb(xyz_file, lookup_id)
            elif self.engine == "orca":
                energy, status, out_file, log_file = self._run_orca(xyz_file, lookup_id)
            else:
                raise ValueError(f"Unknown optimisation engine: {self.engine}")

            results.append({
                "dataset": self.dataset,
                "generation_engine": self.generation_engine,
                "pruning_method": self.pruning_method,
                "optimisation_engine": f"{self.engine}_{self.params.get('level','GFN2')}",
                "lookup_id": lookup_id,
                "energy": energy,
                "status": status,
                "xyz_file": out_file,
                "log_file": log_file
            })

            # Print concise progress line with counter
            print(f"[{idx}/{total_jobs}] {lookup_id} optimisation {status}.")

        summary_path = os.path.join(self.output_dir, "_optimisation_summary.csv")
        pd.DataFrame(results).to_csv(summary_path, index=False)
        print(f"Optimisation summary written to {summary_path}")
        return summary_path

    # ------------------ Helpers ------------------

    def _find_xyz_file(self, lookup_id):
        """
        Locate the original XYZ file for a conformer.
        Assumes it was generated in 3_conformer_xyz/<engine>/...
        """
        # Example: ABCDEFG_conf0 → search in pipeline_data/3_conformer_xyz
        base_dir = "pipeline_data/3_conformer_xyz"
        for root, _, files in os.walk(base_dir):
            for f in files:
                if f.startswith(lookup_id) and f.endswith(".xyz"):
                    return os.path.join(root, f)
        return None

    def _run_engine(self, xyz_file, out_file):
        """
        Run optimisation with the chosen engine.
        Returns (energy, status).
        """
        if self.engine == "gxtb":
            return self._run_gxtb(xyz_file, out_file)
        elif self.engine == "orca":
            return self._run_orca(xyz_file, out_file)
        else:
            raise ValueError(f"Unknown optimisation engine: {self.engine}")

    def _run_gxtb(self, xyz_file, lookup_id):
        """
        Run GFN-xTB optimisation quietly, redirecting logs to log/ subfolder.
        Returns (energy, status, out_file, log_file).
        """
        level = self.params.get("level", "GFN2")
        max_iter = self.params.get("max_iter", 250)

        # Subfolders
        log_dir = os.path.join(self.output_dir, "log")
        xyz_dir = os.path.join(self.output_dir, "xyz")
        os.makedirs(log_dir, exist_ok=True)
        os.makedirs(xyz_dir, exist_ok=True)

        # File paths
        log_file = os.path.join(log_dir, f"{lookup_id}_xtb.log")
        out_file = os.path.join(xyz_dir, f"{lookup_id}_opt.xyz")

        cmd = ["xtb", xyz_file, "--opt", "--gfn", level.replace("GFN", ""), "--iterations", str(max_iter)]

        try:
            with open(log_file, "w") as log:
                subprocess.run(cmd, stdout=log, stderr=log, check=True)

            if os.path.exists("xtbopt.xyz"):
                os.rename("xtbopt.xyz", out_file)

            # Parse energy from log
            energy = None
            with open(log_file) as log:
                for line in log:
                    if "TOTAL ENERGY" in line.upper():
                        parts = line.split()
                        try:
                            energy = float(parts[-2])  # Hartree
                        except Exception:
                            pass

            status = "converged" if energy is not None else "unknown"
            #print(f"{lookup_id} optimised successfully.")
        except subprocess.CalledProcessError:
            energy, status = None, "failed"
            print(f"{lookup_id} optimisation failed.")

        return energy, status, out_file, log_file


    def _run_orca(self, xyz_file, out_file):
        """
        Run ORCA DFT optimisation.
        """
        functional = self.params.get("functional", "B3LYP")
        basis = self.params.get("basis", "def2-SVP")
        max_iter = self.params.get("max_iter", 250)

        # Write ORCA input file
        inp_file = xyz_file.replace(".xyz", ".inp")
        with open(inp_file, "w") as f:
            f.write(f"! {functional} {basis} Opt MaxIter {max_iter}\n")
            f.write("%output PrintLevel Mini\n")
            f.write("* xyz 0 1\n")
            with open(xyz_file) as xyz:
                lines = xyz.readlines()[2:]  # skip atom count + comment
                f.writelines(lines)
            f.write("*\n")

        cmd = ["orca", inp_file]
        try:
            subprocess.run(cmd, check=True)
            # Assume ORCA writes optimised geometry to <xyz_file>.xyz (placeholder)
            os.rename(xyz_file.replace(".xyz", "_opt.xyz"), out_file)
            energy = None
            status = "converged"
        except subprocess.CalledProcessError:
            energy = None
            status = "failed"
        return energy, status
