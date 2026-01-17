import os
import json
import time
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter


GENERATION_BACKENDS = {
    "rdkit": "_backend_rdkit",
    "babel": "_backend_babel",
    "crest": "_backend_crest",
}


class GenerationStage(BaseStage):
    """
    Clean, unified conformer generation stage.

    Responsibilities:
      - Load + validate cleaned dataset (summary_file)
      - Canonicalise SMILES
      - Generate conformers using selected backend
      - Write:
          - xyz/ files
          - summary.csv
          - energies.json
    """

    # ------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------
    def execute(self):
        self.log_header("Starting Generation Stage")

        params = self.parameters

        # ------------------------------------------------------------
        # NEW: Use summary_file from CleaningStage
        # ------------------------------------------------------------
        summary_file = params.get("summary_file")
        if summary_file is None:
            self.fail("GenerationStage requires summary_file produced by CleaningStage.")

        if not os.path.isfile(summary_file):
            self.fail(f"summary_file does not exist: {summary_file}")

        input_csvs = [summary_file]  # keep interface consistent

        num_confs = params.get("num_confs", 5)
        seed = params.get("seed", 42)
        backend = params.get("engine", "rdkit").lower()

        self.log(f"Backend: {backend}")
        self.log(f"Input summary file: {summary_file}")
        self.log(f"Conformers per molecule: {num_confs}")

        # Prepare directories
        dirs = self._prepare_directories()

        # Load cleaned dataset
        df_all = self._load_and_merge_csvs(input_csvs)

        # Validate molecules
        valid_rows = self._validate_and_canonicalise(df_all)

        # Initialise job items (one per valid molecule)
        self.set_items([idx for idx, _ in valid_rows])

        # Dispatch backend
        results = self._dispatch_backend(
            backend=backend,
            valid_rows=valid_rows,
            num_confs=num_confs,
            seed=seed,
            xyz_out_dir=dirs["xyz"]
        )

        # Write outputs
        self._write_outputs(
            results=results,
            input_csvs=input_csvs,
            dirs=dirs,
            valid_rows=valid_rows,
            total_rows=len(df_all)
        )

        self.log_header("Generation Stage Complete")
        self.job.mark_complete()

    # ------------------------------------------------------------
    # Helpers: Input handling
    # ------------------------------------------------------------
    def _normalise_csv_list(self, csvs):
        if isinstance(csvs, str):
            return [csvs]
        if not csvs:
            self.fail("GenerationStage requires at least one input CSV.")
        return csvs

    def _prepare_directories(self):
        dirs = {
            "inputs": self.inputs_dir,
            "outputs": self.outputs_dir,
            "xyz": os.path.join(self.outputs_dir, "xyz"),
            "raw": os.path.join(self.outputs_dir, "raw"),
        }
        for d in dirs.values():
            os.makedirs(d, exist_ok=True)
        return dirs

    def _load_and_merge_csvs(self, csv_paths):
        dfs = []

        for path in csv_paths:
            if not os.path.exists(path):
                self.log(f"[WARNING] Missing CSV: {path}")
                continue
            try:
                df = pd.read_csv(path)
                dfs.append(df)
            except Exception as e:
                self.log(f"[WARNING] Failed to load {path}: {e}")

        if not dfs:
            self.fail("No valid CSVs loaded.")

        merged = pd.concat(dfs, ignore_index=True)
        self.log(f"Loaded {len(merged)} rows from CSVs")
        return merged

    # ------------------------------------------------------------
    # Helpers: Molecule validation
    # ------------------------------------------------------------
    def _validate_and_canonicalise(self, df):
        valid = []

        for idx, row in df.iterrows():
            smiles = row.get("smiles") or row.get("SMILES") or row.get("mol")
            if not smiles:
                continue

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue

            smiles = Chem.MolToSmiles(mol, canonical=True)
            valid.append((idx, smiles))

        self.log(f"Valid molecules: {len(valid)}")
        return valid

    # ------------------------------------------------------------
    # Backend dispatch
    # ------------------------------------------------------------
    def _dispatch_backend(self, backend, valid_rows, num_confs, seed, xyz_out_dir):

        if backend not in GENERATION_BACKENDS:
            self.fail(f"Unknown backend: {backend}")

        runner_name = GENERATION_BACKENDS[backend]
        runner = getattr(self, runner_name)

        self.log_section(f"Running backend: {backend}")
        return runner(valid_rows, num_confs, seed, xyz_out_dir)

    # ------------------------------------------------------------
    # Backend: RDKit
    # ------------------------------------------------------------
    def _backend_rdkit(self, valid_rows, num_confs, seed, xyz_out_dir):
        results = []

        for idx, smiles in valid_rows:
            self.log_section(f"Generating for SMILES: {smiles}")

            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            params = AllChem.ETKDGv3()
            params.randomSeed = seed

            try:
                conf_ids = AllChem.EmbedMultipleConfs(mol, num_confs, params)
            except Exception as e:
                self.log(f"[WARNING] RDKit embedding failed: {e}")
                self.update_progress(idx, success=False)
                continue

            for conf_id in conf_ids:
                try:
                    AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
                    energy = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id).CalcEnergy()
                except Exception as e:
                    self.log(f"[WARNING] UFF optimisation failed: {e}")
                    continue

                inchi_key = Chem.inchi.MolToInchiKey(mol)
                lookup_id = f"{inchi_key}_conf{conf_id:03d}"

                xyz_path = os.path.join(xyz_out_dir, f"{lookup_id}.xyz")
                self._write_xyz(mol, conf_id, xyz_path)

                results.append({
                    "lookup_id": lookup_id,
                    "energy": energy,
                    "xyz_path": xyz_path,
                    "metadata": {
                        "smiles": smiles,
                        "source_row": int(idx),
                    }
                })

            self.update_progress(idx)

        return results

    # ------------------------------------------------------------
    # Backend placeholders
    # ------------------------------------------------------------
    def _backend_babel(self, *args, **kwargs):
        self.fail("OpenBabel backend not implemented yet.")

    def _backend_crest(self, *args, **kwargs):
        self.fail("CREST backend not implemented yet.")

    # ------------------------------------------------------------
    # Helpers: Output writing
    # ------------------------------------------------------------
    def _write_outputs(self, results, input_csvs, dirs, valid_rows, total_rows):

        if not results:
            self.fail("Generation produced zero conformers.")

        # summary.csv
        summary_path = os.path.join(dirs["outputs"], "summary.csv")
        with AtomicWriter(summary_path) as f:
            pd.DataFrame(results).to_csv(f, index=False)

        # energies.json
        energies_path = os.path.join(dirs["outputs"], "energies.json")
        with AtomicWriter(energies_path) as f:
            json.dump(results, f, indent=2)

        # job_state.json
        job_state_path = os.path.join(self.job.job_dir, "job_state.json")
        with AtomicWriter(job_state_path) as f:
            json.dump(
                {
                    "stage": "generation",
                    "num_input_csvs": len(input_csvs),
                    "num_valid_molecules": len(valid_rows),
                    "num_total_rows": total_rows,
                    "num_conformers": len(results),
                },
                f,
                indent=2,
            )

    # ------------------------------------------------------------
    # Helpers: XYZ writer
    # ------------------------------------------------------------
    def _write_xyz(self, mol, conf_id, path):
        atoms = mol.GetAtoms()
        conf = mol.GetConformer(conf_id)

        with open(path, "w") as f:
            f.write(f"{len(atoms)}\n")
            f.write("Generated by GenerationStage\n")
            for atom in atoms:
                pos = conf.GetAtomPosition(atom.GetIdx())
                f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")
