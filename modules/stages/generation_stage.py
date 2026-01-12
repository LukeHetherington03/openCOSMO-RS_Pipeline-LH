import os
import json
import time
import shutil
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from modules.utils.atomic_write import AtomicWriter


GENERATION_METHODS = {
    "rdkit": {"runner": "_generate_rdkit"},
    "babel": {"runner": "_generate_babel"},
    "crest": {"runner": "_generate_crest"},
}


class GenerationStage:
    """
    Generation stage:
      - Accepts one or more CSVs
      - Copies them into inputs/
      - Loads + concatenates them
      - Validates molecules
      - Generates conformers using selected engine
      - Writes energies.json, summary.csv, xyz/
    """

    # ------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------
    def run(self, job):
        params = job.parameters

        input_csvs = params.get("input_csv")
        num_confs = params.get("num_confs", 5)
        seed = params.get("seed", 42)
        engine = params.get("engine", "rdkit").lower()

        # Normalise CSV list
        if isinstance(input_csvs, str):
            input_csvs = [input_csvs]
        if not input_csvs:
            raise ValueError("GenerationStage requires at least one input CSV.")

        job.log_header("Starting generation stage")
        job.log(f"Engine: {engine}")
        job.log(f"Input CSVs: {input_csvs}")
        job.log(f"Requested: {num_confs} conformers per molecule")

        # Prepare directories
        inputs_dir, outputs_dir, xyz_out_dir, raw_out_dir = self._prepare_dirs(job)

        # Copy CSVs
        copied_csvs, missing_csvs = self._copy_csvs(job, input_csvs, inputs_dir)

        # Load CSVs
        df_all = self._load_csvs(job, copied_csvs)

        # Validate molecules
        valid_rows, invalid_rows = self._validate_molecules(job, df_all)

        # Dispatch engine
        results = self._dispatch_engine(
            job, engine, valid_rows, num_confs, seed, xyz_out_dir
        )

        # Write outputs
        self._write_outputs(
            job, results, outputs_dir,
            input_csvs, copied_csvs, missing_csvs,
            valid_rows, invalid_rows
        )

        job.log_header("Generation complete")
        job.log(f"{len(valid_rows)} valid molecules processed")
        job.log(f"{len(results)} conformers generated")
        job.log(f"Outputs written to: {outputs_dir}")

        job.mark_complete()

    # ------------------------------------------------------------
    # Directory setup
    # ------------------------------------------------------------
    def _prepare_dirs(self, job):
        inputs_dir = job.inputs_dir
        outputs_dir = job.outputs_dir

        xyz_out_dir = os.path.join(outputs_dir, "xyz")
        raw_out_dir = os.path.join(outputs_dir, "raw")

        os.makedirs(xyz_out_dir, exist_ok=True)
        os.makedirs(raw_out_dir, exist_ok=True)

        return inputs_dir, outputs_dir, xyz_out_dir, raw_out_dir

    # ------------------------------------------------------------
    # Copy CSVs
    # ------------------------------------------------------------
    def _copy_csvs(self, job, input_csvs, inputs_dir):
        copied_csvs = []
        missing_csvs = []

        for csv_path in input_csvs:
            if os.path.exists(csv_path):
                dst = os.path.join(inputs_dir, os.path.basename(csv_path))
                shutil.copy(csv_path, dst)
                copied_csvs.append(dst)
            else:
                missing_csvs.append(csv_path)

        if len(copied_csvs) == 0:
            raise ValueError(f"None of the provided CSVs exist. Missing: {missing_csvs}")

        if missing_csvs:
            job.log("[WARNING] Missing CSVs:")
            for m in missing_csvs:
                job.log(f"  - {m}", indent=2)

        return copied_csvs, missing_csvs

    # ------------------------------------------------------------
    # Load CSVs
    # ------------------------------------------------------------
    def _load_csvs(self, job, copied_csvs):
        dfs = []
        for csv in copied_csvs:
            try:
                df = pd.read_csv(csv)
                dfs.append(df)
            except Exception as e:
                job.log(f"[WARNING] Failed to load {csv}: {e}")

        if len(dfs) == 0:
            raise ValueError("No CSVs could be loaded successfully.")

        df_all = pd.concat(dfs, ignore_index=True)
        job.log(f"Loaded {len(df_all)} rows from CSVs")

        return df_all

    # ------------------------------------------------------------
    # Validate molecules
    # ------------------------------------------------------------
    def _validate_molecules(self, job, df_all):
        valid_rows = []
        invalid_rows = []

        for idx, row in df_all.iterrows():
            smiles = row.get("smiles") or row.get("SMILES") or row.get("mol")
            if not smiles:
                invalid_rows.append(idx)
                continue

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                invalid_rows.append(idx)
                continue

            valid_rows.append((idx, smiles))

        job.log(f"Valid molecules: {len(valid_rows)}")
        if invalid_rows:
            job.log(f"[WARNING] {len(invalid_rows)} invalid molecule rows skipped")

        return valid_rows, invalid_rows

    # ------------------------------------------------------------
    # Engine dispatch
    # ------------------------------------------------------------
    def _dispatch_engine(self, job, engine, valid_rows, num_confs, seed, xyz_out_dir):
        if engine not in GENERATION_METHODS:
            raise ValueError(f"Unknown generation engine: {engine}")

        runner = getattr(self, GENERATION_METHODS[engine]["runner"])
        return runner(job, valid_rows, num_confs, seed, xyz_out_dir)

    # ------------------------------------------------------------
    # RDKit backend (default)
    # ------------------------------------------------------------
    def _generate_rdkit(self, job, valid_rows, num_confs, seed, xyz_out_dir):
        results = []
        total_mols = len(valid_rows)
        global_start = time.perf_counter()
        total_confs = 0

        for mol_idx, (row_idx, smiles) in enumerate(valid_rows, start=1):
            job.log_section(f"Molecule {mol_idx}/{total_mols}: {smiles}")

            mol_start = time.perf_counter()
            mol_results = self._generate_rdkit_for_molecule(
                job, smiles, row_idx, num_confs, seed, xyz_out_dir
            )

            elapsed = time.perf_counter() - mol_start
            total_confs += len(mol_results)

            job.log(f"Generated {len(mol_results)} conformers in {elapsed:.2f} seconds", indent=1)
            job.log(f"Total so far: {total_confs}", indent=1)

            results.extend(mol_results)

        stage_elapsed = time.perf_counter() - global_start
        job.log(f"RDKit generation time: {stage_elapsed:.2f} seconds")

        return results

    def _generate_rdkit_for_molecule(self, job, smiles, row_idx, num_confs, seed, xyz_out_dir):
        results = []

        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))

        params = AllChem.ETKDGv3()
        params.randomSeed = seed

        try:
            conf_ids = AllChem.EmbedMultipleConfs(mol, num_confs, params)
        except Exception as e:
            job.log(f"[WARNING] Failed to generate conformers for {smiles}: {e}", indent=1)
            return results

        for conf_id in conf_ids:
            try:
                AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
                energy = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id).CalcEnergy()
            except Exception as e:
                job.log(f"[WARNING] UFF optimisation failed for {smiles}: {e}", indent=2)
                continue

            inchi_key = Chem.inchi.MolToInchiKey(mol)
            lookup_id = f"{inchi_key}_conf{conf_id}"

            xyz_path = os.path.join(xyz_out_dir, f"{lookup_id}.xyz")
            self._write_xyz(mol, conf_id, xyz_path)

            results.append({
                "lookup_id": lookup_id,
                "energy": energy,
                "xyz_path": xyz_path,
                "log_path": None,
                "metadata": {
                    "smiles": smiles,
                    "source_row": int(row_idx),
                }
            })

        return results

    # ------------------------------------------------------------
    # OpenBabel backend (placeholder)
    # ------------------------------------------------------------
    def _generate_babel(self, job, valid_rows, num_confs, seed, xyz_out_dir):
        job.log("[ERROR] OpenBabel generation not implemented yet.")
        raise NotImplementedError("OpenBabel generation not implemented yet.")

    # ------------------------------------------------------------
    # CREST backend (placeholder)
    # ------------------------------------------------------------
    def _generate_crest(self, job, valid_rows, num_confs, seed, xyz_out_dir):
        job.log("[ERROR] CREST generation not implemented yet.")
        raise NotImplementedError("CREST generation not implemented yet.")

    # ------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------
    def _write_outputs(self, job, results, outputs_dir,
                       input_csvs, copied_csvs, missing_csvs,
                       valid_rows, invalid_rows):

        if len(results) == 0:
            raise ValueError("Generation produced zero conformers.")

        # summary.csv
        summary_path = os.path.join(outputs_dir, "summary.csv")
        pd.DataFrame(results).to_csv(summary_path, index=False)

        # energies.json
        energies_path = os.path.join(outputs_dir, "energies.json")
        with AtomicWriter(energies_path) as f:
            json.dump(results, f, indent=2)

        # job_state.json
        job_state_path = os.path.join(job.job_dir, "job_state.json")
        with AtomicWriter(job_state_path) as f:
            json.dump(
                {
                    "stage": "generation",
                    "num_input_csvs": len(input_csvs),
                    "num_copied_csvs": len(copied_csvs),
                    "missing_csvs": missing_csvs,
                    "num_valid_molecules": len(valid_rows),
                    "num_invalid_molecules": len(invalid_rows),
                    "num_conformers": len(results),
                },
                f,
                indent=2,
            )

    # ------------------------------------------------------------
    # XYZ writer
    # ------------------------------------------------------------
    def _write_xyz(self, mol, conf_id, path):
        atoms = mol.GetAtoms()
        conf = mol.GetConformer(conf_id)

        with open(path, "w") as f:
            f.write(f"{len(atoms)}\n")
            f.write("Generated by GenerationStage\n")
            for atom in atoms:
                pos = conf.GetAtomPosition(atom.GetIdx())
                f.write(
                    f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n"
                )
