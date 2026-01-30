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
          - molecule_metadata/<InChIKey>.json (job-local)
          - GLOBAL_METADATA_DIR/<InChIKey>.json (global, no overwrite)
    """

    # ------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------
    def execute(self):
        self.log_header("Starting Generation Stage")

        cfg = self.config
        if not cfg:
            self.fail("GenerationStage: missing config in job.config")

        self.GLOBAL_METADATA_DIR = cfg["constant_files"]["metadata_dir"]
        os.makedirs(self.GLOBAL_METADATA_DIR, exist_ok=True)

        params = self.parameters

        summary_file = params.get("summary_file")
        if summary_file is None:
            self.fail("GenerationStage requires summary_file produced by CleaningStage.")

        if not os.path.isfile(summary_file):
            self.fail(f"summary_file does not exist: {summary_file}")

        input_csvs = [summary_file]

        num_confs = params.get("num_confs", 5)
        seed = params.get("seed", 42)
        backend = params.get("engine", "rdkit").lower()

        self.log(f"Backend: {backend}")
        self.log(f"Input summary file: {summary_file}")
        self.log(f"Conformers per molecule: {num_confs}")

        dirs = self._prepare_directories()

        df_all = self._load_and_merge_csvs(input_csvs)

        valid_rows = self._validate_and_canonicalise(df_all)

        self.set_items([idx for idx, *_ in valid_rows])

        results = self._dispatch_backend(
            backend=backend,
            valid_rows=valid_rows,
            num_confs=num_confs,
            seed=seed,
            xyz_out_dir=dirs["xyz"]
        )

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
    def _prepare_directories(self):
        dirs = {
            "inputs": self.inputs_dir,
            "outputs": self.outputs_dir,
            "xyz": os.path.join(self.outputs_dir, "xyz"),
            "raw": os.path.join(self.outputs_dir, "raw"),
            "metadata": os.path.join(self.outputs_dir, "molecule_metadata"),
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

            # ----------------------------------------------------
            # NEW: Normalise melting_temp and melting_temp_source
            # ----------------------------------------------------
            melting_temp = row.get("melting_temp")
            if pd.notna(melting_temp):
                try:
                    melting_temp = float(str(melting_temp).strip())
                except Exception:
                    melting_temp = None
            else:
                melting_temp = None

            melting_temp_source = row.get("melting_temp_source")
            if isinstance(melting_temp_source, str):
                melting_temp_source = melting_temp_source.strip()
                if melting_temp_source == "":
                    melting_temp_source = None

            # Store back into row so backend sees correct values
            row["melting_temp"] = melting_temp
            row["melting_temp_source"] = melting_temp_source

            valid.append((idx, smiles, row))

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

        for idx, smiles, row in valid_rows:
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

            inchi = Chem.inchi.MolToInchi(mol)
            inchi_key = Chem.inchi.InchiToInchiKey(inchi)

            melting_temp = row.get("melting_temp")
            melting_temp_source = row.get("melting_temp_source")

            for conf_id in conf_ids:
                try:
                    AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
                    energy = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id).CalcEnergy()
                except Exception as e:
                    self.log(f"[WARNING] UFF optimisation failed: {e}")
                    continue

                lookup_id = f"{inchi_key}_conf{conf_id:03d}"

                xyz_path = os.path.join(xyz_out_dir, f"{lookup_id}.xyz")
                self._write_xyz(mol, conf_id, xyz_path)

                results.append({
                    "lookup_id": lookup_id,
                    "inchi": inchi,
                    "inchi_key": inchi_key,
                    "energy": energy,
                    "xyz_path": xyz_path,
                    "metadata": {
                        "smiles": smiles,
                        "source_row": int(idx),
                        "melting_temp": melting_temp,
                        "melting_temp_source": melting_temp_source,
                    }
                })

            self.update_progress(idx)

        return results

    # ------------------------------------------------------------
    # Helpers: Output writing
    # ------------------------------------------------------------
    def _write_outputs(self, results, input_csvs, dirs, valid_rows, total_rows):

        if not results:
            self.fail("Generation produced zero conformers.")

        summary_path = os.path.join(dirs["outputs"], "summary.csv")
        with AtomicWriter(summary_path) as f:
            pd.DataFrame(results).to_csv(f, index=False)

        energies_path = os.path.join(dirs["outputs"], "energies.json")
        with AtomicWriter(energies_path) as f:
            json.dump(results, f, indent=2)

        self._write_metadata_templates(results, dirs["metadata"])

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
    # Write molecule metadata templates
    # ------------------------------------------------------------
    def _write_metadata_templates(self, results, job_metadata_dir):

        molecules = {}
        for entry in results:
            key = entry["inchi_key"]
            molecules.setdefault(key, entry)

        for inchi_key, entry in molecules.items():

            md = entry["metadata"]

            template = {
                "lookup_id": inchi_key,
                "inchi_key": inchi_key,
                "inchi": entry["inchi"],
                "smiles": md["smiles"],

                "melting_temp": md.get("melting_temp", None),
                "melting_temp_source": md.get("melting_temp_source", None),

                "Hfus": "N/A",
                "Gfus_model": "MyrdalYalkowsky",

                "notes": "Fill in melting point, Hfus, and Gfus model if known."
            }

            job_path = os.path.join(job_metadata_dir, f"{inchi_key}.json")
            with AtomicWriter(job_path) as f:
                json.dump(template, f, indent=2)

            global_path = os.path.join(self.GLOBAL_METADATA_DIR, f"{inchi_key}.json")

            if not os.path.exists(global_path):
                with AtomicWriter(global_path) as f:
                    json.dump(template, f, indent=2)
            else:
                self.log(f"[INFO] Global metadata exists, not overwriting: {global_path}")

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
