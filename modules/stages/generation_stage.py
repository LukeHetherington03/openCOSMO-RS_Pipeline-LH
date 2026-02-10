import os
import json
import pandas as pd
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import AllChem

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter
from modules.utils.conformers import ConformerRecord, ConformerSet


GENERATION_BACKENDS = {
    "rdkit": "_backend_rdkit",
    # future: "babel": "_backend_babel",
    # future: "crest": "_backend_crest",
}


class GenerationStage(BaseStage):
    """
    Conformer generation stage.

    Input:
        stage_input = cleaned.csv from CleaningStage

    Output:
        energies.json  (canonical stage output)
        summary.csv
        xyz/ directory
    """

    # ------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------
    def execute(self):
        # Strict mode
        self.strict_mode = self.strict("generation")

        # Declare canonical output file
        self.set_stage_output("energies.json")

        # Retrieve stage input
        stage_input = self.get_stage_input()
        self.require_file(stage_input, "cleaned dataset")

        params = self.parameters
        num_confs = params.get("num_confs", 5)
        seed = params.get("seed", 42)
        backend = params.get("engine", "rdkit").lower()

        self.log(f"Backend: {backend}")
        self.log(f"Stage input: {stage_input}")
        self.log(f"Conformers per molecule: {num_confs}")
        self.log(f"Strict mode: {self.strict_mode}")

        # Prepare directories
        dirs = self._prepare_directories()

        # Load cleaned dataset
        df_all = self._load_csv(stage_input)

        # Validate rows
        valid_rows = self._validate_rows(df_all)
        self.set_items([row_idx for row_idx, *_ in valid_rows])

        self.warnings = {
            "embedding_failures": 0,
            "optimisation_failures": 0,
            "zero_conformers": 0,
        }

        # Run backend
        conformer_set = self._dispatch_backend(
            backend=backend,
            valid_rows=valid_rows,
            num_confs=num_confs,
            seed=seed,
            xyz_out_dir=dirs["xyz"],
        )

        if len(conformer_set) == 0:
            self.fail("Generation produced zero conformers.")

        # Write outputs
        self._write_outputs(
            conformer_set=conformer_set,
            dirs=dirs,
            num_valid=len(valid_rows),
            total_rows=len(df_all),
            backend=backend,
            num_confs=num_confs,
            seed=seed,
        )

        self._log_warning_summary()

    # ------------------------------------------------------------
    # Helpers: Input handling
    # ------------------------------------------------------------
    def _prepare_directories(self):
        dirs = {
            "inputs": self.inputs_dir,
            "outputs": self.outputs_dir,
            "xyz": os.path.join(self.outputs_dir, "xyz"),
        }
        for d in dirs.values():
            os.makedirs(d, exist_ok=True)
        return dirs

    def _load_csv(self, path: str) -> pd.DataFrame:
        try:
            df = pd.read_csv(path)
        except Exception as e:
            self.fail(f"Failed to load stage_input {path}: {e}")
        self.log(f"Loaded {len(df)} rows from {path}")
        return df

    # ------------------------------------------------------------
    # Helpers: Molecule validation
    # ------------------------------------------------------------
    def _validate_rows(self, df: pd.DataFrame):
        """
        Validate presence of inchi_key and smiles, canonicalise smiles.
        Returns a list of tuples: (row_idx, inchi_key, smiles, row)
        """
        valid = []

        for idx, row in df.iterrows():
            inchi_key = row.get("inchi_key")
            smiles = row.get("smiles")

            if not inchi_key or not isinstance(inchi_key, str):
                msg = f"Row {idx}: missing or invalid inchi_key"
                self.log(f"[ERROR] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                continue

            if not smiles:
                msg = f"Row {idx}: missing SMILES"
                self.log(f"[ERROR] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                continue

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                msg = f"Row {idx}: invalid SMILES '{smiles}'"
                self.log(f"[WARNING] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                continue

            smiles_canon = Chem.MolToSmiles(mol, canonical=True)
            df.at[idx, "smiles"] = smiles_canon

            valid.append((idx, inchi_key, smiles_canon, row))

        self.log(f"Valid molecules for generation: {len(valid)}")
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
        conformer_set = ConformerSet()
        timestamp = datetime.utcnow().isoformat() + "Z"

        for idx, inchi_key, smiles, row in valid_rows:
            self.log_section(f"Generating conformers for {inchi_key} ({smiles})")

            base_mol = Chem.MolFromSmiles(smiles)
            if base_mol is None:
                msg = f"Row {idx}: RDKit failed to parse SMILES '{smiles}'"
                self.log(f"[WARNING] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.warnings["embedding_failures"] += 1
                self.update_progress(idx, success=False)
                continue

            mol = Chem.AddHs(base_mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = seed

            try:
                conf_ids = AllChem.EmbedMultipleConfs(mol, num_confs, params)
            except Exception as e:
                msg = f"RDKit embedding failed for {inchi_key}: {e}"
                self.log(f"[WARNING] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.warnings["embedding_failures"] += 1
                self.update_progress(idx, success=False)
                continue

            if len(conf_ids) == 0:
                msg = f"No conformers generated for {inchi_key}"
                self.log(f"[WARNING] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.warnings["zero_conformers"] += 1
                self.update_progress(idx, success=False)
                continue

            for local_conf_idx, conf_id in enumerate(conf_ids):
                conf_num = int(local_conf_idx)

                try:
                    AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
                    ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                    energy = float(ff.CalcEnergy())
                except Exception as e:
                    msg = f"UFF optimisation failed for {inchi_key} conf {conf_num}: {e}"
                    self.log(f"[WARNING] {msg}")
                    if self.strict_mode:
                        self.fail(msg)
                    self.warnings["optimisation_failures"] += 1
                    continue

                xyz_path = os.path.join(
                    xyz_out_dir,
                    f"{inchi_key}_conf{conf_num:03d}.xyz",
                )

                record = ConformerRecord(
                    inchi_key=inchi_key,
                    conf_num=conf_num,
                    xyz_path=xyz_path,
                    energy=energy,
                    smiles=smiles,
                    provenance={
                        "backend": "rdkit",
                        "seed": seed,
                        "generation_timestamp": timestamp,
                        "source_row": int(idx),
                    },
                )

                # 0th optimisation entry (generation-level geometry + energy)
                record.optimisation_history.append({
                    "stage": "generation",
                    "engine": "rdkit",
                    "energy": energy,
                    "xyz_path": xyz_path,
                    "timestamp": timestamp,
                })

                self._write_xyz(mol, conf_id, xyz_path)
                conformer_set.add(record)

            self.update_progress(idx)

        return conformer_set

    # ------------------------------------------------------------
    # Helpers: Output writing
    # ------------------------------------------------------------
    def _write_outputs(
        self,
        conformer_set: ConformerSet,
        dirs,
        num_valid,
        total_rows,
        backend,
        num_confs,
        seed,
    ):
        # summary.csv
        summary_path = os.path.join(dirs["outputs"], "summary.csv")
        with AtomicWriter(summary_path) as f:
            pd.DataFrame(conformer_set.to_list()).to_csv(f, index=False)

        # energies.json (canonical output)
        energies_path = self.get_stage_output()
        conformer_set.save(energies_path)

        # Optional: lightweight metadata file (not job_state.json)
        meta_path = os.path.join(self.outputs_dir, "generation_metadata.json")
        with AtomicWriter(meta_path) as f:
            json.dump(
                {
                    "stage": "generation",
                    "backend": backend,
                    "num_valid_molecules": num_valid,
                    "num_total_rows": total_rows,
                    "num_conformers": len(conformer_set),
                    "num_confs_requested": num_confs,
                    "seed": seed,
                    "timestamp": datetime.utcnow().isoformat() + "Z",
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

    # ------------------------------------------------------------
    # Warning summary
    # ------------------------------------------------------------
    def _log_warning_summary(self):
        total = sum(self.warnings.values())
        if total == 0:
            self.log("Generation completed with no warnings.")
            return

        self.log("Generation completed with warnings:")
        if self.warnings["embedding_failures"] > 0:
            self.log(
                f" - Embedding failures: {self.warnings['embedding_failures']} molecule(s)."
            )
        if self.warnings["optimisation_failures"] > 0:
            self.log(
                f" - Optimisation failures: {self.warnings['optimisation_failures']} conformer(s)."
            )
        if self.warnings["zero_conformers"] > 0:
            self.log(
                f" - Zero conformers generated for "
                f"{self.warnings['zero_conformers']} molecule(s)."
            )
