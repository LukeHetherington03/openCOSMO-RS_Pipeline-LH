"""
conformer_generation.py

Adaptive conformer generator with interchangeable backends (RDKit, CREST, Balloon).
Each backend generates conformers and returns them. A shared save_conformers()
method handles writing conformers to disk with consistent naming.

Usage:
    from conformer_generation import ConformerGenerator

    gen = ConformerGenerator(method="rdkit", output_dir="xyz_output")
    mol, conf_ids = gen.generate(smiles="CCO", charge=0, num_confs=50)
"""

import os
import csv
import subprocess as spr
from rdkit import Chem
from rdkit.Chem import AllChem
from .molecule_utils import MoleculeUtils


class ConformerGenerator:
    def __init__(self, dataset: str, method: str = "rdkit", output_dir: str = "xyz_output"):
        """
        Parameters
        ----------
        dataset : str
            Name of the dataset (e.g. 'acr_t5').
        method : str
            Backend engine ('rdkit', 'crest', 'balloon').
        output_dir : str
            Directory to write conformer files and master CSV.
        """
        self.dataset = dataset
        self.method = method.lower()
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def generate(self, smiles, charge: int = 0, num_confs: int = None, seed: int = 42):
        """
        Generate conformers using the selected backend.
        Returns: dict with master_csv path, xyz folder, and conformer records.
        """
        if self.method == "rdkit":
            records = self._generate_rdkit(smiles, charge, num_confs, seed)
        elif self.method == "crest":
            records = self._generate_crest(smiles, charge)
        elif self.method == "balloon":
            records = self._generate_balloon(smiles, charge, num_confs)
        else:
            raise ValueError(f"Unknown conformer generation method: {self.method}")

        # Save conformers and master CSV
        self._save_conformers(records)

        # Build master CSV path with new naming convention
        master_csv_name = f"_master_{self.dataset}_{self.method}_{num_confs}.csv"
        master_csv_path = os.path.join(self.output_dir, master_csv_name)

        # Return structured outputs
        return {
            "dataset": self.dataset,
            "engine": self.method,
            "conformer_xyz_folder": self.output_dir,
            "master_csv": master_csv_path,
        }

    # ------------------ Shared Save ------------------

    def _save_conformers(self, conformer_records):
        """
        Save conformers to disk with consistent naming and initialise master CSV.
        """
        num_confs = len(conformer_records)
        # New naming convention: _master_<dataset>_<engine>_<n_confs>.csv
        csv_name = f"_master_{self.dataset}_{self.method}_{num_confs}.csv"
        csv_path = os.path.join(self.output_dir, csv_name)

        write_header = not os.path.exists(csv_path)

        with open(csv_path, "a", newline="") as f_csv:
            writer = csv.writer(f_csv)
            if write_header:
                writer.writerow([
                    "lookup_id", "inchi_key", "conf_id",
                    "energy_gen", "status_gen",
                    "method", "charge", "rot_bonds", "filename",
                    "energy_xtb", "status_xtb",
                    "energy_gxtb", "status_gxtb",
                    "energy_dft", "status_dft"
                ])

            for record in conformer_records:
                inchi_key = record["inchi_key"]
                conf_id   = record["conf_id"]
                method    = record["method"]
                energy    = record.get("energy", "NA")
                charge    = record.get("charge", "NA")
                rot_bonds = record.get("rot_bonds", "NA")

                rel_name = f"{inchi_key}_conf{conf_id}.xyz"
                xyz_path = os.path.join(self.output_dir, rel_name)

                if record.get("mol") is not None:
                    MoleculeUtils.save_xyz(
                        record["mol"], conf_id, xyz_path,
                        inchi_key, method
                    )
                elif record.get("file") is not None:
                    os.rename(record["file"], xyz_path)

                lookup_id = f"{inchi_key}_conf{conf_id}"

                writer.writerow([
                    lookup_id, inchi_key, conf_id,
                    energy, "generated",
                    method, charge, rot_bonds, rel_name,
                    "", "", "", "", "", ""
                ])

        # ------------------ Backends ------------------
    def _generate_rdkit(self, smiles, charge, num_confs, seed):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        nrot = MoleculeUtils.rotatable_bonds(mol)
        num_confs = num_confs or self._choose_num_confs(nrot)

        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        params.pruneRmsThresh = -1
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)

        inchi_str = Chem.inchi.MolToInchi(mol)
        inchi_key = Chem.inchi.InchiToInchiKey(inchi_str)

        records = []
        for conf_id in conf_ids:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
            energy = ff.CalcEnergy()
            records.append({
                "inchi_key": inchi_key,
                "conf_id": conf_id,
                "mol": mol,
                "file": None,
                "energy": energy,
                "method": self.method,
                "charge": charge,
                "rot_bonds": nrot
            })
        return records


    def _generate_crest(self, smiles, charge):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        xyz_file = os.path.join(self.output_dir, f"{smiles}.xyz")
        MoleculeUtils.save_xyz(mol, 0, xyz_file)

        spr.run(["crest", xyz_file], cwd=self.output_dir, check=False)
        crest_out = os.path.join(self.output_dir, "crest_conformers.xyz")

        return mol, [crest_out]

    def _generate_balloon(self, smiles, charge, num_confs):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        balloon_exe = "balloon"
        sdf_out = os.path.join(self.output_dir, f"{smiles}.sdf")
        mol_file = os.path.join(self.output_dir, f"{smiles}.mol")
        Chem.MolToMolFile(mol, mol_file)

        spr.run([balloon_exe, "--nconfs", str(num_confs), "--output-format", "sdf",
                 mol_file, sdf_out], check=False)

        return mol, [sdf_out]




    # ------------------ Utility ------------------

    def _choose_num_confs(self, nrot, override=None):
        if override is not None:
            return override
        if nrot <= 7:
            return 50
        elif 8 <= nrot <= 12:
            return 200
        else:
            return 300
