"""
data_cleaning.py

Module to clean and standardise chemical datasets for use in the conformer/COSMO pipeline.
Contains one class: DataCleaner

Features:
- Fuzzy searches for 'name' or 'mol_name' column and maps it to 'mol_name_iupac'
- Ensures required columns: key_inchi, smiles, charge
- Generates InChIKey from SMILES using RDKit
- Defaults charge to 0 if missing
- Fuzzy maps other headers to a standardised schema
- Outputs a cleaned CSV with consistent column order
- Writes a twin .inp file for conformer generation

Usage:
    from data_cleaning import DataCleaner

    cleaner = DataCleaner("raw_dataset.csv", "clean_dataset.csv")
    cleaner.clean()
"""

import pandas as pd
import difflib
from rdkit import Chem
from rdkit.Chem import inchi
import os

STANDARD_NAMES = {
    "experimental solubility /mol frac": "experimental_solubility_mol_frac",
    "melting point": "melting_point",
    "melting point source": "melting_point_source",
    "mol_id": "mol_id",
    "pubchem_sid": "pubchem_sid",
    "mol_name": "mol_name_iupac",
    "name": "mol_name_iupac",
    "smiles": "smiles",
    "formula_calcd": "formula_calcd",
    "anionic": "anionic"
}

class DataCleaner:
    def __init__(self, input_csv, output_csv):
        self.input_csv = input_csv
        self.output_csv = output_csv
        self.df = None

    def _fuzzy_match(self, col):
        matches = difflib.get_close_matches(col.lower(), STANDARD_NAMES.keys(), n=1, cutoff=0.6)
        if matches:
            return STANDARD_NAMES[matches[0]]
        return col.lower().replace(" ", "_")

    def _smiles_to_inchikey(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            inchi_str = inchi.MolToInchi(mol)
            return inchi.InchiToInchiKey(inchi_str)
        except Exception:
            return None

    def clean(self):
        # Read input
        self.df = pd.read_csv(self.input_csv)

        # Fuzzy map headers
        self.df.rename(columns=lambda c: self._fuzzy_match(c), inplace=True)

        # Ensure required columns
        if "smiles" not in self.df.columns:
            raise ValueError("SMILES column missing — cannot proceed.")
        if "charge" not in self.df.columns:
            print("Charge column missing. Defaulting all charges to 0.")
            self.df["charge"] = 0

        # Generate InChIKey
        self.df.insert(0, "key_inchi", self.df["smiles"].apply(self._smiles_to_inchikey))

        # Ensure mol_name_iupac exists
        if "mol_name_iupac" not in self.df.columns:
            raise ValueError("No name column found to map into mol_name_iupac.")

        # Reorder columns
        first_cols = ["key_inchi", "smiles", "charge", "mol_name_iupac"]
        other_cols = [c for c in self.df.columns if c not in first_cols]
        self.df = self.df[first_cols + other_cols]

        # Save cleaned CSV
        self.df.to_csv(self.output_csv, index=False)
        print(f"Saved cleaned dataset to {self.output_csv}")

        # Save twin .inp file
        inp_file = self.output_csv.replace(".csv", ".inp")
        self._write_inp_file(inp_file)

        # Return the .inp file path so it can be chained
        return inp_file


    def _write_inp_file(self, inp_file):
        with open(inp_file, "w") as f:
            for _, row in self.df.iterrows():
                name = str(row["mol_name_iupac"])
                smiles = str(row["smiles"])
                charge = str(row["charge"])
                f.write(f"{name}\t{smiles}\t\t{charge}\n")
        print(f"Saved conformer generator input file to {inp_file}")
