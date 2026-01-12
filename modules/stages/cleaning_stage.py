import os
import pandas as pd
import difflib
from rdkit import Chem
from rdkit.Chem import inchi


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


class CleaningStage:
    """
    Data cleaning stage.

    Behaviour:
      - Each raw file in job.input_folder corresponds to one lookup_id
      - All cleaned rows are appended into one combined dataset
      - Produces a single output CSV: clean_dataset.csv
      - Updates job progress for each processed file
    """

    def __init__(self):
        pass

    # ------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # Stage entrypoint
    # ------------------------------------------------------------
    def run(self, job):
        """
        Execute the cleaning stage.

        Expected:
            job.parameters["input_folder"]  -> folder containing raw CSVs
            job.items                       -> stems of raw filenames
        """

        input_folder = job.parameters["input_folder"]

        # Output path
        output_csv = job.resolve_output_path("clean_dataset.csv")

        # Combined dataframe
        combined = []

        # --------------------------------------------------------
        # Process each raw file
        # --------------------------------------------------------
        for lookup_id in list(job.pending_items):

            raw_path = os.path.join(input_folder, f"{lookup_id}.csv")
            if not os.path.isfile(raw_path):
                job.write_pipeline_log(f"WARNING: Missing raw file for {lookup_id}: {raw_path}")
                job.update_progress(lookup_id)
                continue

            df = pd.read_csv(raw_path)

            # Fuzzy map headers
            df.rename(columns=lambda c: self._fuzzy_match(c), inplace=True)

            # Ensure required columns
            if "smiles" not in df.columns:
                raise ValueError(f"SMILES column missing in {raw_path}")

            if "charge" not in df.columns:
                df["charge"] = 0

            # Generate InChIKey
            df.insert(0, "key_inchi", df["smiles"].apply(self._smiles_to_inchikey))

            # Ensure mol_name_iupac exists
            if "mol_name_iupac" not in df.columns:
                raise ValueError(f"No name column found in {raw_path}")

            # Reorder columns
            first_cols = ["key_inchi", "smiles", "charge", "mol_name_iupac"]
            other_cols = [c for c in df.columns if c not in first_cols]
            df = df[first_cols + other_cols]

            # Append to combined dataset
            combined.append(df)

            # Mark progress
            job.update_progress(lookup_id)

        # --------------------------------------------------------
        # Combine and write output
        # --------------------------------------------------------
        if combined:
            final_df = pd.concat(combined, ignore_index=True)
            final_df.to_csv(output_csv, index=False)
        else:
            raise ValueError("No valid raw files were processed — cannot create cleaned dataset.")

        # --------------------------------------------------------
        # Return updated job
        # --------------------------------------------------------
        return job.with_outputs({"clean_dataset": output_csv})
