import os
import pandas as pd
import difflib
from rdkit import Chem
from rdkit.Chem import inchi

from modules.stages.base_stage import BaseStage


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
    "anionic": "anionic",
}


class CleaningStage(BaseStage):
    """
    Cleans raw CSV files into a unified dataset.

    Behaviour:
      - Accepts input_folder(s) OR input_csv(s)
      - Collects all CSVs into a unified list
      - Standardises column names
      - Validates SMILES and generates InChIKey
      - Produces clean_dataset.csv for the next stage
    """

    # ------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------
    def execute(self):
        self.log_header("Starting Cleaning Stage")

        # ------------------------------------------------------------
        # Resolve input sources
        # ------------------------------------------------------------
        csv_list = self._resolve_input_sources()

        # Derive lookup_ids
        lookup_ids = [
            os.path.splitext(os.path.basename(p))[0]
            for p in csv_list
        ]

        output_csv = self.output_path("clean_dataset.csv")

        # Initialise items
        self.set_items(lookup_ids)

        combined_rows = []

        for lookup_id in list(self.job.pending_items):
            df = self._load_raw_csv(lookup_id, csv_list)
            if df is None:
                self.update_progress(lookup_id, success=False)
                continue

            df = self._standardise_headers(df, lookup_id)
            if not self._validate_required_fields(df, lookup_id):
                self.update_progress(lookup_id, success=False)
                continue

            df = self._generate_inchikeys(df, lookup_id)
            df = self._reorder_columns(df)

            combined_rows.append(df)
            self.update_progress(lookup_id)

        if not combined_rows:
            self.fail("No valid raw files processed — cannot create cleaned dataset.")

        final_df = pd.concat(combined_rows, ignore_index=True)
        final_df.to_csv(output_csv, index=False)

        self.log(f"Clean dataset written to: {output_csv}")
        self.job.mark_complete()

    # ------------------------------------------------------------
    # Resolve input sources
    # ------------------------------------------------------------
    def _resolve_input_sources(self):
        """
        Accepts:
          - input_folder: "path/"
          - input_folders: ["path1/", "path2/"]
          - input_csv: "file.csv"
          - input_csv: ["file1.csv", "file2.csv"]

        Returns:
          - csv_list: list of absolute CSV paths
        """

        csv_list = []

        # --- Folder mode ---
        if "input_folder" in self.parameters:
            folders = [self.parameters["input_folder"]]
        elif "input_folders" in self.parameters:
            folders = self.parameters["input_folders"]
            if isinstance(folders, str):
                folders = [folders]
        else:
            folders = []

        # Collect CSVs from folders
        for folder in folders:
            if not os.path.isdir(folder):
                self.fail(f"input_folder is not a directory: {folder}")

            for fname in os.listdir(folder):
                if fname.lower().endswith(".csv"):
                    csv_list.append(os.path.join(folder, fname))

        # --- CSV list mode ---
        if "input_csv" in self.parameters:
            raw = self.parameters["input_csv"]
            if isinstance(raw, str):
                raw = [raw]

            for p in raw:
                if not os.path.isfile(p):
                    self.fail(f"input_csv file does not exist: {p}")
                csv_list.append(p)

        # Final validation
        if not csv_list:
            self.fail("No CSV files detected. Provide input_folder(s) or input_csv(s).")

        csv_list = sorted(set(csv_list))
        self.log(f"Resolved {len(csv_list)} CSV files for cleaning")

        return csv_list

    # ------------------------------------------------------------
    # Load CSV
    # ------------------------------------------------------------
    def _load_raw_csv(self, lookup_id, csv_list):
        matches = [
            p for p in csv_list
            if os.path.splitext(os.path.basename(p))[0] == lookup_id
        ]

        if not matches:
            self.log(f"[WARNING] No CSV found for lookup_id {lookup_id}")
            return None

        path = matches[0]

        try:
            return pd.read_csv(path)
        except Exception as e:
            self.log(f"[ERROR] Failed to read {path}: {e}")
            return None

    # ------------------------------------------------------------
    # Helpers: header standardisation
    # ------------------------------------------------------------
    def _normalise_header(self, col):
        return col.strip().lower().replace(" ", "_")

    def _fuzzy_match(self, col):
        col_norm = self._normalise_header(col)
        matches = difflib.get_close_matches(
            col_norm, STANDARD_NAMES.keys(), n=1, cutoff=0.65
        )
        if matches:
            return STANDARD_NAMES[matches[0]]
        return col_norm

    def _standardise_headers(self, df, lookup_id):
        original_cols = df.columns.tolist()
        df = df.rename(columns=lambda c: self._fuzzy_match(c))

        for old, new in zip(original_cols, df.columns):
            if old != new and new not in STANDARD_NAMES.values():
                self.log(f"[INFO] Column '{old}' normalised to '{new}'")

        return df

    # ------------------------------------------------------------
    # Helpers: validation
    # ------------------------------------------------------------
    def _validate_required_fields(self, df, lookup_id):
        if "smiles" not in df.columns:
            self.log(f"[ERROR] SMILES column missing in {lookup_id}.csv")
            return False

        if "mol_name_iupac" not in df.columns:
            self.log(f"[ERROR] mol_name_iupac column missing in {lookup_id}.csv")
            return False

        if "charge" not in df.columns:
            df["charge"] = 0

        return True

    # ------------------------------------------------------------
    # Helpers: InChIKey generation
    # ------------------------------------------------------------
    def _smiles_to_inchikey(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return inchi.InchiToInchiKey(inchi.MolToInchi(mol))
        except Exception:
            return None

    def _generate_inchikeys(self, df, lookup_id):
        df["key_inchi"] = df["smiles"].apply(self._smiles_to_inchikey)

        invalid = df["key_inchi"].isna()
        if invalid.any():
            bad = df.loc[invalid, "smiles"].tolist()
            self.log(f"[WARNING] Invalid SMILES in {lookup_id}: {bad}")

        return df

    # ------------------------------------------------------------
    # Helpers: column ordering
    # ------------------------------------------------------------
    def _reorder_columns(self, df):
        first = ["key_inchi", "smiles", "charge", "mol_name_iupac"]
        others = [c for c in df.columns if c not in first]
        return df[first + others]
