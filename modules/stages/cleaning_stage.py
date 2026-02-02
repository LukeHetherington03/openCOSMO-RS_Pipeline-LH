import os
import json
import difflib
import pandas as pd
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import inchi

from modules.stages.base_stage import BaseStage


STANDARD_NAMES = {
    "experimental solubility /mol frac": "experimental_solubility_mol_frac",
    "melting point": "melting_temp",
    "melting point source": "melting_temp_source",
    "mol_id": "mol_id",
    "pubchem_sid": "pubchem_sid",
    "mol_name": "mol_name_iupac",
    "name": "mol_name",  # original/AKA name
    "smiles": "smiles",
    "formula_calcd": "formula_calcd",
    "anionic": "anionic",
    "hfus": "Hfus",
    "gfus_mode": "Gfus_mode",
}


DEFAULT_METADATA = {
    "melting_temp": "N/A",
    "melting_temp_source": "N/A",
    "Hfus": "N/A",
    "Gfus_mode": "MyrdalYalkowsky",
    "formula_calcd": "",
    "anionic": False,
    "mol_name": "N/A",
    "mol_name_iupac": "N/A",
}

METADATA_VERSION = 1


class CleaningStage(BaseStage):
    """
    Cleans raw CSV files into a unified dataset and produces molecule metadata.

    Behaviour:
      - Accepts input_folder(s) OR input_csv(s)
      - Collects all CSVs into a unified list
      - Standardises column names
      - Validates SMILES and generates InChIKey
      - Canonicalises SMILES
      - Produces:
            clean_dataset.csv
            molecule_metadata/<inchi_key>.json
    """

    # ------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------
    def execute(self):
        self.log_header("Starting Cleaning Stage")

        self.strict_mode = bool(self.config.get("cleaning", {}).get("strict", False))
        self.warnings = {
            "missing_melting_temp": 0,
            "invalid_smiles": 0,
            "default_charge": 0,
        }

        csv_list = self._resolve_input_sources()
        output_csv = self.output_path("clean_dataset.csv")
        metadata_dir = self.output_path("molecule_metadata")
        os.makedirs(metadata_dir, exist_ok=True)

        # We track items by source file stem, but do not expose lookup_id downstream
        lookup_ids = [
            os.path.splitext(os.path.basename(p))[0]
            for p in csv_list
        ]
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

            df = self._canonicalise_smiles(df, lookup_id)
            df = self._generate_inchikeys(df, lookup_id)

            # Write per-molecule metadata (molecule-level, no lookup_id/conf_num)
            source_file = self._find_source_file(lookup_id, csv_list)
            self._write_metadata_for_frame(df, metadata_dir, source_file)

            df = self._reorder_columns(df)
            combined_rows.append(df)
            self.update_progress(lookup_id)

        if not combined_rows:
            self.fail("No valid raw files processed — cannot create cleaned dataset.")

        final_df = pd.concat(combined_rows, ignore_index=True)
        final_df.to_csv(output_csv, index=False)

        self._log_warning_summary()
        self.log(f"Clean dataset written to: {output_csv}")
        self.job.mark_complete()

    # ------------------------------------------------------------
    # Resolve input sources
    # ------------------------------------------------------------
    def _resolve_input_sources(self):
        csv_list = []

        if "input_folder" in self.parameters:
            folders = [self.parameters["input_folder"]]
        elif "input_folders" in self.parameters:
            folders = self.parameters["input_folders"]
            if isinstance(folders, str):
                folders = [folders]
        else:
            folders = []

        for folder in folders:
            if not os.path.isdir(folder):
                self.fail(f"input_folder is not a directory: {folder}")

            for fname in os.listdir(folder):
                if fname.lower().endswith(".csv"):
                    csv_list.append(os.path.join(folder, fname))

        if "input_csv" in self.parameters:
            raw = self.parameters["input_csv"]
            if isinstance(raw, str):
                raw = [raw]

            for p in raw:
                if not os.path.isfile(p):
                    self.fail(f"input_csv file does not exist: {p}")
                csv_list.append(p)

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

    def _find_source_file(self, lookup_id, csv_list):
        for p in csv_list:
            if os.path.splitext(os.path.basename(p))[0] == lookup_id:
                return os.path.abspath(p)
        return None

    # ------------------------------------------------------------
    # Header standardisation
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
    # Validation
    # ------------------------------------------------------------
    def _validate_required_fields(self, df, lookup_id):
        if "smiles" not in df.columns:
            msg = f"SMILES column missing in {lookup_id}.csv"
            self.log(f"[ERROR] {msg}")
            if self.strict_mode:
                self.fail(msg)
            return False

        has_name = ("mol_name" in df.columns) or ("mol_name_iupac" in df.columns)
        if not has_name:
            msg = f"No name column (Name/mol_name_iupac) in {lookup_id}.csv"
            self.log(f"[ERROR] {msg}")
            if self.strict_mode:
                self.fail(msg)
            return False

        if "charge" not in df.columns:
            self.log(f"[INFO] 'charge' column missing in {lookup_id}.csv — defaulting to 0")
            df["charge"] = 0
            self.warnings["default_charge"] += 1

        return True

    # ------------------------------------------------------------
    # SMILES canonicalisation
    # ------------------------------------------------------------
    def _canonicalise_smiles(self, df, lookup_id):
        def _canon(smi):
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    return smi
                return Chem.MolToSmiles(mol, canonical=True)
            except Exception:
                return smi

        df["smiles"] = df["smiles"].astype(str).apply(_canon)
        return df

    # ------------------------------------------------------------
    # InChIKey generation
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
            self.warnings["invalid_smiles"] += len(bad)

        return df

    # ------------------------------------------------------------
    # Metadata writing
    # ------------------------------------------------------------
    def _write_metadata_for_frame(self, df, metadata_dir, source_file):
        """
        For each unique InChIKey in the dataframe, write a molecule-level
        metadata JSON file. Metadata is molecule-level, keyed by InChIKey.
        No lookup_id or conf_num is stored here.
        """
        from rdkit.Chem import Lipinski, rdMolDescriptors, Crippen

        pipeline_version = self.config.get("pipeline_version", "unknown")
        timestamp = datetime.utcnow().isoformat() + "Z"

        for inchi_key, group in df.groupby("key_inchi"):
            if not inchi_key:
                continue

            row = group.iloc[0]

            smiles = row.get("smiles", "")
            mol = Chem.MolFromSmiles(smiles)

            # Basic fields
            charge_val = int(row.get("charge", 0))
            multiplicity = 1 if charge_val == 0 else 2

            mol_name = row.get("mol_name", DEFAULT_METADATA["mol_name"])
            mol_name_iupac = row.get("mol_name_iupac", DEFAULT_METADATA["mol_name_iupac"])

            # --- NEW DESCRIPTORS ---
            if mol is not None:
                rotatable_bonds = Lipinski.NumRotatableBonds(mol)
                molecular_weight = rdMolDescriptors.CalcExactMolWt(mol)
                heavy_atom_count = mol.GetNumHeavyAtoms()
                hbond_donors = Lipinski.NumHDonors(mol)
                hbond_acceptors = Lipinski.NumHAcceptors(mol)
                tpsa = rdMolDescriptors.CalcTPSA(mol)
                logp = Crippen.MolLogP(mol)
                aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
            else:
                # Fallbacks if SMILES invalid
                rotatable_bonds = None
                molecular_weight = None
                heavy_atom_count = None
                hbond_donors = None
                hbond_acceptors = None
                tpsa = None
                logp = None
                aromatic_rings = None

            meta = {
                "metadata_version": METADATA_VERSION,
                "inchi_key": inchi_key,
                "smiles": smiles,
                "mol_name": mol_name,
                "mol_name_iupac": mol_name_iupac,
                "charge": charge_val,
                "multiplicity": multiplicity,

                # --- NEW DESCRIPTORS ---
                "rotatable_bonds": rotatable_bonds,
                "molecular_weight": molecular_weight,
                "heavy_atom_count": heavy_atom_count,
                "hbond_donors": hbond_donors,
                "hbond_acceptors": hbond_acceptors,
                "tpsa": tpsa,
                "logp": logp,
                "aromatic_rings": aromatic_rings,

                "provenance": {
                    "source_file": source_file,
                    "cleaning_timestamp": timestamp,
                    "pipeline_version": pipeline_version,
                },
            }

            # Add default metadata fields (melting point, etc.)
            for key, default in DEFAULT_METADATA.items():
                if key in ("mol_name", "mol_name_iupac"):
                    continue
                value = row.get(key, default)
                if isinstance(value, float) and pd.isna(value):
                    value = default
                meta[key] = value

            if meta["melting_temp"] == "N/A":
                self.log(
                    f"[WARNING] melting_temp not provided for {inchi_key}. "
                    "This may reduce solubility prediction accuracy downstream."
                )
                self.warnings["missing_melting_temp"] += 1

            out_path = os.path.join(metadata_dir, f"{inchi_key}.json")
            with open(out_path, "w") as f:
                json.dump(meta, f, indent=2)


    # ------------------------------------------------------------
    # Column ordering
    # ------------------------------------------------------------
    def _reorder_columns(self, df):
        first = ["key_inchi", "smiles", "charge"]
        if "mol_name" in df.columns:
            first.append("mol_name")
        if "mol_name_iupac" in df.columns:
            first.append("mol_name_iupac")

        others = [c for c in df.columns if c not in first]
        return df[first + others]

    # ------------------------------------------------------------
    # Warning summary
    # ------------------------------------------------------------
    def _log_warning_summary(self):
        total = sum(self.warnings.values())
        if total == 0:
            self.log("Cleaning completed with no warnings.")
            return

        self.log("Cleaning completed with warnings:")
        if self.warnings["missing_melting_temp"] > 0:
            self.log(
                f" - {self.warnings['missing_melting_temp']} molecule(s) missing melting_temp "
                "(may reduce solubility prediction accuracy)."
            )
        if self.warnings["invalid_smiles"] > 0:
            self.log(
                f" - {self.warnings['invalid_smiles']} invalid SMILES encountered (rows retained but flagged)."
            )
        if self.warnings["default_charge"] > 0:
            self.log(
                f" - 'charge' defaulted to 0 for {self.warnings['default_charge']} file(s)."
            )
