import os
import json
import difflib
import pandas as pd
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem import Descriptors
from modules.stages.base_stage import BaseStage


# ---------------------------------------------------------------------------
# Column name normalisation
# ---------------------------------------------------------------------------

STANDARD_NAMES = {
    "experimental solubility /mol frac": "experimental_solubility_mol_frac",
    "melting point": "melting_temp",
    "melting point source": "melting_temp_source",
    "mol_id": "mol_id",
    "pubchem_sid": "pubchem_sid",
    "mol_name": "mol_name_iupac",
    "name": "mol_name",
    "smiles": "smiles",
    "formula_calcd": "formula_calcd",
    "anionic": "anionic",
    "hfus": "Hfus",
    "gfus_mode": "Gfus_mode",
}

# --- NEW: AqSol format support (non-destructive extension) ---
STANDARD_NAMES.update({
    "solubility": "experimental_solubility_mol_frac",
    "predicted_solubility": "aqsol_predicted_solubility",
})


DEFAULT_METADATA = {
    "melting_temp": "N/A",
    "melting_temp_source": "N/A",
    "Hfus": "N/A",
    "Gfus_mode": "MyrdalYalkowsky",
    "formula_calcd": "",
    "anionic": False,
    "mol_name": "N/A",
    "mol_name_iupac": "N/A",
    "experimental_solubility_mol_frac": None,
}

METADATA_VERSION = 3


class CleaningStage(BaseStage):
    """
    Cleans raw CSV files into a unified dataset and produces molecule metadata.

    Supports:
        - Legacy experimental CSVs
        - AqSol ML CSVs (Name, Solubility, SMILES, Predicted_Solubility)

    Produces:
        cleaned.csv
        molecule_metadata/<inchi_key>.json
    """

    # ----------------------------------------------------------------------
    # Execute
    # ----------------------------------------------------------------------
    def execute(self):
        self.strict_mode = self.strict("cleaning")

        output_csv = self.set_stage_output("cleaned.csv")
        metadata_dir = self.output_path("molecule_metadata")
        os.makedirs(metadata_dir, exist_ok=True)

        self.warnings = {
            "missing_melting_temp": 0,
            "invalid_smiles": 0,
            "default_charge": 0,
        }

        csv_list = self._resolve_input_sources()
        self._write_combined_raw_input(csv_list)

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

            source_file = self._find_source_file(lookup_id, csv_list)
            self._write_metadata_for_frame(df, metadata_dir, source_file)

            df = self._minimal_cleaned_frame(df)
            combined_rows.append(df)

            self.update_progress(lookup_id)

        if not combined_rows:
            self.fail("No valid raw files processed — cannot create cleaned dataset.")

        final_df = pd.concat(combined_rows, ignore_index=True)
        final_df.to_csv(output_csv, index=False)

        self._log_warning_summary()
        self.log(f"Clean dataset written to: {output_csv}")

    # ----------------------------------------------------------------------
    # Resolve input sources
    # ----------------------------------------------------------------------
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
                self.require_file(p, "input_csv")
                csv_list.append(p)

        if not csv_list:
            self.fail("No CSV files detected. Provide input_folder(s) or input_csv(s).")

        csv_list = sorted(set(csv_list))
        self.log(f"Resolved {len(csv_list)} CSV files for cleaning")
        return csv_list

    # ----------------------------------------------------------------------
    # Write combined raw input
    # ----------------------------------------------------------------------
    def _write_combined_raw_input(self, csv_list):
        frames = []
        for path in csv_list:
            try:
                df = pd.read_csv(path)
                df["__source_file"] = os.path.abspath(path)
                frames.append(df)
            except Exception as e:
                self.log(f"[WARNING] Failed to read {path}: {e}")

        if not frames:
            self.fail("No valid CSVs could be read for raw input merge.")

        combined = pd.concat(frames, ignore_index=True)
        out_path = self.input_path("raw_combined.csv")
        combined.to_csv(out_path, index=False)
        self.log(f"Merged raw input written to: {out_path}")

    # ----------------------------------------------------------------------
    # Load CSV
    # ----------------------------------------------------------------------
    def _load_raw_csv(self, lookup_id, csv_list):
        matches = [
            p for p in csv_list
            if os.path.splitext(os.path.basename(p))[0] == lookup_id
        ]
        if not matches:
            self.log(f"[WARNING] No CSV found for lookup_id {lookup_id}")
            return None

        try:
            return pd.read_csv(matches[0])
        except Exception as e:
            self.log(f"[ERROR] Failed to read {matches[0]}: {e}")
            return None

    def _find_source_file(self, lookup_id, csv_list):
        for p in csv_list:
            if os.path.splitext(os.path.basename(p))[0] == lookup_id:
                return os.path.abspath(p)
        return None

    # ----------------------------------------------------------------------
    # Header standardisation
    # ----------------------------------------------------------------------
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

    # ----------------------------------------------------------------------
    # Validation
    # ----------------------------------------------------------------------
    def _validate_required_fields(self, df, lookup_id):
        has_smiles = "smiles" in df.columns
        has_name = ("mol_name" in df.columns) or ("mol_name_iupac" in df.columns)
        has_sol = "experimental_solubility_mol_frac" in df.columns

        if not (has_smiles and has_name and has_sol):
            msg = (
                f"Missing required fields in {lookup_id}.csv — "
                f"need SMILES, Name, Solubility"
            )
            self.log(f"[ERROR] {msg}")
            if self.strict_mode:
                self.fail(msg)
            return False

        if "charge" not in df.columns:
            df["charge"] = 0
            self.warnings["default_charge"] += 1

        return True

    # ----------------------------------------------------------------------
    # SMILES canonicalisation
    # ----------------------------------------------------------------------
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

    # ----------------------------------------------------------------------
    # InChIKey generation
    # ----------------------------------------------------------------------
    def _smiles_to_inchikey(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return inchi.InchiToInchiKey(inchi.MolToInchi(mol))
        except Exception:
            return None

    def _generate_inchikeys(self, df, lookup_id):
        df["inchi_key"] = df["smiles"].apply(self._smiles_to_inchikey)

        invalid = df["inchi_key"].isna()
        if invalid.any():
            bad = df.loc[invalid, "smiles"].tolist()
            self.log(f"[WARNING] Invalid SMILES in {lookup_id}: {bad}")
            self.warnings["invalid_smiles"] += len(bad)

        return df
    
    # ----------------------------------------------------------------------
    # Metadata writer
    # ----------------------------------------------------------------------
    def _write_metadata_for_frame(self, df, metadata_dir, source_file):
        from rdkit.Chem import Lipinski, rdMolDescriptors, Crippen
        from modules.utils.molecule_utils import MoleculeUtils

        pipeline_version = self.config.get("pipeline_version", "unknown")
        timestamp = datetime.utcnow().isoformat() + "Z"
        overwrite = self.parameters.get("overwrite_metadata", False)

        for inchi_key, group in df.groupby("inchi_key"):
            if not inchi_key:
                continue

            row = group.iloc[0]
            smiles = row.get("smiles", "")
            mol = Chem.MolFromSmiles(smiles)

            charge_val = int(row.get("charge", 0))
            multiplicity = 1 if charge_val == 0 else 2

            mol_name = row.get("mol_name", DEFAULT_METADATA["mol_name"])
            mol_name_iupac = row.get("mol_name_iupac", DEFAULT_METADATA["mol_name_iupac"])

            # Physchem descriptors
            if mol is not None:
                rotatable_bonds = Lipinski.NumRotatableBonds(mol)
                molecular_weight = rdMolDescriptors.CalcExactMolWt(mol)
                heavy_atom_count = mol.GetNumHeavyAtoms()
                hbond_donors = Lipinski.NumHDonors(mol)
                hbond_acceptors = Lipinski.NumHAcceptors(mol)
                tpsa = rdMolDescriptors.CalcTPSA(mol)
                logp = Crippen.MolLogP(mol)
                aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
                Fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
                Bertz = Descriptors.BertzCT(mol)
                molar_refractivity = Crippen.MolMR(mol)
            else:
                rotatable_bonds = molecular_weight = heavy_atom_count = None
                hbond_donors = hbond_acceptors = None
                tpsa = logp = aromatic_rings = None
                Fsp3 = Bertz = molar_refractivity = None

            meta = {
                "metadata_version": METADATA_VERSION,
                "inchi_key": inchi_key,
                "smiles": smiles,
                "mol_name": mol_name,
                "mol_name_iupac": mol_name_iupac,
                "charge": charge_val,
                "multiplicity": multiplicity,

                # Physchem
                "rotatable_bonds": rotatable_bonds,
                "molecular_weight": molecular_weight,
                "heavy_atom_count": heavy_atom_count,
                "hbond_donors": hbond_donors,
                "hbond_acceptors": hbond_acceptors,
                "tpsa": tpsa,
                "logp": logp,
                "aromatic_rings": aromatic_rings,
                "Fsp3": Fsp3,
                "Bertz_complexity": Bertz,
                "molar_refractivity": molar_refractivity,

                "provenance": {
                    "source_file": source_file,
                    "cleaning_timestamp": timestamp,
                    "pipeline_version": pipeline_version,
                    "generated_by_request": self.job.request_id,
                    "generated_by_job": self.job.job_id,
                },
            }

            # Add solubility fields (legacy + AqSol)
            for field in ["experimental_solubility_mol_frac", "aqsol_predicted_solubility"]:
                if field in row and not pd.isna(row[field]):
                    meta[field] = row[field]

            # Fill defaults
            for key, default in DEFAULT_METADATA.items():
                if key not in meta:
                    meta[key] = default

            # --------------------------------------------------------------
            # NEW: PubChem fallback for melting_temp
            # --------------------------------------------------------------
            mp_info = MoleculeUtils.get_melting_point(inchi_key)

            # Only fill if missing or "N/A"
            if meta["melting_temp"] in (None, "N/A", "", "nan"):
                meta["melting_temp"] = mp_info["melting_temp"]
                meta["melting_temp_source"] = mp_info["melting_temp_source"]

                if mp_info["melting_temp_source"] == "pubchem":
                    self.log(f"[MP] PubChem fallback used for {inchi_key}: {mp_info['melting_temp']} °C")
                else:
                    self.log(f"[MP] No melting point found for {inchi_key} (missing).")
            else:
                # Existing value preserved
                meta["melting_temp_source"] = "existing"

            # Warn if still missing
            if meta["melting_temp"] == "N/A":
                self.log(
                    f"[WARNING] melting_temp not provided for {inchi_key}. "
                    "This may reduce solubility prediction accuracy downstream."
                )
                self.warnings["missing_melting_temp"] += 1

            # Write local metadata
            local_path = os.path.join(metadata_dir, f"{inchi_key}.json")
            with open(local_path, "w") as f:
                f.write(json.dumps(meta, indent=2))

            # Write global metadata
            global_meta_dir = self.config.get("constant_files", {}).get("metadata_dir")
            if not global_meta_dir:
                self.log("[WARNING] No global metadata_dir defined in config; skipping global write.")
                continue

            os.makedirs(global_meta_dir, exist_ok=True)
            global_path = os.path.join(global_meta_dir, f"{inchi_key}.json")

            if overwrite or not os.path.exists(global_path):
                with open(global_path, "w") as f:
                    f.write(json.dumps(meta, indent=2))

    # ----------------------------------------------------------------------
    # Minimal cleaned CSV
    # ----------------------------------------------------------------------
    def _minimal_cleaned_frame(self, df):
        keep = [
            "inchi_key",
            "smiles",
            "charge",
            "mol_name",
            "mol_name_iupac",
            "experimental_solubility_mol_frac",
            "aqsol_predicted_solubility",
        ]
        return df[[c for c in keep if c in df.columns]]

    # ----------------------------------------------------------------------
    # Warning summary
    # ----------------------------------------------------------------------
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
