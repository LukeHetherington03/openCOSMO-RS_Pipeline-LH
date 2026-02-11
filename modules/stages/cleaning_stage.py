import os
import json
import difflib
import pandas as pd
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem import Descriptors
import hashlib

from modules.stages.base_stage import BaseStage


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

DEFAULT_METADATA = {
    "melting_temp": "N/A",
    "melting_temp_source": "N/A",
    "Hfus": "N/A",
    "Gfus_mode": "MyrdalYalkowsky",
    "formula_calcd": "",
    "anionic": False,
    "mol_name": "N/A",
    "mol_name_iupac": "N/A",
    "experimental_solubility_mol_frac": None
}

METADATA_VERSION = 1


class CleaningStage(BaseStage):
    """
    Cleans raw CSV files into a unified dataset and produces molecule metadata.

    Produces:
        cleaned.csv
        molecule_metadata/<inchi_key>.json
    """

    def execute(self):
        # Strict mode
        self.strict_mode = self.strict("cleaning")

        # Declare canonical output
        output_csv = self.set_stage_output("cleaned.csv")

        metadata_dir = self.output_path("molecule_metadata")
        os.makedirs(metadata_dir, exist_ok=True)

        self.warnings = {
            "missing_melting_temp": 0,
            "invalid_smiles": 0,
            "default_charge": 0,
        }

        # Resolve input CSVs
        csv_list = self._resolve_input_sources()

        # NEW: merge and store raw input
        self._write_combined_raw_input(csv_list)

        # Track items by filename stem
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

            df = self._reorder_columns(df)
            combined_rows.append(df)

            self.update_progress(lookup_id)

        if not combined_rows:
            self.fail("No valid raw files processed — cannot create cleaned dataset.")

        final_df = pd.concat(combined_rows, ignore_index=True)
        final_df.to_csv(output_csv, index=False)

        self._log_warning_summary()
        self.log(f"Clean dataset written to: {output_csv}")

    # ------------------------------------------------------------
    # Resolve input sources
    # ------------------------------------------------------------
    def _resolve_input_sources(self):
        csv_list = []

        # input_folder(s)
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

        # input_csv(s)
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

    def _write_combined_raw_input(self, csv_list):
        """
        Merge all input CSVs into a single raw_combined.csv file
        stored under the job's inputs directory.
        """
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

        # Write to job inputs directory
        out_path = self.input_path("raw_combined.csv")
        combined.to_csv(out_path, index=False)

        self.log(f"Merged raw input written to: {out_path}")


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
        df["inchi_key"] = df["smiles"].apply(self._smiles_to_inchikey)

        invalid = df["inchi_key"].isna()
        if invalid.any():
            bad = df.loc[invalid, "smiles"].tolist()
            self.log(f"[WARNING] Invalid SMILES in {lookup_id}: {bad}")
            self.warnings["invalid_smiles"] += len(bad)

        return df
    
    def _write_metadata_for_frame(self, df, metadata_dir, source_file):
        from rdkit.Chem import (
            Lipinski, rdMolDescriptors, Crippen, Descriptors, rdmolops
        )
        import hashlib

        # ------------------------------------------------------------
        # Functional group SMARTS patterns
        # ------------------------------------------------------------
        FG_SMARTS = {
            # Existing groups
            "is_amide": "[NX3][CX3](=[OX1])[#6]",
            "is_ester": "[CX3](=O)[OX2H0][#6]",
            "is_alcohol": "[OX2H][CX4]",
            "is_amine": "[NX3;H2,H1,H0;!$(NC=O)]",
            "is_carboxylic_acid": "C(=O)[OH]",
            "is_ketone": "[CX3](=O)[#6]",
            "is_aldehyde": "[CX3H1](=O)[#6]",
            "is_nitrile": "C#N",

            # Polymer-relevant groups
            "is_acrylate": "C=CC(=O)O[*]",
            "is_methacrylate": "C=C(C)C(=O)O[*]",
            "is_acrylamide": "C=CC(=O)N[*]",
            "is_methacrylamide": "C=C(C)C(=O)N[*]",
            "is_vinyl_ether": "C=CO[*]",
            "is_vinyl_ester": "C=COC(=O)[*]",
            "is_styrenic": "c1ccccc1C=C[*]",

            # Polar polymer groups
            "is_ether": "[#6]-O-[#6]",
            "is_carbonate": "O=C(O)O",
            "is_urethane": "O=C(N)O",
            "is_urea": "N-C(=O)-N",
            "is_sulfonamide": "S(=O)(=O)N",
            "is_sulfonate": "S(=O)(=O)O",
            "is_phosphate": "P(=O)(O)(O)",

            # Halogen flags
            "is_fluorinated": "[F]",
            "is_chlorinated": "[Cl]",
            "is_brominated": "[Br]",
        }

        FG_PATTERNS = {k: Chem.MolFromSmarts(v) for k, v in FG_SMARTS.items()}

        # ------------------------------------------------------------
        # Hash helper
        # ------------------------------------------------------------
        def _json_hash(obj):
            data = json.dumps(obj, sort_keys=True).encode("utf-8")
            return hashlib.sha256(data).hexdigest()

        pipeline_version = self.config.get("pipeline_version", "unknown")
        timestamp = datetime.utcnow().isoformat() + "Z"
        overwrite = self.parameters.get("overwrite_metadata", False)

        # ------------------------------------------------------------
        # Process each molecule
        # ------------------------------------------------------------
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

            # ------------------------------------------------------------
            # Physchem descriptors
            # ------------------------------------------------------------
            if mol is not None:
                rotatable_bonds = Lipinski.NumRotatableBonds(mol)
                molecular_weight = rdMolDescriptors.CalcExactMolWt(mol)
                heavy_atom_count = mol.GetNumHeavyAtoms()
                hbond_donors = Lipinski.NumHDonors(mol)
                hbond_acceptors = Lipinski.NumHAcceptors(mol)
                tpsa = rdMolDescriptors.CalcTPSA(mol)
                logp = Crippen.MolLogP(mol)
                aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

                # ------------------------------------------------------------
                # Structural counts
                # ------------------------------------------------------------
                N_count = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N")
                O_count = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "O")
                S_count = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "S")
                halogen_count = sum(1 for a in mol.GetAtoms() if a.GetSymbol() in ["F", "Cl", "Br", "I"])

                ring_count = mol.GetRingInfo().NumRings()
                heterocycle_count = sum(
                    1 for ring in mol.GetRingInfo().AtomRings()
                    if any(mol.GetAtomWithIdx(i).GetAtomicNum() not in (6, 1) for i in ring)
                )

                double_bond_count = sum(1 for b in mol.GetBonds() if b.GetBondType().name == "DOUBLE")
                triple_bond_count = sum(1 for b in mol.GetBonds() if b.GetBondType().name == "TRIPLE")

                # ------------------------------------------------------------
                # Functional group flags & counts
                # ------------------------------------------------------------
                fg_flags = {}
                fg_counts = {}

                for name, patt in FG_PATTERNS.items():
                    matches = mol.GetSubstructMatches(patt)
                    fg_flags[name] = len(matches) > 0
                    fg_counts[name.replace("is_", "") + "_count"] = len(matches)

                # ------------------------------------------------------------
                # Global shape descriptors
                # ------------------------------------------------------------
                Fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
                Bertz = Descriptors.BertzCT(mol)
                molar_refractivity = Crippen.MolMR(mol)

            else:
                # All None if SMILES invalid
                rotatable_bonds = molecular_weight = heavy_atom_count = None
                hbond_donors = hbond_acceptors = None
                tpsa = logp = aromatic_rings = None
                N_count = O_count = S_count = halogen_count = None
                ring_count = heterocycle_count = None
                double_bond_count = triple_bond_count = None
                fg_flags = {}
                fg_counts = {}
                Fsp3 = Bertz = molar_refractivity = None

            # ------------------------------------------------------------
            # Build metadata block
            # ------------------------------------------------------------
            meta = {
                "metadata_version": 3,
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

                # Structural counts
                "N_count": N_count,
                "O_count": O_count,
                "S_count": S_count,
                "halogen_count": halogen_count,
                "ring_count": ring_count,
                "heterocycle_count": heterocycle_count,
                "double_bond_count": double_bond_count,
                "triple_bond_count": triple_bond_count,

                # Functional groups
                **fg_flags,
                **fg_counts,

                # Global shape
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

            # Optional metadata fields
            for key, default in DEFAULT_METADATA.items():
                if key in ("mol_name", "mol_name_iupac"):
                    continue
                value = row.get(key, default)
                if isinstance(value, float) and pd.isna(value):
                    value = default
                meta[key] = value

            # Warn if melting temp missing
            if meta["melting_temp"] == "N/A":
                self.log(
                    f"[WARNING] melting_temp not provided for {inchi_key}. "
                    "This may reduce solubility prediction accuracy downstream."
                )
                self.warnings["missing_melting_temp"] += 1

            # Serialise
            meta_json = json.dumps(meta, indent=2)
            meta_hash = _json_hash(meta)

            # ------------------------------------------------------------
            # Write local metadata
            # ------------------------------------------------------------
            local_path = os.path.join(metadata_dir, f"{inchi_key}.json")
            with open(local_path, "w") as f:
                f.write(meta_json)

            # ------------------------------------------------------------
            # Write global metadata (with overwrite option)
            # ------------------------------------------------------------
            global_meta_dir = self.config.get("constant_files", {}).get("metadata_dir")
            if not global_meta_dir:
                self.log("[WARNING] No global metadata_dir defined in config; skipping global write.")
                continue

            os.makedirs(global_meta_dir, exist_ok=True)
            global_path = os.path.join(global_meta_dir, f"{inchi_key}.json")

            write_global = True

            if not overwrite:
                if os.path.exists(global_path):
                    try:
                        with open(global_path) as f:
                            existing = json.load(f)
                        if _json_hash(existing) == meta_hash:
                            write_global = False
                    except Exception:
                        write_global = True

            if write_global:
                with open(global_path, "w") as f:
                    f.write(meta_json)
                if overwrite:
                    self.log(f"[INFO] Global metadata overwritten: {global_path}")
                else:
                    self.log(f"[INFO] Global metadata updated: {global_path}")
            else:
                self.log(f"[INFO] Global metadata unchanged: {global_path}")


    # ------------------------------------------------------------
    # Column ordering
    # ------------------------------------------------------------
    def _reorder_columns(self, df):
        first = ["inchi_key", "smiles", "charge"]
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
