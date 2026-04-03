import os
import json
import difflib
import pandas as pd
import numpy as np
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import inchi, Descriptors
from modules.stages.base_stage import BaseStage


# ─────────────────────────────────────────────────────────────────────────────
# Mode constants
# ─────────────────────────────────────────────────────────────────────────────

MODE_PIPELINE = "pipeline"
# Output of filter.py.  Already has canonical column names, solubility already
# in mol fraction, melting_temp already in Kelvin.  No conversion needed.
# Detected by: 'experimental_solubility_mol_frac' column present
#           OR 'melting_temp' column present with median value > 100 (Kelvin)

MODE_CUSTOM = "custom"
# Any hand-curated or third-party CSV.
# Required  : SMILES + a name column.
# Optional  : solubility (assumed mol frac unless `solubility_unit` param set)
#             melting point (any column whose name contains "melting" or "mp")
# If melting point is absent the stage logs clearly and falls back to PubChem.


# ─────────────────────────────────────────────────────────────────────────────
# Solubility unit options  (set via parameter `solubility_unit`)
# ─────────────────────────────────────────────────────────────────────────────

UNIT_MOL_FRAC  = "mol_frac"   # 0–1, dimensionless          (default)
UNIT_MOL_PER_L = "mol_per_L"  # linear mol/L concentration
UNIT_LOGS      = "logS"       # log₁₀(mol/L), may be negative

WATER_MOLARITY = 55.51        # mol/L at 25 °C, for mol/L → mol frac conversion


# ─────────────────────────────────────────────────────────────────────────────
# Melting point plausibility window (Kelvin)
# Values outside this range are flagged and not used — the molecule is kept;
# only the MP datum is marked suspect and PubChem is tried as fallback.
# ─────────────────────────────────────────────────────────────────────────────

MP_MIN_K = 100.0    # ≈ −173 °C
MP_MAX_K = 1200.0   # ≈  927 °C


# ─────────────────────────────────────────────────────────────────────────────
# Custom-mode fuzzy column matching
# ─────────────────────────────────────────────────────────────────────────────
# Non-overlapping alias lists — each raw column maps to at most one target.

_CUSTOM_ALIASES = {
    "smiles": [
        "smiles", "smi", "canonical_smiles", "structure",
    ],
    "mol_name": [
        "name", "mol_name", "compound", "compound_name",
        "molecule", "chemical_name", "substance_name",
        "mol_name_iupac", "iupac_name", "iupac",
    ],
    "sol_raw": [
        "solubility", "experimental_solubility_/mol_frac",
        "experimental_solubility", "sol", "concentration",
    ],
    "melting_temp_kelvin": [
        "melting_point", "melting_temp", "meltingpoint", "mp",
        "tm", "melting", "mpt", "melting_point_k", "melting_temp_k",
    ],
    "melting_temp_source": [
        "melting_point_source", "melting_temp_source",
        "mp_source", "source", "mp_reference",
    ],
}

DEFAULT_METADATA = {
    "melting_temp":                     "N/A",
    "melting_temp_source":              "N/A",
    "Hfus":                             "N/A",
    "Gfus_mode":                        "MyrdalYalkowsky",
    "formula_calcd":                    "",
    "anionic":                          False,
    "mol_name":                         "N/A",
    "mol_name_iupac":                   "N/A",
    "experimental_solubility_mol_frac": None,
}

METADATA_VERSION = 4

# ─────────────────────────────────────────────────────────────────────────────
# Functional group SMARTS patterns
# Counts are stored as functional_groups: {name: int} in metadata.
# Only groups with ≥5% prevalence in a typical drug-like / solubility dataset
# are included — rarer groups produce near-zero columns that L1 regularisation
# zeroes out anyway, and pollute the metadata.
# ─────────────────────────────────────────────────────────────────────────────
FUNCTIONAL_GROUP_SMARTS = {
    # Acids / esters / carbonyl
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "ester":           "[CX3](=O)[OX2][#6]",
    "amide":           "[CX3](=O)[NX3]",
    "ketone":          "[CX3H0;!R](=O)([#6])[#6]",
    "aldehyde":        "[CX3H1](=O)",
    "lactam":          "[NX3R][CX3R](=O)",
    "lactone":         "[OX2R][CX3R](=O)",
    # Hydroxyl
    "alcohol":         "[OX2H;!$(OC=O);!$(Oc1ccccc1)]",
    "phenol":          "[OX2H]c1ccccc1",
    # Ethers
    "ether":           "[OD2;!$(OC=O)]([#6])[#6]",
    # Amines
    "amine_primary":   "[NX3;H2;!$(NC=O)]",
    "amine_secondary": "[NX3;H1;!$(NC=O)]",
    "amine_tertiary":  "[NX3;H0;!$(NC=O);!$([N+])]",
    # Sulfur / nitrogen special
    "nitro":           "[NX3+](=O)[O-]",
    "sulfonamide":     "[SX4](=O)(=O)[NX3]",
    "urea":            "[NX3][CX3](=O)[NX3]",
    # Halogens
    "halogen_F":       "[F]",
    "halogen_Cl":      "[Cl]",
    "halogen_Br":      "[Br]",
}


class CleaningStage(BaseStage):
    """
    Normalise an input CSV into the unified pipeline format and write
    per-molecule metadata JSON files.

    Two modes
    ---------
    pipeline
        Output of filter.py.  Canonical column names, solubility in mol
        fraction, melting_temp in Kelvin.  No conversion performed.

    custom
        Any other CSV.  Only SMILES + a name column are required.
        Solubility and melting point are optional.

    Parameters
    ----------
    solubility_unit : str  (custom mode only, default: auto-detect)
        "mol_frac"   — values already in 0–1 range
        "mol_per_L"  — linear concentration, converted: x = c / (c + 55.51)
        "logS"       — log₁₀(mol/L), converted: c = 10^logS then mol frac
        If omitted and values are unambiguously in [0, 1] → mol_frac assumed.
        If omitted and values are ambiguous → column is OMITTED with a warning.

    overwrite_metadata : bool  (default: false)
        Whether to overwrite existing global metadata JSON files.
    """

    # ─────────────────────────────────────────────────────────────────────────
    # Execute
    # ─────────────────────────────────────────────────────────────────────────
    def execute(self):
        self.strict_mode = self.strict("cleaning")

        output_csv   = self.set_stage_output("cleaned.csv")
        metadata_dir = self.output_path("molecule_metadata")
        os.makedirs(metadata_dir, exist_ok=True)

        self.warnings = {
            "missing_melting_temp":    0,
            "implausible_mp":          0,
            "invalid_smiles":          0,
            "default_charge":          0,
            "default_multiplicity":    0,
            "duplicate_inchikey":      0,
            "impossible_solubility":   0,
            "zero_solubility":         0,
            "ambiguous_sol_unit":      0,
        }


        csv_list = self._resolve_input_sources()
        self._write_combined_raw_input(csv_list)

        lookup_ids = [os.path.splitext(os.path.basename(p))[0] for p in csv_list]
        self.set_items(lookup_ids)
        combined_rows = []

        for lookup_id in list(self.job.pending_items):
            self.log(
                f"[FILE]  ── {lookup_id}.csv "
                f"{'─' * max(1, 50 - len(lookup_id))}"
            )

            df = self._load_raw_csv(lookup_id, csv_list)
            if df is None:
                self.update_progress(lookup_id, success=False)
                continue

            mode = self._detect_mode(df, lookup_id)
            df   = self._normalise_columns(df, lookup_id, mode)

            if not self._validate_required_fields(df, lookup_id):
                self.update_progress(lookup_id, success=False)
                continue

            df = self._convert_solubility(df, lookup_id, mode)
            df = self._validate_solubility(df, lookup_id)

            if df.empty:
                self.log(
                    f"[ERROR] {lookup_id}: no rows remain after solubility "
                    "validation — skipping"
                )
                self.update_progress(lookup_id, success=False)
                continue

            df = self._canonicalise_smiles(df, lookup_id)
            df = self._generate_inchikeys(df, lookup_id)
            df = self._resolve_duplicates(df, lookup_id)

            source_file = self._find_source_file(lookup_id, csv_list)
            self._write_metadata(df, metadata_dir, source_file, lookup_id, mode)

            df = self._minimal_frame(df)
            combined_rows.append(df)

            self.log(
                f"[FILE]  {lookup_id} complete — "
                f"{len(df)} rows, {df['inchi_key'].nunique()} unique molecules"
            )
            self.update_progress(lookup_id)

        if not combined_rows:
            self.fail("No valid raw files processed — cannot create cleaned dataset.")

        final_df = pd.concat(combined_rows, ignore_index=True)
        final_df.to_csv(output_csv, index=False)

        self._log_warning_summary()
        self.log(
            f"[OUTPUT] cleaned.csv — "
            f"{len(final_df)} rows, "
            f"{final_df['inchi_key'].nunique()} unique InChIKeys"
        )

    # ─────────────────────────────────────────────────────────────────────────
    # Input resolution
    # ─────────────────────────────────────────────────────────────────────────
    def _resolve_input_sources(self):
        csv_list = []
        folders  = []

        if "input_folder" in self.parameters:
            folders = [self.parameters["input_folder"]]
        elif "input_folders" in self.parameters:
            raw     = self.parameters["input_folders"]
            folders = [raw] if isinstance(raw, str) else list(raw)

        for folder in folders:
            if not os.path.isdir(folder):
                self.fail(f"input_folder is not a directory: {folder}")
            for fname in sorted(os.listdir(folder)):
                if fname.lower().endswith(".csv"):
                    csv_list.append(os.path.join(folder, fname))

        if "input_csv" in self.parameters:
            raw = self.parameters["input_csv"]
            for p in ([raw] if isinstance(raw, str) else raw):
                self.require_file(p, "input_csv")
                csv_list.append(p)

        if not csv_list:
            self.fail("No CSV files found. Set input_folder(s) or input_csv.")

        csv_list = sorted(set(csv_list))
        self.log(f"[INPUT] {len(csv_list)} CSV file(s):")
        for p in csv_list:
            self.log(f"        {p}")
        return csv_list

    def _write_combined_raw_input(self, csv_list):
        frames = []
        for path in csv_list:
            try:
                df = pd.read_csv(path)
                df["__source_file"] = os.path.abspath(path)
                frames.append(df)
            except Exception as e:
                self.log(f"[WARNING] Could not read {path}: {e}")
        if not frames:
            self.fail("No valid CSVs could be read.")
        out = self.input_path("raw_combined.csv")
        pd.concat(frames, ignore_index=True).to_csv(out, index=False)
        self.log(f"[RAW]   Merged raw input → {out}")

    def _load_raw_csv(self, lookup_id, csv_list):
        matches = [p for p in csv_list
                   if os.path.splitext(os.path.basename(p))[0] == lookup_id]
        if not matches:
            self.log(f"[ERROR] No CSV found for '{lookup_id}'")
            return None
        try:
            df = pd.read_csv(matches[0])
            self.log(f"[LOAD]  {len(df)} rows, {len(df.columns)} columns")
            return df
        except Exception as e:
            self.log(f"[ERROR] Cannot read {matches[0]}: {e}")
            return None

    def _find_source_file(self, lookup_id, csv_list):
        for p in csv_list:
            if os.path.splitext(os.path.basename(p))[0] == lookup_id:
                return os.path.abspath(p)
        return None

    # ─────────────────────────────────────────────────────────────────────────
    # Mode detection
    # ─────────────────────────────────────────────────────────────────────────
    def _detect_mode(self, df: pd.DataFrame, lookup_id: str) -> str:
        """
        Pipeline signals (either is sufficient):
          • Column 'experimental_solubility_mol_frac' present
          • Column 'melting_temp' present AND median value > 100 K

        Everything else is custom.
        """
        cols = set(df.columns)

        if "experimental_solubility_mol_frac" in cols:
            self.log(
                f"[MODE]  pipeline — has 'experimental_solubility_mol_frac'"
            )
            return MODE_PIPELINE

        if "melting_temp" in cols:
            mt = pd.to_numeric(df["melting_temp"], errors="coerce").dropna()
            if len(mt) and mt.median() > 100:
                self.log(
                    f"[MODE]  pipeline — has 'melting_temp' "
                    f"(median {mt.median():.1f} K)"
                )
                return MODE_PIPELINE

        self.log(f"[MODE]  custom")
        return MODE_CUSTOM

    # ─────────────────────────────────────────────────────────────────────────
    # Column normalisation
    # ─────────────────────────────────────────────────────────────────────────
    @staticmethod
    def _norm(col: str) -> str:
        return col.strip().lower().replace(" ", "_")

    def _best_match(self, col: str) -> tuple[str | None, float]:
        """
        Return (target_key, score) for the best alias match, or (None, 0).
        Picks the single best scoring target across all alias lists.
        """
        col_norm    = self._norm(col)
        best_target = None
        best_score  = 0.0

        for target, aliases in _CUSTOM_ALIASES.items():
            hits = difflib.get_close_matches(col_norm, aliases, n=1, cutoff=0.60)
            if hits:
                from difflib import SequenceMatcher
                score = SequenceMatcher(None, col_norm, hits[0]).ratio()
                if score > best_score:
                    best_score  = score
                    best_target = target

        return (best_target, best_score) if best_score >= 0.60 else (None, 0.0)

    def _normalise_columns(
        self, df: pd.DataFrame, lookup_id: str, mode: str
    ) -> pd.DataFrame:
        """
        Pipeline mode  — minimal tidy-up; canonical names already present.
        Custom mode    — fuzzy-map each column to an internal target.

        Collision handling: if two raw columns map to the same target the
        higher-scoring one wins; the other is left under its original name
        (it will simply be ignored downstream).
        """
        df = df.copy()

        if mode == MODE_PIPELINE:
            df.columns = [c.strip() for c in df.columns]
            renames = {}
            for c in df.columns:
                n = self._norm(c)
                if n == "smiles" and c != "smiles":
                    renames[c] = "smiles"
                elif n in ("name", "mol_name") and "mol_name" not in df.columns:
                    renames[c] = "mol_name"
                elif n in ("inchikey", "inchi_key") and "inchi_key" not in df.columns:
                    renames[c] = "inchi_key"
            if renames:
                df = df.rename(columns=renames)
                for old, new in renames.items():
                    self.log(f"[HEADER] '{old}' → '{new}'")

            # If the file has mol_name_iupac but no mol_name, use it as the
            # display name.  mol_name_iupac is retained as-is so both columns
            # are available downstream.
            if "mol_name" not in df.columns and "mol_name_iupac" in df.columns:
                df["mol_name"] = df["mol_name_iupac"]
                self.log(
                    f"[HEADER] 'mol_name' absent — using 'mol_name_iupac' as display name"
                )

            return df

        # Custom mode
        TARGET_TO_CANON = {
            "smiles":              "smiles",
            "mol_name":            "mol_name",
            "sol_raw":             "sol_raw",           # intermediate name
            "melting_temp_kelvin": "melting_temp",
            "melting_temp_source": "melting_temp_source",
        }

        winner: dict[str, tuple[str, float]] = {}  # target -> (raw_col, score)

        for col in df.columns:
            target, score = self._best_match(col)
            if target is None:
                continue
            if target not in winner or score > winner[target][1]:
                if target in winner:
                    self.log(
                        f"[HEADER] '{col}' ({score:.2f}) replaces "
                        f"'{winner[target][0]}' ({winner[target][1]:.2f}) "
                        f"for target '{target}'"
                    )
                winner[target] = (col, score)

        rename_map = {}
        for target, (raw_col, _) in winner.items():
            canon = TARGET_TO_CANON.get(target)
            if canon and raw_col != canon:
                rename_map[raw_col] = canon
                self.log(f"[HEADER] '{raw_col}' → '{canon}'")

        if rename_map:
            df = df.rename(columns=rename_map)

        return df

    # ─────────────────────────────────────────────────────────────────────────
    # Required field validation
    # ─────────────────────────────────────────────────────────────────────────
    def _validate_required_fields(self, df: pd.DataFrame, lookup_id: str) -> bool:
        has_smiles = "smiles" in df.columns
        # mol_name_iupac is accepted as a name column — it is used as the
        # display name when mol_name is absent (see _normalise_columns)
        has_name   = "mol_name" in df.columns or "mol_name_iupac" in df.columns

        missing = []
        if not has_smiles:
            missing.append("SMILES")
        if not has_name:
            missing.append("Name / mol_name (or mol_name_iupac)")

        if missing:
            msg = (
                f"{lookup_id}: missing required field(s): "
                f"{', '.join(missing)}.  "
                "Minimum: a SMILES column and a name column "
                "(mol_name, Name, or mol_name_iupac)."
            )
            self.log(f"[ERROR] {msg}")
            if self.strict_mode:
                self.fail(msg)
            return False

        if "charge" not in df.columns:
            df["charge"] = None  # sentinel; resolved per molecule below (same as multiplicity)
            self.warnings["default_charge"] += 1
            self.log_warning(f"[HEADER] 'charge' column absent — will use RDKit or default 0 per molecule")

        if "multiplicity" not in df.columns:
            df["multiplicity"] = None  # sentinel; resolved per molecule below

        # ── Melting point status (logged prominently) ──────────────────────
        has_mp = "melting_temp" in df.columns
        if has_mp:
            mt      = pd.to_numeric(df["melting_temp"], errors="coerce")
            n_valid = mt.notna().sum()
            n_total = len(df)
            n_impl  = ((mt < MP_MIN_K) | (mt > MP_MAX_K)).sum()
            self.log(
                f"[MP]    Melting point column present — "
                f"{n_valid}/{n_total} values populated"
                + (f", {n_impl} outside plausibility window [{MP_MIN_K}–{MP_MAX_K} K]"
                   if n_impl else "")
            )
        else:
            self.log_warning(
                f"[MP]    No melting point column detected in "
                f"{lookup_id}. PubChem will be queried for every molecule. "
                "To avoid this, include a 'Melting point' column (Kelvin)."
            )

        # ── Solubility status ──────────────────────────────────────────────
        has_sol = (
            "experimental_solubility_mol_frac" in df.columns
            or "sol_raw" in df.columns
        )
        self.log(
            f"[VALIDATE] SMILES ✓  Name ✓  "
            f"Solubility {'✓' if has_sol else '— absent (will be omitted)'}  "
            f"Melting pt {'✓' if has_mp else '— absent (PubChem fallback)'}"
        )
        return True

    # ─────────────────────────────────────────────────────────────────────────
    # Solubility conversion
    # ─────────────────────────────────────────────────────────────────────────
    def _convert_solubility(
        self, df: pd.DataFrame, lookup_id: str, mode: str
    ) -> pd.DataFrame:
        """
        Ensure 'experimental_solubility_mol_frac' is in mol fraction.

        Pipeline mode  — already correct, nothing to do.
        Custom mode    — convert 'sol_raw' according to solubility_unit parameter.

        Auto-detection (no parameter)
        ─────────────────────────────
        All values in [0, 1] → assume mol_frac, log the assumption.
        Any value > 1 or < 0 → cannot safely assume a unit; log a clear
        warning and OMIT the column rather than store corrupt data.
        """
        target = "experimental_solubility_mol_frac"
        df     = df.copy()

        if mode == MODE_PIPELINE:
            if target not in df.columns:
                self.log(
                    f"[SOL]   pipeline mode — '{target}' not found, "
                    "solubility absent in output"
                )
            else:
                self.log(f"[SOL]   pipeline mode — mol frac already present")
            return df

        # Custom mode
        if "sol_raw" not in df.columns:
            self.log(f"[SOL]   No solubility column found — omitting from output")
            return df

        sol_raw = pd.to_numeric(df["sol_raw"], errors="coerce")
        unit    = self.parameters.get("solubility_unit", "").strip().lower()

        if not unit:
            has_neg  = sol_raw.lt(0).any()
            has_over = sol_raw.gt(1.0).any()

            if not has_neg and not has_over:
                unit = UNIT_MOL_FRAC
                self.log(
                    f"[SOL]   solubility_unit not set — "
                    "values all in [0, 1], assuming mol_frac"
                )
            else:
                clues = []
                if has_neg:
                    clues.append(
                        f"{sol_raw.lt(0).sum()} negative value(s) — suggests logS"
                    )
                if has_over:
                    clues.append(
                        f"{sol_raw.gt(1).sum()} value(s) > 1.0 — suggests mol/L"
                    )
                self.log(
                    f"[SOL]   WARNING: solubility unit ambiguous "
                    f"({'; '.join(clues)}). "
                    f"Set parameter 'solubility_unit' to one of: "
                    f"'{UNIT_MOL_FRAC}', '{UNIT_MOL_PER_L}', '{UNIT_LOGS}'. "
                    "Solubility column OMITTED from output until resolved."
                )
                self.warnings["ambiguous_sol_unit"] += 1
                return df.drop(columns=["sol_raw"])

        if unit == UNIT_MOL_FRAC:
            df[target] = sol_raw
            self.log(f"[SOL]   mol_frac — used directly")

        elif unit == UNIT_MOL_PER_L:
            df[target] = sol_raw / (sol_raw + WATER_MOLARITY)
            self.log(f"[SOL]   mol/L → mol frac")

        elif unit == UNIT_LOGS:
            c = 10.0 ** sol_raw
            df[target] = c / (c + WATER_MOLARITY)
            self.log(
                f"[SOL]   logS → mol frac  "
                f"(logS range [{sol_raw.min():.2f}, {sol_raw.max():.2f}])"
            )

        else:
            self.fail(
                f"Unknown solubility_unit '{unit}'. "
                f"Valid: '{UNIT_MOL_FRAC}', '{UNIT_MOL_PER_L}', '{UNIT_LOGS}'"
            )

        return df.drop(columns=["sol_raw"])

    # ─────────────────────────────────────────────────────────────────────────
    # Solubility validation  (values always in mol frac at this point)
    # ─────────────────────────────────────────────────────────────────────────
    def _validate_solubility(self, df: pd.DataFrame, lookup_id: str) -> pd.DataFrame:
        col = "experimental_solubility_mol_frac"
        if col not in df.columns:
            return df

        df[col] = pd.to_numeric(df[col], errors="coerce")
        before  = len(df)
        s       = df[col]

        null_mask = s.isna()
        zero_mask = s == 0
        neg_mask  = s.lt(0) & ~null_mask
        over_mask = s.gt(1.0)
        bad       = null_mask | zero_mask | neg_mask | over_mask

        if bad.any():
            for mask, label, wkey in [
                (null_mask, "null",                   None),
                (zero_mask, "== 0",                   "zero_solubility"),
                (neg_mask,  "< 0",                    None),
                (over_mask, "> 1.0 — check solubility_unit", "impossible_solubility"),
            ]:
                if not mask.any():
                    continue
                names = df.loc[mask, "mol_name"].tolist() if "mol_name" in df.columns else []
                self.log(
                    f"[SOL]   Excluding {mask.sum()} row(s) solubility {label}: "
                    f"{names}"
                )
                if wkey:
                    self.warnings[wkey] += mask.sum()
            df = df[~bad].reset_index(drop=True)
            self.log(f"[SOL]   {before} → {len(df)} rows after validation")
        else:
            self.log(f"[SOL]   All {len(df)} values valid ✓")
        return df

    # ─────────────────────────────────────────────────────────────────────────
    # SMILES canonicalisation
    # ─────────────────────────────────────────────────────────────────────────
    def _canonicalise_smiles(self, df: pd.DataFrame, lookup_id: str) -> pd.DataFrame:
        def _canon(smi):
            try:
                mol = Chem.MolFromSmiles(str(smi))
                return Chem.MolToSmiles(mol, canonical=True) if mol else smi
            except Exception:
                return smi

        orig         = df["smiles"].astype(str).copy()
        df["smiles"] = orig.apply(_canon)
        changed      = (df["smiles"] != orig).sum()
        if changed:
            self.log(f"[SMILES] Canonicalised {changed} SMILES string(s)")

        if "mol_name" in df.columns:
            df["mol_name"] = df["mol_name"].astype(str).str.strip()

        return df

    # ─────────────────────────────────────────────────────────────────────────
    # InChIKey generation
    # ─────────────────────────────────────────────────────────────────────────
    def _smiles_to_inchikey(self, smi: str):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return None
            return inchi.InchiToInchiKey(inchi.MolToInchi(mol))
        except Exception:
            return None

    def _generate_inchikeys(self, df: pd.DataFrame, lookup_id: str) -> pd.DataFrame:
        # Always regenerate from SMILES for consistency (canonical SMILES
        # guarantees the InChIKey is correct even if one was pre-supplied)
        df["inchi_key"] = df["smiles"].apply(self._smiles_to_inchikey)
        invalid         = df["inchi_key"].isna()
        if invalid.any():
            self.log(
                f"[SMILES] {invalid.sum()} SMILES failed InChIKey generation "
                f"(retained, flagged): {df.loc[invalid, 'smiles'].tolist()}"
            )
            self.warnings["invalid_smiles"] += invalid.sum()
        else:
            self.log(f"[SMILES] All {len(df)} InChIKeys generated ✓")
        return df

    # ─────────────────────────────────────────────────────────────────────────
    # Duplicate InChIKey resolution
    # ─────────────────────────────────────────────────────────────────────────
    def _resolve_duplicates(self, df: pd.DataFrame, lookup_id: str) -> pd.DataFrame:
        dupe_mask = df["inchi_key"].notna() & df["inchi_key"].duplicated(keep=False)
        if not dupe_mask.any():
            self.log(f"[DEDUP]  No duplicate InChIKeys ✓")
            return df

        sol      = "experimental_solubility_mol_frac"
        before   = len(df)
        resolved = []

        for key, group in df.groupby("inchi_key", sort=False):
            if len(group) == 1:
                resolved.append(group)
                continue

            keeper = group.iloc[[0]].copy()

            if sol in group.columns:
                vals = group[sol].dropna()
                if len(vals) == 0:
                    pass
                elif vals.nunique() == 1:
                    keeper[sol] = vals.iloc[0]
                else:
                    mean_val = vals.mean()
                    pct      = (vals.max() - vals.min()) / vals.mean() * 100
                    self.warnings["duplicate_inchikey"] += len(group) - 1
                    self.log(
                        f"[DEDUP]  {key}: {len(group)} rows, "
                        f"spread {pct:.1f}% of mean → mean {mean_val:.3e}. "
                        f"Names: {group.get('mol_name', pd.Series(['?'])).tolist()}"
                    )
                    keeper[sol] = mean_val

            resolved.append(keeper)

        df_out = pd.concat(resolved, ignore_index=True)
        self.log(
            f"[DEDUP]  {before} → {len(df_out)} "
            f"({before - len(df_out)} duplicate(s) resolved)"
        )
        return df_out

    # ─────────────────────────────────────────────────────────────────────────
    # Melting temperature resolution  (per row)
    # ─────────────────────────────────────────────────────────────────────────
    def _resolve_mp(self, row) -> tuple:
        """
        Return (kelvin: float|None, source: str).

        Tries the 'melting_temp' column first; if absent or implausible,
        returns (None, "missing") and the caller will use PubChem.
        """
        def _scalar(val):
            if isinstance(val, pd.Series):
                v = val.dropna()
                return float(v.iloc[0]) if len(v) else None
            try:
                f = float(val)
                return None if pd.isna(f) else f
            except (TypeError, ValueError):
                return None

        kt = _scalar(row.get("melting_temp"))
        if kt is None:
            return None, "missing"

        if MP_MIN_K <= kt <= MP_MAX_K:
            raw_src = row.get("melting_temp_source")
            if isinstance(raw_src, pd.Series):
                raw_src = raw_src.dropna()
                raw_src = raw_src.iloc[0] if len(raw_src) else None
            if raw_src is not None and str(raw_src).strip() not in ("", "nan", "None", "N/A"):
                src = str(raw_src).strip()
            else:
                src = "user_input"
            return kt, src

        name = row.get("mol_name", "?")
        self.log(
            f"[MP]    Implausible value for '{name}': {kt:.1f} K "
            f"({kt - 273.15:.1f} °C) — outside [{MP_MIN_K}, {MP_MAX_K}] K. "
            "PubChem will be tried."
        )
        self.warnings["implausible_mp"] += 1
        return None, "missing"

    # ─────────────────────────────────────────────────────────────────────────
    # Metadata writer
    # ─────────────────────────────────────────────────────────────────────────
    def _write_metadata(
        self,
        df:           pd.DataFrame,
        metadata_dir: str,
        source_file,
        lookup_id:    str,
        mode:         str,
    ):
        from rdkit.Chem import Lipinski as L, rdMolDescriptors as R, Crippen as C
        from modules.utils.melting_point_finder import MeltingPointFinder
        _mp_finder = MeltingPointFinder()

        pipeline_version = self.config.get("pipeline_version", "unknown")
        timestamp        = datetime.utcnow().isoformat() + "Z"
        overwrite        = self.parameters.get("overwrite_metadata", False)
        mp_counts        = {"provided": 0, "pubchem": 0, "missing": 0}
        charge_counts    = {"user_input": 0, "rdkit_computed": 0, "default": 0}
        mult_counts      = {"user_input": 0, "rdkit_computed": 0, "default": 0}

        def _safe_float(val):
            if isinstance(val, pd.Series):
                val = val.dropna()
                val = val.iloc[0] if len(val) else None
            try:
                f = float(val)
                return None if pd.isna(f) else f
            except (TypeError, ValueError):
                return None

        for ik, group in df.groupby("inchi_key"):
            if not ik:
                continue

            row    = group.iloc[0]
            smiles = str(row.get("smiles", ""))
            mol    = Chem.MolFromSmiles(smiles)

            # ── Charge (3-case) ────────────────────────────────────────────
            raw_charge = row.get("charge")
            if raw_charge is not None and str(raw_charge).strip() not in ("", "nan"):
                charge        = int(raw_charge)
                charge_source = "user_input"
            elif mol is not None:
                charge        = Chem.GetFormalCharge(mol)
                charge_source = "rdkit_computed"
            else:
                charge        = 0
                charge_source = "default"
                self.log_warning(
                    f"[META]  {ik}: charge not provided and RDKit cannot parse "
                    f"SMILES — defaulted to 0"
                )
            charge_counts[charge_source] += 1

            # ── Multiplicity (3-case) ──────────────────────────────────────
            raw_mult = row.get("multiplicity")
            if raw_mult is not None and str(raw_mult).strip() not in ("", "nan", "None"):
                multiplicity        = int(raw_mult)
                multiplicity_source = "user_input"
            elif mol is not None and Descriptors.NumRadicalElectrons(mol) == 0:
                multiplicity        = 1
                multiplicity_source = "rdkit_computed"
            else:
                multiplicity        = 1
                multiplicity_source = "default"
                self.log_warning(
                    f"[META]  {ik}: multiplicity not provided and cannot be "
                    f"reliably determined — defaulted to 1"
                )
            mult_counts[multiplicity_source] += 1

            # Physchem descriptors
            if mol is not None:
                phys = {
                    "rotatable_bonds":    L.NumRotatableBonds(mol),
                    "molecular_weight":   R.CalcExactMolWt(mol),
                    "heavy_atom_count":   mol.GetNumHeavyAtoms(),
                    "hbond_donors":       L.NumHDonors(mol),
                    "hbond_acceptors":    L.NumHAcceptors(mol),
                    "tpsa":               R.CalcTPSA(mol),
                    "logp":               C.MolLogP(mol),
                    "aromatic_rings":     R.CalcNumAromaticRings(mol),
                    "Fsp3":               R.CalcFractionCSP3(mol),
                    "Bertz_complexity":   Descriptors.BertzCT(mol),
                    "molar_refractivity": C.MolMR(mol),
                }
                # Functional group counts (SMARTS-based)
                fg_counts = {}
                for fg_name, smarts in FUNCTIONAL_GROUP_SMARTS.items():
                    pat = Chem.MolFromSmarts(smarts)
                    fg_counts[fg_name] = len(mol.GetSubstructMatches(pat)) if pat is not None else 0
                phys["functional_groups"] = fg_counts
            else:
                phys = {k: None for k in (
                    "rotatable_bonds", "molecular_weight", "heavy_atom_count",
                    "hbond_donors", "hbond_acceptors", "tpsa", "logp",
                    "aromatic_rings", "Fsp3", "Bertz_complexity",
                    "molar_refractivity",
                )}
                phys["functional_groups"] = {fg: None for fg in FUNCTIONAL_GROUP_SMARTS}
                self.log(
                    f"[META]  {ik}: RDKit cannot parse SMILES — descriptors null"
                )

            meta = {
                "metadata_version": METADATA_VERSION,
                "inchi_key":        ik,
                "smiles":           smiles,
                "mol_name":         str(row.get("mol_name", DEFAULT_METADATA["mol_name"])),
                "mol_name_iupac":   str(row.get("mol_name_iupac",
                                                  DEFAULT_METADATA["mol_name_iupac"])),
                "charge":               charge,
                "charge_source":        charge_source,
                "multiplicity":         multiplicity,
                "multiplicity_source":  multiplicity_source,
                **phys,
                "provenance": {
                    "source_file":          source_file,
                    "input_mode":           mode,
                    "cleaning_timestamp":   timestamp,
                    "pipeline_version":     pipeline_version,
                    "generated_by_request": self.job.request_id,
                    "generated_by_job":     self.job.job_id,
                },
            }

            # Solubility (experimental + AqSol predicted if present)
            for field in (
                "experimental_solubility_mol_frac",
                "aqsol_predicted_solubility_mol_frac",
            ):
                v = _safe_float(row.get(field))
                if v is not None:
                    meta[field] = v

            # Fill any remaining defaults
            for k, default in DEFAULT_METADATA.items():
                if k not in meta:
                    meta[k] = default

            # ── Melting temperature ────────────────────────────────────────
            kt, src = self._resolve_mp(row)

            if kt is not None:
                meta["melting_temp"]        = kt
                meta["melting_temp_c"]      = round(kt - 273.15, 4)
                meta["melting_temp_source"] = src
                mp_counts["provided"] += 1
            else:
                # PubChem + NIST fallback via MeltingPointFinder
                mp_result = _mp_finder.get_best(ik)
                mt_raw    = mp_result.melting_temp
                mt_src    = mp_result.source

                if mt_raw is not None:
                    try:
                        kt_pub = float(mt_raw)
                        if MP_MIN_K <= kt_pub <= MP_MAX_K:
                            meta["melting_temp"]        = kt_pub
                            meta["melting_temp_c"]      = round(kt_pub - 273.15, 4)
                            meta["melting_temp_source"] = mt_src
                            mp_counts["pubchem"] += 1
                        else:
                            self.log_warning(
                                f"[MP]    {ik}: implausible PubChem value "
                                f"{kt_pub:.1f} K — marked N/A, no melting point available"
                            )
                            self.warnings["implausible_mp"] += 1
                            meta["melting_temp"]        = "N/A"
                            meta["melting_temp_source"] = "implausible_pubchem"
                            mp_counts["missing"] += 1
                    except (TypeError, ValueError):
                        meta["melting_temp"] = "N/A"
                        mp_counts["missing"] += 1
                else:
                    self.log_warning(
                        f"[MP]    {ik}: no melting point in CSV or PubChem — "
                        "COSMO-RS accuracy may be reduced"
                    )
                    meta["melting_temp"]        = "N/A"
                    meta["melting_temp_source"] = "N/A"
                    mp_counts["missing"] += 1
                    self.warnings["missing_melting_temp"] += 1

            # ── Write JSON files ───────────────────────────────────────────
            local_path = os.path.join(metadata_dir, f"{ik}.json")
            with open(local_path, "w") as f:
                f.write(json.dumps(meta, indent=2))

            global_dir = self.config.get("constant_files", {}).get("metadata_dir")
            if not global_dir:
                continue
            os.makedirs(global_dir, exist_ok=True)
            global_path = os.path.join(global_dir, f"{ik}.json")
            if overwrite or not os.path.exists(global_path):
                with open(global_path, "w") as f:
                    f.write(json.dumps(meta, indent=2))

        self.log(
            f"[MP]    {lookup_id}: "
            f"provided={mp_counts['provided']}, "
            f"pubchem={mp_counts['pubchem']}, "
            f"missing={mp_counts['missing']}"
        )

        # ── Charge/multiplicity source breakdown ───────────────────────────
        n_mols = sum(charge_counts.values())
        if n_mols:
            non_user_charge = charge_counts["rdkit_computed"] + charge_counts["default"]
            if non_user_charge > 0:
                self.log_warning(
                    f"[META]  {lookup_id}: charge not fully user-provided "
                    f"({n_mols} molecule(s)) — "
                    f"user_input={charge_counts['user_input']}, "
                    f"rdkit_computed={charge_counts['rdkit_computed']}, "
                    f"default={charge_counts['default']}"
                )
                self.warnings["default_charge"] += charge_counts["default"]

            non_user_mult = mult_counts["rdkit_computed"] + mult_counts["default"]
            if non_user_mult > 0:
                self.log_warning(
                    f"[META]  {lookup_id}: multiplicity not fully user-provided "
                    f"({n_mols} molecule(s)) — "
                    f"user_input={mult_counts['user_input']}, "
                    f"rdkit_computed={mult_counts['rdkit_computed']} (closed-shell), "
                    f"default={mult_counts['default']} (radical/open-shell assumed singlet)"
                )
                self.warnings["default_multiplicity"] += mult_counts["default"]

    # ─────────────────────────────────────────────────────────────────────────
    # Minimal cleaned CSV  (downstream stages read this)
    # ─────────────────────────────────────────────────────────────────────────
    def _minimal_frame(self, df: pd.DataFrame) -> pd.DataFrame:
        keep = [
            "inchi_key",
            "smiles",
            "charge",
            "multiplicity",
            "mol_name",
            "mol_name_iupac",
            "experimental_solubility_mol_frac",
            "aqsol_predicted_solubility_mol_frac",
        ]
        return df[[c for c in keep if c in df.columns]]

    # ─────────────────────────────────────────────────────────────────────────
    # Warning summary
    # ─────────────────────────────────────────────────────────────────────────
    def _log_warning_summary(self):
        w     = self.warnings
        total = sum(w.values())
        if total == 0:
            self.log("[SUMMARY] Cleaning completed with no warnings.")
            return

        self.log("[SUMMARY] Cleaning completed with warnings:")
        rows = [
            ("missing_melting_temp",
             "molecule(s) with no melting point — may reduce COSMO-RS accuracy"),
            ("implausible_mp",
             f"MP value(s) outside [{MP_MIN_K}–{MP_MAX_K} K] — marked N/A"),
            ("invalid_smiles",
             "SMILES that failed InChIKey generation — retained, flagged"),
            ("default_charge",
             "molecule(s) where charge could not be resolved (SMILES invalid) — defaulted to 0"),
            ("default_multiplicity",
             "molecule(s) where multiplicity could not be determined (possible radical) — defaulted to 1"),
            ("duplicate_inchikey",
             "row(s) collapsed — solubility averaged within duplicate groups"),
            ("impossible_solubility",
             "row(s) excluded — solubility > 1.0; check solubility_unit parameter"),
            ("zero_solubility",
             "row(s) excluded — solubility == 0"),
            ("ambiguous_sol_unit",
             "file(s) where solubility unit was ambiguous — column omitted; "
             "set 'solubility_unit' parameter"),
        ]
        for key, desc in rows:
            if w[key]:
                self.log(f"          {key:<28}: {w[key]:>4}  {desc}")