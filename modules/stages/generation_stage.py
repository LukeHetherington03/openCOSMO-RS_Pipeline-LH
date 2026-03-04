#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
modules/stages/generation_stage.py

Conformer generation stage.  Strictly generation only — no optimisation,
no pruning, no forcefield minimisation beyond what is intrinsic to the
generator (e.g. RDKit MMFF energy evaluation for sorting purposes).

=========================================================================
INPUT
=========================================================================

  stage_input : cleaned.csv from CleaningStage
    Required columns: inchi_key, smiles
    Optional: charge, mol_name

=========================================================================
OUTPUT  (canonical stage output)
=========================================================================

  energies.json   — ConformerSet JSON
    One entry per generated conformer.

=========================================================================
AUXILIARY OUTPUTS
=========================================================================

  summary.csv
  generation_metadata.json
  xyz/<inchi_key>_conf<NNN>.xyz

=========================================================================
BACKENDS
=========================================================================

  rdkit      — ETKDGv3 embedding + MMFF94/UFF energy evaluation
               (energy is the FF energy at the generated geometry, not an
                optimisation — geometry is not updated)
  crest      — Single RDKit seed → CREST conformer search
               (note: CREST is not validated on the current HPC env;
                tested separately before production use)
  openbabel  — RDKit seed → OBConformerSearch GA
               (note: OpenBabel backend is not validated on the current
                HPC env; tested separately before production use)

=========================================================================
ITEM TRACKING
=========================================================================

  Items are inchi_keys (one per unique molecule in the CSV).
  Multiple rows with the same inchi_key are collapsed to one item —
  the first row's SMILES is used.

  update_progress(inchi_key) is called after each molecule completes.
  This stage is serial — no run_parallel().

=========================================================================
STRICT MODE
=========================================================================

  self.strict("generation") reads from config["generation"]["strict"].
  When True, any per-molecule failure raises and aborts the stage.

"""

import json
import os
import shutil
import subprocess
from datetime import datetime

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

from modules.stages.base_stage import BaseStage
from modules.utils.atomic_write import AtomicWriter
from modules.utils.conformers import ConformerRecord, ConformerSet

try:
    from openbabel import openbabel as ob
    HAVE_OPENBABEL = True
except ImportError:
    HAVE_OPENBABEL = False


GENERATION_BACKENDS = {
    "rdkit":     "_backend_rdkit",
    "crest":     "_backend_crest",
    "openbabel": "_backend_openbabel",
}


class GenerationStage(BaseStage):
    """
    Conformer generation stage.

    Serial execution — item granularity is one molecule (inchi_key).
    Backends: rdkit (default), crest, openbabel.

    execute() sets items = unique inchi_keys from the CSV, then iterates
    calling update_progress() after each molecule.

    Stages never call mark_complete() — BaseStage.run() owns lifecycle.
    """

    # =========================================================================
    # Entry point
    # =========================================================================

    def execute(self):
        self.strict_mode = self.strict("generation")
        self.set_stage_output("energies.json")

        stage_input = self.get_stage_input()
        self.require_file(stage_input, "cleaned dataset")

        # Copy cleaned.csv into inputs/ for reproducibility
        inputs_cleaned = os.path.join(self.inputs_dir, "cleaned.csv")
        os.makedirs(self.inputs_dir, exist_ok=True)
        if os.path.abspath(stage_input) != os.path.abspath(inputs_cleaned):
            shutil.copy(stage_input, inputs_cleaned)
        self.log_info(f"Stage input: {stage_input}")

        params    = self.parameters
        engine    = params.get("engine", "rdkit").lower()
        num_confs = params.get("n")
        seed      = params.get("seed", 42)
        threads   = params.get("threads", 1)

        self.log_config(f"engine={engine}  seed={seed}  threads={threads}")
        self.log_config(f"strict={self.strict_mode}")

        # Prepare directories
        xyz_dir = os.path.join(self.outputs_dir, "xyz")
        os.makedirs(xyz_dir, exist_ok=True)

        # Load cleaned dataset
        df_all = self._load_csv(stage_input)

        # Load metadata (charge, spin per molecule)
        metadata_index = self._load_metadata_index()

        # Validate rows and deduplicate to one item per inchi_key
        valid_rows = self._validate_rows(df_all)

        if not valid_rows:
            self.fail("No valid molecules found in cleaned dataset.")

        # Items are inchi_keys (strings) — one per unique molecule
        inchi_keys = list(dict.fromkeys(ik for _, ik, _, _ in valid_rows))
        self.set_items(inchi_keys)

        # Conformer count
        if num_confs is None:
            num_confs = self._default_num_confs(valid_rows)
            self.log_info(f"Conformer count heuristic: n={num_confs}")
        else:
            self.log_info(f"Conformer count (explicit): n={num_confs}")

        self.warnings = {
            "embedding_failures":    0,
            "energy_eval_failures":  0,
            "zero_conformers":       0,
            "crest_failures":        0,
            "openbabel_failures":    0,
        }

        # Version strings for provenance
        rdkit_version     = getattr(Chem, "__version__", None)
        crest_version     = self._get_crest_version()
        openbabel_version = self._get_openbabel_version() if HAVE_OPENBABEL else None

        # Dispatch backend
        if engine not in GENERATION_BACKENDS:
            self.fail(f"Unknown generation backend: '{engine}'. "
                      f"Valid: {list(GENERATION_BACKENDS)}")

        runner = getattr(self, GENERATION_BACKENDS[engine])
        self.log_section(f"Backend: {engine}")

        conformer_set = runner(
            valid_rows        = valid_rows,
            num_confs         = num_confs,
            seed              = seed,
            threads           = threads,
            xyz_out_dir       = xyz_dir,
            metadata_index    = metadata_index,
            rdkit_version     = rdkit_version,
            crest_version     = crest_version,
            openbabel_version = openbabel_version,
        )

        if len(conformer_set) == 0:
            self.fail("Generation produced zero conformers.")

        self._write_outputs(
            conformer_set = conformer_set,
            xyz_dir       = xyz_dir,
            num_valid     = len(inchi_keys),
            total_rows    = len(df_all),
            backend       = engine,
            num_confs     = num_confs,
            seed          = seed,
        )

        self._log_warning_summary()

    # =========================================================================
    # Input loading
    # =========================================================================

    def _load_csv(self, path: str) -> pd.DataFrame:
        try:
            df = pd.read_csv(path)
        except Exception as e:
            self.fail(f"Failed to load stage_input {path}: {e}")
        self.log_info(f"Loaded {len(df)} rows")
        return df

    # =========================================================================
    # Metadata index
    # =========================================================================

    def _load_metadata_index(self) -> dict:
        """
        Load per-molecule metadata from constant_files.metadata_dir.
        Returns dict of {inchi_key: metadata_dict}.
        Falls back to empty dict (charge=0, spin=1 defaults) if absent.
        """
        metadata_index = {}
        const = self.config.get("constant_files", {})
        meta_dir = const.get("metadata_dir")

        if not meta_dir or not os.path.isdir(meta_dir):
            self.log_info(
                "No metadata_dir found — using defaults: charge=0, spin=1"
            )
            return metadata_index

        for fname in os.listdir(meta_dir):
            if not fname.endswith(".json"):
                continue
            ik = fname[:-5]
            try:
                with open(os.path.join(meta_dir, fname)) as f:
                    metadata_index[ik] = json.load(f)
            except Exception as e:
                self.log_warning(f"Could not read metadata for {ik}: {e}")

        self.log_info(f"Loaded metadata for {len(metadata_index)} molecules")
        return metadata_index

    def _get_charge_spin(
        self, inchi_key: str, metadata_index: dict
    ) -> tuple[int, int]:
        """
        Resolve charge and spin.
        Hierarchy: parameters > metadata JSON > defaults (0, 1).
        """
        params = self.parameters
        meta   = metadata_index.get(inchi_key, {})

        charge = (
            params.get("charge")
            if params.get("charge") is not None
            else meta.get("charge", 0)
        )
        spin = (
            params.get("spin")
            if params.get("spin") is not None
            else meta.get("spin", 1)
        )
        return int(charge), int(spin)

    # =========================================================================
    # Row validation
    # =========================================================================

    def _validate_rows(self, df: pd.DataFrame) -> list:
        """
        Validate each row has a valid inchi_key and parseable SMILES.
        Returns list of (row_idx, inchi_key, smiles_canon, row).

        Deduplication: if multiple rows share an inchi_key, only the first
        valid row for that key is kept.  Generation is per-molecule.
        """
        seen   = {}   # inchi_key -> first valid row_idx
        valid  = []
        params = self.parameters

        for idx, row in df.iterrows():
            ik     = row.get("inchi_key")
            smiles = row.get("smiles")

            if not ik or not isinstance(ik, str) or not ik.strip():
                self.log_warning(f"Row {idx}: missing inchi_key — skipped")
                if self.strict_mode:
                    self.fail(f"Row {idx}: missing inchi_key")
                continue

            if not smiles:
                self.log_warning(f"Row {idx} ({ik}): missing SMILES — skipped")
                if self.strict_mode:
                    self.fail(f"Row {idx}: missing SMILES")
                continue

            mol = Chem.MolFromSmiles(str(smiles))
            if mol is None:
                self.log_warning(
                    f"Row {idx} ({ik}): invalid SMILES '{smiles}' — skipped"
                )
                if self.strict_mode:
                    self.fail(f"Row {idx}: invalid SMILES")
                continue

            if ik in seen:
                continue   # deduplicate — first occurrence wins

            smiles_canon = Chem.MolToSmiles(mol, canonical=True)
            seen[ik] = idx
            valid.append((idx, ik, smiles_canon, row))

        self.log_info(f"Valid unique molecules: {len(valid)}")
        return valid

    # =========================================================================
    # Conformer count heuristic
    # =========================================================================

    def _default_num_confs(self, valid_rows: list) -> int:
        """
        O'Boyle 2011 heuristic by rotatable bond count.
        Returns the max across all molecules (safe upper bound).
        """
        max_n = 0
        for _, _, smiles, _ in valid_rows:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            try:
                mol = Chem.RemoveHs(mol)
            except Exception:
                pass
            try:
                n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
            except Exception:
                n_rot = 8
            n = 50 if n_rot <= 7 else (200 if n_rot <= 12 else 300)
            max_n = max(max_n, n)
        return max_n if max_n > 0 else 50

    # =========================================================================
    # Version detection
    # =========================================================================

    def _get_crest_version(self) -> str | None:
        exe = self.config.get("crest", {}).get("executable")
        if not exe or not os.path.isfile(exe):
            return None
        try:
            r = subprocess.run(
                [exe, "--version"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            )
            if r.returncode == 0:
                return r.stdout.strip().splitlines()[0].strip()
        except Exception:
            pass
        return None

    def _get_openbabel_version(self) -> str | None:
        if not HAVE_OPENBABEL:
            return None
        try:
            return ob.OBReleaseVersion()
        except Exception:
            return None

    # =========================================================================
    # Backend: RDKit
    # =========================================================================

    def _backend_rdkit(
        self,
        valid_rows,
        num_confs,
        seed,
        threads,
        xyz_out_dir,
        metadata_index,
        rdkit_version,
        **kwargs,
    ) -> ConformerSet:
        """
        ETKDGv3 embedding.  MMFF94 (or UFF fallback) energy evaluated at
        the generated geometry — no geometry update is performed here.
        Sorting/pruning by energy is for PruningStage.
        """
        conformer_set = ConformerSet()
        timestamp     = datetime.utcnow().isoformat() + "Z"

        for idx, inchi_key, smiles, row in valid_rows:
            self.log_section(f"RDKit  {inchi_key}")

            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))

            emb = AllChem.ETKDGv3()
            emb.randomSeed            = seed
            emb.numThreads            = threads
            emb.useSmallRingTorsions  = True
            emb.useMacrocycleTorsions = True
            emb.enforceChirality      = True

            try:
                conf_ids = AllChem.EmbedMultipleConfs(mol, num_confs, emb)
            except Exception as e:
                self.log_warning(f"{inchi_key}: embedding failed — {e}")
                if self.strict_mode:
                    self.fail(f"Embedding failed for {inchi_key}: {e}")
                self.warnings["embedding_failures"] += 1
                self.update_progress(inchi_key, success=False)
                continue

            if len(conf_ids) == 0:
                self.log_warning(f"{inchi_key}: zero conformers embedded")
                if self.strict_mode:
                    self.fail(f"No conformers embedded for {inchi_key}")
                self.warnings["zero_conformers"] += 1
                self.update_progress(inchi_key, success=False)
                continue

            # MMFF94 energy evaluation (no geometry update)
            mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
            use_mmff   = mmff_props is not None
            forcefield = "MMFF94" if use_mmff else "UFF"

            charge, spin = self._get_charge_spin(inchi_key, metadata_index)

            for local_idx, conf_id in enumerate(conf_ids):
                try:
                    if use_mmff:
                        ff = AllChem.MMFFGetMoleculeForceField(
                            mol, mmff_props, confId=conf_id
                        )
                    else:
                        ff = AllChem.UFFGetMoleculeForceField(
                            mol, confId=conf_id
                        )
                    energy = float(ff.CalcEnergy())
                except Exception as e:
                    self.log_warning(
                        f"{inchi_key} conf{local_idx:03d}: "
                        f"energy evaluation failed — {e}"
                    )
                    if self.strict_mode:
                        self.fail(str(e))
                    self.warnings["energy_eval_failures"] += 1
                    continue

                xyz_path = os.path.join(
                    xyz_out_dir, f"{inchi_key}_conf{local_idx:03d}.xyz"
                )
                self._write_xyz(mol, conf_id, xyz_path)

                record = ConformerRecord(
                    inchi_key = inchi_key,
                    conf_num  = local_idx,
                    xyz_path  = xyz_path,
                    energy    = energy,
                    smiles    = smiles,
                    provenance = {
                        "backend":              "rdkit",
                        "rdkit_version":        rdkit_version,
                        "seed":                 seed,
                        "threads":              threads,
                        "forcefield":           forcefield,
                        "generation_timestamp": timestamp,
                        "source_row":           int(idx),
                        "charge":               charge,
                        "spin":                 spin,
                    },
                )
                record.optimisation_history.append({
                    "stage":     "generation",
                    "engine":    "rdkit",
                    "energy":    energy,
                    "xyz_path":  xyz_path,
                    "timestamp": timestamp,
                })
                conformer_set.add(record)

            self.log_info(
                f"{inchi_key}: {len(conf_ids)} conformers embedded, "
                f"{len([r for r in conformer_set.records if r.inchi_key == inchi_key])} "
                f"with valid energy"
            )
            self.update_progress(inchi_key)

        return conformer_set

    # =========================================================================
    # Backend: CREST
    # =========================================================================

    def _backend_crest(
        self,
        valid_rows,
        num_confs,
        seed,
        threads,
        xyz_out_dir,
        metadata_index,
        crest_version,
        rdkit_version,
        **kwargs,
    ) -> ConformerSet:
        """
        Single RDKit seed → CREST conformer search.
        NOTE: CREST backend is not validated on the current HPC environment.
        Test independently before production use.
        """
        conformer_set = ConformerSet()
        timestamp     = datetime.utcnow().isoformat() + "Z"

        exe = self.config.get("crest", {}).get("executable")
        if not exe or not os.path.isfile(exe):
            self.fail(f"CREST executable not found: {exe}")

        for idx, inchi_key, smiles, row in valid_rows:
            self.log_section(f"CREST  {inchi_key}")

            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            emb = AllChem.ETKDGv3()
            emb.randomSeed           = seed
            emb.numThreads           = threads
            emb.useSmallRingTorsions = True
            emb.enforceChirality     = True

            try:
                conf_ids = AllChem.EmbedMultipleConfs(mol, 1, emb)
            except Exception as e:
                self.log_warning(f"{inchi_key}: seed embedding failed — {e}")
                if self.strict_mode:
                    self.fail(str(e))
                self.warnings["embedding_failures"] += 1
                self.update_progress(inchi_key, success=False)
                continue

            if len(conf_ids) == 0:
                self.log_warning(f"{inchi_key}: no seed conformer embedded")
                if self.strict_mode:
                    self.fail(f"No seed conformer for {inchi_key}")
                self.warnings["zero_conformers"] += 1
                self.update_progress(inchi_key, success=False)
                continue

            charge, spin = self._get_charge_spin(inchi_key, metadata_index)

            workdir   = os.path.join(self.outputs_dir, "crest_work", inchi_key)
            os.makedirs(workdir, exist_ok=True)
            input_xyz = os.path.join(workdir, "input.xyz")
            self._write_xyz(mol, conf_ids[0], input_xyz)

            gfn_mode = self._resolve_gfn_mode()
            gfn_flag = "--gfn0" if gfn_mode == "gfn0" else "--gfn2"

            cmd = [
                exe,
                input_xyz,
                gfn_flag,
                "--nci",
                "--nconf", str(num_confs),
                "--ewin",  "6",
                "--chrg",  str(charge),
                "--uhf",   str(max(0, spin - 1)),
                "--nthreads", str(threads),
            ]
            self.log_info(f"CREST cmd: {' '.join(cmd)}")

            try:
                result = subprocess.run(
                    cmd, cwd=workdir,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                )
            except Exception as e:
                self.log_error(f"{inchi_key}: CREST execution error — {e}")
                if self.strict_mode:
                    self.fail(str(e))
                self.warnings["crest_failures"] += 1
                self.update_progress(inchi_key, success=False)
                continue

            if result.returncode != 0:
                self.log_error(f"{inchi_key}: CREST non-zero exit")
                self.log_debug(f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}")
                if self.strict_mode:
                    self.fail(f"CREST failed for {inchi_key}")
                self.warnings["crest_failures"] += 1
                self.update_progress(inchi_key, success=False)
                continue

            crest_xyz      = os.path.join(workdir, "crest_conformers.xyz")
            crest_energies = os.path.join(workdir, "crest.energies")

            if not os.path.isfile(crest_xyz) or not os.path.isfile(crest_energies):
                self.log_error(f"{inchi_key}: CREST output files missing")
                if self.strict_mode:
                    self.fail(f"CREST outputs missing for {inchi_key}")
                self.warnings["crest_failures"] += 1
                self.update_progress(inchi_key, success=False)
                continue

            energies = self._parse_crest_energies(crest_energies)
            records  = self._split_crest_xyz(
                crest_xyz    = crest_xyz,
                energies     = energies,
                inchi_key    = inchi_key,
                smiles       = smiles,
                xyz_out_dir  = xyz_out_dir,
                timestamp    = timestamp,
                source_row   = int(idx),
                charge       = charge,
                spin         = spin,
                gfn_mode     = gfn_mode,
                crest_version = crest_version,
            )

            if not records:
                self.log_warning(f"{inchi_key}: no conformers parsed from CREST output")
                if self.strict_mode:
                    self.fail(f"No conformers from CREST for {inchi_key}")
                self.warnings["zero_conformers"] += 1
                self.update_progress(inchi_key, success=False)
                continue

            for rec in records:
                conformer_set.add(rec)

            self.log_info(f"{inchi_key}: {len(records)} CREST conformers")
            self.update_progress(inchi_key)

        return conformer_set

    def _resolve_gfn_mode(self) -> str:
        """gfn0 | gfn2.  Reads parameters["gfn"] or parameters["level"]."""
        params = self.parameters
        gfn    = params.get("gfn")
        if isinstance(gfn, str) and gfn.lower() in ("gfn0", "gfn2"):
            return gfn.lower()
        level = params.get("level", "")
        if isinstance(level, str) and level.lower() == "fast":
            return "gfn0"
        return "gfn2"

    def _parse_crest_energies(self, path: str) -> list:
        energies = []
        try:
            with open(path) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            energies.append(float(parts[1]))
                        except ValueError:
                            pass
        except Exception as e:
            self.log_warning(f"Failed to parse crest.energies: {e}")
        return energies

    def _split_crest_xyz(
        self,
        crest_xyz,
        energies,
        inchi_key,
        smiles,
        xyz_out_dir,
        timestamp,
        source_row,
        charge,
        spin,
        gfn_mode,
        crest_version,
    ) -> list:
        records = []
        try:
            with open(crest_xyz) as f:
                lines = [l.rstrip("\n") for l in f]
        except Exception as e:
            self.log_warning(f"Failed to read crest_conformers.xyz: {e}")
            return records

        i        = 0
        conf_num = 0
        n_lines  = len(lines)

        while i < n_lines:
            try:
                natoms = int(lines[i].strip())
            except ValueError:
                break

            if i + 1 >= n_lines:
                break

            energy_line = lines[i + 1].strip()
            atom_block  = lines[i + 2: i + 2 + natoms]

            if len(atom_block) < natoms:
                break

            xyz_path = os.path.join(
                xyz_out_dir, f"{inchi_key}_conf{conf_num:03d}.xyz"
            )
            try:
                with open(xyz_path, "w") as out:
                    out.write(f"{natoms}\n{energy_line}\n")
                    for al in atom_block:
                        out.write(al + "\n")
            except Exception as e:
                self.log_warning(
                    f"{inchi_key} conf{conf_num:03d}: write failed — {e}"
                )
                break

            energy = energies[conf_num] if conf_num < len(energies) else None

            record = ConformerRecord(
                inchi_key  = inchi_key,
                conf_num   = conf_num,
                xyz_path   = xyz_path,
                energy     = energy,
                smiles     = smiles,
                provenance = {
                    "backend":              "crest",
                    "crest_version":        crest_version,
                    "gfn":                  gfn_mode,
                    "generation_timestamp": timestamp,
                    "source_row":           source_row,
                    "charge":               charge,
                    "spin":                 spin,
                },
            )
            record.optimisation_history.append({
                "stage":     "generation",
                "engine":    "crest",
                "energy":    energy,
                "xyz_path":  xyz_path,
                "timestamp": timestamp,
            })
            records.append(record)
            conf_num += 1
            i += 2 + natoms

        return records

    # =========================================================================
    # Backend: OpenBabel
    # =========================================================================

    def _backend_openbabel(
        self,
        valid_rows,
        num_confs,
        seed,
        threads,
        xyz_out_dir,
        metadata_index,
        rdkit_version,
        openbabel_version,
        **kwargs,
    ) -> ConformerSet:
        """
        RDKit seed → OBConformerSearch GA.
        MMFF94 energy evaluated after conformer generation (no geometry update).
        NOTE: OpenBabel backend is not validated on the current HPC environment.
        Test independently before production use.
        """
        if not HAVE_OPENBABEL:
            self.fail(
                "OpenBabel backend requested but openbabel is not installed."
            )

        import numpy as np
        from openbabel import openbabel as ob

        conformer_set = ConformerSet()
        timestamp     = datetime.utcnow().isoformat() + "Z"

        for idx, inchi_key, smiles, row in valid_rows:
            self.log_section(f"OpenBabel  {inchi_key}")

            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            emb = AllChem.ETKDGv3()
            emb.randomSeed           = seed
            emb.numThreads           = threads
            emb.useSmallRingTorsions = True
            emb.enforceChirality     = True

            try:
                conf_ids = AllChem.EmbedMultipleConfs(mol, 1, emb)
            except Exception as e:
                self.log_warning(f"{inchi_key}: seed embedding failed — {e}")
                if self.strict_mode:
                    self.fail(str(e))
                self.update_progress(inchi_key, success=False)
                continue

            if len(conf_ids) == 0:
                self.log_warning(f"{inchi_key}: no seed conformer for OpenBabel")
                self.update_progress(inchi_key, success=False)
                continue

            workdir  = os.path.join(self.outputs_dir, "openbabel_work", inchi_key)
            os.makedirs(workdir, exist_ok=True)
            seed_xyz = os.path.join(workdir, "input.xyz")
            self._write_xyz(mol, conf_ids[0], seed_xyz)

            # Load seed into OB
            obmol = ob.OBMol()
            conv  = ob.OBConversion()
            conv.SetInAndOutFormats("xyz", "xyz")
            conv.ReadFile(obmol, seed_xyz)

            # GA conformer search
            search = ob.OBConformerSearch()
            search.Setup(obmol, int(num_confs))
            search.SetScore(ob.OBRMSDConformerScore())
            search.SetFilter(ob.OBStericConformerFilter())
            search.Search()
            search.GetConformers(obmol)

            n_generated = obmol.NumConformers()
            self.log_info(f"{inchi_key}: {n_generated} OB conformers generated")

            if n_generated == 0:
                self.log_warning(f"{inchi_key}: OpenBabel produced zero conformers")
                if self.strict_mode:
                    self.fail(f"Zero conformers from OpenBabel for {inchi_key}")
                self.warnings["openbabel_failures"] += 1
                self.update_progress(inchi_key, success=False)
                continue

            # Write OB conformers to temp XYZs
            ob_xyz_dir = os.path.join(workdir, "xyz")
            os.makedirs(ob_xyz_dir, exist_ok=True)
            for cidx in range(n_generated):
                obmol.SetConformer(cidx)
                conv.WriteFile(obmol, os.path.join(ob_xyz_dir, f"c{cidx:03d}.xyz"))

            # Rebuild RDKit mol with OB geometries, evaluate MMFF94 energies
            rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            for cidx in range(n_generated):
                m2 = Chem.MolFromXYZFile(os.path.join(ob_xyz_dir, f"c{cidx:03d}.xyz"))
                if m2 is not None:
                    rdkit_mol.AddConformer(m2.GetConformer(), assignId=True)

            # CalcEnergy only — no geometry update
            mmff_results = AllChem.MMFFOptimizeMoleculeConfs(
                rdkit_mol, numThreads=threads
            )
            energies = np.array([r[1] for r in mmff_results], dtype=float)

            charge, spin = self._get_charge_spin(inchi_key, metadata_index)

            for conf_num in range(n_generated):
                xyz_path = os.path.join(
                    xyz_out_dir, f"{inchi_key}_conf{conf_num:03d}.xyz"
                )
                self._write_xyz(rdkit_mol, conf_num, xyz_path)

                record = ConformerRecord(
                    inchi_key  = inchi_key,
                    conf_num   = conf_num,
                    xyz_path   = xyz_path,
                    energy     = float(energies[conf_num]),
                    smiles     = smiles,
                    provenance = {
                        "backend":              "openbabel",
                        "openbabel_version":    openbabel_version,
                        "rdkit_version":        rdkit_version,
                        "seed":                 seed,
                        "threads":              threads,
                        "generation_timestamp": timestamp,
                        "source_row":           int(idx),
                        "charge":               charge,
                        "spin":                 spin,
                    },
                )
                record.optimisation_history.append({
                    "stage":     "generation",
                    "engine":    "openbabel",
                    "energy":    float(energies[conf_num]),
                    "xyz_path":  xyz_path,
                    "timestamp": timestamp,
                })
                conformer_set.add(record)

            self.update_progress(inchi_key)

        return conformer_set

    # =========================================================================
    # XYZ writer
    # =========================================================================

    def _write_xyz(self, mol, conf_id: int, path: str):
        conf  = mol.GetConformer(conf_id)
        atoms = list(mol.GetAtoms())
        with open(path, "w") as f:
            f.write(f"{len(atoms)}\nGenerated by GenerationStage\n")
            for atom in atoms:
                pos = conf.GetAtomPosition(atom.GetIdx())
                f.write(
                    f"{atom.GetSymbol()} "
                    f"{pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n"
                )

    # =========================================================================
    # Output writing
    # =========================================================================

    def _write_outputs(
        self,
        conformer_set: ConformerSet,
        xyz_dir:       str,
        num_valid:     int,
        total_rows:    int,
        backend:       str,
        num_confs:     int,
        seed:          int,
    ):
        # summary.csv
        summary_path = os.path.join(self.outputs_dir, "summary.csv")
        with AtomicWriter(summary_path) as f:
            pd.DataFrame(conformer_set.to_list()).to_csv(f, index=False)

        # energies.json  (canonical stage output)
        energies_path = self.get_stage_output()
        conformer_set.save(energies_path)

        # generation_metadata.json
        meta_path = os.path.join(self.outputs_dir, "generation_metadata.json")
        with AtomicWriter(meta_path) as f:
            json.dump(
                {
                    "stage":                 "generation",
                    "backend":               backend,
                    "num_valid_molecules":   num_valid,
                    "num_total_rows":        total_rows,
                    "num_conformers":        len(conformer_set),
                    "num_confs_requested":   num_confs,
                    "seed":                  seed,
                    "timestamp":             datetime.utcnow().isoformat() + "Z",
                },
                f, indent=2,
            )

        self.log_info(
            f"Outputs: energies.json ({len(conformer_set)} conformers), "
            f"summary.csv, generation_metadata.json"
        )

    # =========================================================================
    # Warning summary
    # =========================================================================

    def _log_warning_summary(self):
        total = sum(self.warnings.values())
        if total == 0:
            self.log_info("Generation completed with no warnings.")
            return

        self.log_info("Generation completed with warnings:")
        labels = {
            "embedding_failures":   "Embedding failures",
            "energy_eval_failures": "Energy evaluation failures",
            "zero_conformers":      "Zero conformers generated",
            "crest_failures":       "CREST failures",
            "openbabel_failures":   "OpenBabel failures",
        }
        for key, label in labels.items():
            if self.warnings[key]:
                self.log_warning(f"{label}: {self.warnings[key]}")