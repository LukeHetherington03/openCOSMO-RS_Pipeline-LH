import os
import json
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
    "rdkit": "_backend_rdkit",
    "crest": "_backend_crest",
    "openbabel": "_backend_openbabel",
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

        # Copy cleaned.csv into inputs/ for reproducibility
        inputs_cleaned = os.path.join(self.inputs_dir, "cleaned.csv")
        if os.path.abspath(stage_input) != os.path.abspath(inputs_cleaned):
            shutil.copy(stage_input, inputs_cleaned)
        self.log(f"Copied cleaned.csv → {inputs_cleaned}")

        params = self.parameters
        engine = params.get("engine", "rdkit").lower()
        num_confs = params.get("n")
        seed = params.get("seed", 42)
        threads = params.get("threads", 1)

        self.log(f"Backend: {engine}")
        self.log(f"Stage input: {stage_input}")
        self.log(f"Seed: {seed}")
        self.log(f"Threads: {threads}")
        self.log(f"Strict mode: {self.strict_mode}")

        # Prepare directories
        dirs = self._prepare_directories()

        # Load cleaned dataset
        df_all = self._load_csv(stage_input)

        # Load metadata (charge, spin, etc.)
        metadata_index = self._load_metadata_index()

        # Validate rows
        valid_rows = self._validate_rows(df_all)
        self.set_items([row_idx for row_idx, *_ in valid_rows])

        # If n not provided, use rotatable-bond heuristic
        if num_confs is None:
            num_confs = self._default_num_confs_from_rotatable_bonds(valid_rows)
            self.log(f"Using heuristic conformer count: n={num_confs}")
        else:
            self.log(f"Conformers per molecule (explicit): {num_confs}")

        self.warnings = {
            "embedding_failures": 0,
            "optimisation_failures": 0,
            "zero_conformers": 0,
            "crest_failures": 0,
            "openbabel_failures": 0,
        }

        # Cache versions for provenance
        rdkit_version = getattr(Chem, "__version__", None)
        crest_version = self._get_crest_version()
        openbabel_version = self._get_openbabel_version() if HAVE_OPENBABEL else None

        # Dispatch backend
        conformer_set = self._dispatch_backend(
            backend=engine,
            valid_rows=valid_rows,
            num_confs=num_confs,
            seed=seed,
            threads=threads,
            xyz_out_dir=dirs["xyz"],
            metadata_index=metadata_index,
            rdkit_version=rdkit_version,
            crest_version=crest_version,
            openbabel_version=openbabel_version,
        )

        if len(conformer_set) == 0:
            self.fail("Generation produced zero conformers.")

        # Write outputs
        self._write_outputs(
            conformer_set=conformer_set,
            dirs=dirs,
            num_valid=len(valid_rows),
            total_rows=len(df_all),
            backend=engine,
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
    # Helpers: Metadata loading
    # ------------------------------------------------------------
    def _load_metadata_index(self):
        """
        Load molecule metadata from constant_files.metadata_dir.

        Expects files named {inchi_key}.json containing at least:
            {
              "charge": 0,
              "spin": 1,
              ...
            }
        """
        metadata_index = {}

        constant_files = self.config.get("constant_files", {})
        metadata_dir = constant_files.get("metadata_dir")

        if not metadata_dir or not os.path.isdir(metadata_dir):
            self.log(
                f"[INFO] No metadata_dir found or directory missing: {metadata_dir}. "
                "Using default charge=0, spin=1."
            )
            return metadata_index

        self.log(f"Loading molecule metadata from: {metadata_dir}")

        for fname in os.listdir(metadata_dir):
            if not fname.endswith(".json"):
                continue
            inchi_key = fname[:-5]  # strip .json
            fpath = os.path.join(metadata_dir, fname)
            try:
                with open(fpath, "r") as f:
                    metadata = json.load(f)
                metadata_index[inchi_key] = metadata
            except Exception as e:
                self.log(f"[WARNING] Failed to load metadata for {inchi_key}: {e}")

        self.log(f"Loaded metadata for {len(metadata_index)} molecules")
        return metadata_index

    def _get_charge_spin(self, inchi_key, params, metadata_index):
        """
        Resolve charge and spin for a molecule.

        Precedence:
            1. Stage parameters (global override)
            2. Metadata JSON (per-molecule)
            3. Defaults: charge=0, spin=1
        """
        default_charge = 0
        default_spin = 1

        # Stage-level overrides
        param_charge = params.get("charge")
        param_spin = params.get("spin")

        meta = metadata_index.get(inchi_key, {})
        meta_charge = meta.get("charge")
        meta_spin = meta.get("spin")

        charge = (
            param_charge
            if param_charge is not None
            else (meta_charge if meta_charge is not None else default_charge)
        )
        spin = (
            param_spin
            if param_spin is not None
            else (meta_spin if meta_spin is not None else default_spin)
        )

        return int(charge), int(spin)

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
    # Helpers: default conformer count heuristic
    # ------------------------------------------------------------
    def _default_num_confs_from_rotatable_bonds(self, valid_rows):
        """
        Use O'Boyle 2011 heuristic based on number of rotatable bonds:
            <= 7  -> 50
            8-12  -> 200
            >= 13 -> 300

        We compute this per molecule, but for simplicity we return a single
        global n as the max over all molecules (safe upper bound).
        """
        max_confs = 0
        for _, _, smiles, _ in valid_rows:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            try:
                mol_no_h = Chem.RemoveHs(mol)
            except Exception:
                mol_no_h = mol
            try:
                n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol_no_h)
            except Exception:
                n_rot = 8  # middle heuristic if something goes wrong

            if n_rot <= 7:
                n = 50
            elif 8 <= n_rot <= 12:
                n = 200
            else:
                n = 300

            if n > max_confs:
                max_confs = n

        return max_confs if max_confs > 0 else 50

    # ------------------------------------------------------------
    # Helpers: version detection
    # ------------------------------------------------------------
    def _get_crest_version(self):
        crest_cfg = self.config.get("crest", {})
        crest_exe = crest_cfg.get("executable")
        if not crest_exe or not os.path.isfile(crest_exe):
            return None
        try:
            result = subprocess.run(
                [crest_exe, "--version"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode == 0:
                line = result.stdout.strip().splitlines()[0]
                return line.strip()
        except Exception:
            return None
        return None

    def _get_openbabel_version(self):
        if not HAVE_OPENBABEL:
            return None
        try:
            return ob.OBReleaseVersion()
        except Exception:
            return None

    # ------------------------------------------------------------
    # Backend dispatch
    # ------------------------------------------------------------
    def _dispatch_backend(
        self,
        backend,
        valid_rows,
        num_confs,
        seed,
        threads,
        xyz_out_dir,
        metadata_index,
        rdkit_version,
        crest_version,
        openbabel_version,
    ):
        if backend not in GENERATION_BACKENDS:
            self.fail(f"Unknown backend: {backend}")

        runner_name = GENERATION_BACKENDS[backend]
        runner = getattr(self, runner_name)

        self.log_section(f"Running backend: {backend}")
        return runner(
            valid_rows=valid_rows,
            num_confs=num_confs,
            seed=seed,
            threads=threads,
            xyz_out_dir=xyz_out_dir,
            metadata_index=metadata_index,
            rdkit_version=rdkit_version,
            crest_version=crest_version,
            openbabel_version=openbabel_version,
        )

    # ------------------------------------------------------------
    # Backend: RDKit (pure generation)
    # ------------------------------------------------------------
    def _backend_rdkit(
        self,
        valid_rows,
        num_confs,
        seed,
        threads,
        xyz_out_dir,
        metadata_index,
        rdkit_version,
        crest_version,
        openbabel_version,
    ):
        conformer_set = ConformerSet()
        timestamp = datetime.utcnow().isoformat() + "Z"

        for idx, inchi_key, smiles, row in valid_rows:
            self.log_section(f"Generating conformers (RDKit) for {inchi_key} ({smiles})")

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
            params.numThreads = threads
            params.useSmallRingTorsions = True
            params.useMacrocycleTorsions = True
            params.enforceChirality = True

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

            # RDKit forcefield: MMFF94 by default (pure generation energy)
            mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
            use_mmff = mmff_props is not None

            for local_conf_idx, conf_id in enumerate(conf_ids):
                conf_num = int(local_conf_idx)

                try:
                    if use_mmff:
                        ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf_id)
                    else:
                        AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
                        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)

                    energy = float(ff.CalcEnergy())
                except Exception as e:
                    msg = f"RDKit energy evaluation failed for {inchi_key} conf {conf_num}: {e}"
                    self.log(f"[WARNING] {msg}")
                    if self.strict_mode:
                        self.fail(msg)
                    self.warnings["optimisation_failures"] += 1
                    continue

                xyz_path = os.path.join(
                    xyz_out_dir,
                    f"{inchi_key}_conf{conf_num:03d}.xyz",
                )

                charge, spin = self._get_charge_spin(inchi_key, self.parameters, metadata_index)

                record = ConformerRecord(
                    inchi_key=inchi_key,
                    conf_num=conf_num,
                    xyz_path=xyz_path,
                    energy=energy,
                    smiles=smiles,
                    provenance={
                        "backend": "rdkit",
                        "rdkit_version": rdkit_version,
                        "seed": seed,
                        "threads": threads,
                        "forcefield": "MMFF94" if use_mmff else "UFF",
                        "generation_timestamp": timestamp,
                        "source_row": int(idx),
                        "charge": charge,
                        "spin": spin,
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
    # Backend: CREST (pure generation)
    # ------------------------------------------------------------
    def _backend_crest(
        self,
        valid_rows,
        num_confs,
        seed,
        threads,
        xyz_out_dir,
        metadata_index,
        rdkit_version,
        crest_version,
        openbabel_version,
    ):
        conformer_set = ConformerSet()
        timestamp = datetime.utcnow().isoformat() + "Z"

        crest_cfg = self.config.get("crest", {})
        crest_exe = crest_cfg.get("executable")
        if not crest_exe or not os.path.isfile(crest_exe):
            self.fail(f"CREST executable not found or invalid: {crest_exe}")

        for idx, inchi_key, smiles, row in valid_rows:
            self.log_section(f"Generating conformers (CREST) for {inchi_key} ({smiles})")

            base_mol = Chem.MolFromSmiles(smiles)
            if base_mol is None:
                msg = f"Row {idx}: RDKit failed to parse SMILES '{smiles}' for CREST input"
                self.log(f"[WARNING] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.warnings["embedding_failures"] += 1
                self.update_progress(idx, success=False)
                continue

            mol = Chem.AddHs(base_mol)

            # Single 3D conformer for CREST input (no pre-optimisation here)
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            params.numThreads = threads
            params.useSmallRingTorsions = True
            params.useMacrocycleTorsions = True
            params.enforceChirality = True

            try:
                conf_ids = AllChem.EmbedMultipleConfs(mol, 1, params)
            except Exception as e:
                msg = f"RDKit embedding failed for CREST input {inchi_key}: {e}"
                self.log(f"[WARNING] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.warnings["embedding_failures"] += 1
                self.update_progress(idx, success=False)
                continue

            if len(conf_ids) == 0:
                msg = f"No initial conformer generated for CREST input {inchi_key}"
                self.log(f"[WARNING] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.warnings["zero_conformers"] += 1
                self.update_progress(idx, success=False)
                continue

            # Working directory per molecule
            workdir = os.path.join(self.outputs_dir, "crest_work", inchi_key)
            os.makedirs(workdir, exist_ok=True)

            input_xyz = os.path.join(workdir, "input.xyz")
            self._write_xyz(mol, conf_ids[0], input_xyz)

            # Resolve charge and spin
            charge, spin = self._get_charge_spin(inchi_key, self.parameters, metadata_index)

            # Resolve GFN level
            gfn_mode = self._resolve_gfn_mode(self.parameters)
            if gfn_mode == "gfn0":
                gfn_flag = "--gfn0"
            else:
                gfn_flag = "--gfn2"

            # CREST command (pure generation)
            cmd = [
                crest_exe,
                input_xyz,
                gfn_flag,
                "--nci",
                "--nconf", str(num_confs),
                "--ewin", "6",
                "--chrg", str(charge),
                "--uhf", str(max(0, spin - 1)),
                "--nthreads", str(threads),
            ]

            self.log(f"[CREST] Command: {' '.join(cmd)}")
            self.log(f"[CREST] Working directory: {workdir}")

            try:
                result = subprocess.run(
                    cmd,
                    cwd=workdir,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
            except Exception as e:
                msg = f"CREST execution failed for {inchi_key}: {e}"
                self.log(f"[ERROR] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.warnings["crest_failures"] += 1
                self.update_progress(idx, success=False)
                continue

            if result.returncode != 0:
                self.log(f"[ERROR] CREST non-zero exit for {inchi_key}")
                self.log(f"[CREST stdout]\n{result.stdout}")
                self.log(f"[CREST stderr]\n{result.stderr}")
                if self.strict_mode:
                    self.fail(f"CREST failed for {inchi_key}")
                self.warnings["crest_failures"] += 1
                self.update_progress(idx, success=False)
                continue

            crest_xyz = os.path.join(workdir, "crest_conformers.xyz")
            crest_energies = os.path.join(workdir, "crest.energies")

            if not os.path.isfile(crest_xyz) or not os.path.isfile(crest_energies):
                msg = f"CREST outputs missing for {inchi_key} (crest_conformers.xyz / crest.energies)"
                self.log(f"[ERROR] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.warnings["crest_failures"] += 1
                self.update_progress(idx, success=False)
                continue

            # Parse energies
            energies = self._parse_crest_energies(crest_energies)

            # Parse XYZ and write per-conformer files
            conf_records = self._split_crest_xyz(
                crest_xyz=crest_xyz,
                energies=energies,
                inchi_key=inchi_key,
                smiles=smiles,
                xyz_out_dir=xyz_out_dir,
                timestamp=timestamp,
                idx=idx,
                charge=charge,
                spin=spin,
                gfn_mode=gfn_mode,
                crest_version=crest_version,
            )

            if not conf_records:
                msg = f"No conformers parsed from CREST outputs for {inchi_key}"
                self.log(f"[WARNING] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.warnings["zero_conformers"] += 1
                self.update_progress(idx, success=False)
                continue

            for rec in conf_records:
                conformer_set.add(rec)

            self.update_progress(idx)

        return conformer_set

    def _resolve_gfn_mode(self, params):
        """
        Resolve GFN level for CREST.

        Priority:
            1. params["gfn"] if present ("gfn0" or "gfn2")
            2. params["level"] if present ("fast" -> gfn0, "normal" -> gfn2)
            3. default: gfn2
        """
        gfn = params.get("gfn")
        if isinstance(gfn, str):
            gfn = gfn.lower()
            if gfn in ("gfn0", "gfn2"):
                return gfn

        level = params.get("level")
        if isinstance(level, str):
            level = level.lower()
            if level == "fast":
                return "gfn0"
            if level == "normal":
                return "gfn2"

        return "gfn2"

    def _parse_crest_energies(self, path):
        """
        Parse crest.energies file.

        Typically contains lines with at least:
            index  energy(Ha)  ...
        Returns a list of energies (float).
        """
        energies = []
        try:
            with open(path, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) < 2:
                        continue
                    try:
                        e = float(parts[1])
                        energies.append(e)
                    except ValueError:
                        continue
        except Exception as e:
            self.log(f"[WARNING] Failed to parse crest.energies: {e}")
        return energies

    def _split_crest_xyz(
        self,
        crest_xyz,
        energies,
        inchi_key,
        smiles,
        xyz_out_dir,
        timestamp,
        idx,
        charge,
        spin,
        gfn_mode,
        crest_version,
    ):
        records = []

        try:
            with open(crest_xyz, "r") as f:
                lines = [l.rstrip("\n") for l in f]
        except Exception as e:
            self.log(f"[WARNING] Failed to read crest_conformers.xyz: {e}")
            return records

        i = 0
        conf_num = 0
        n_lines = len(lines)

        while i < n_lines:
            # 1. natoms
            try:
                natoms = int(lines[i].strip())
            except ValueError:
                self.log(f"[WARNING] Unexpected line in CREST XYZ: {lines[i]}")
                break

            # 2. energy line
            if i + 1 >= n_lines:
                break
            energy_line = lines[i + 1].strip()

            # 3. atom block
            start = i + 2
            end = start + natoms
            if end > n_lines:
                break

            atom_block = lines[start:end]

            # Write XYZ
            xyz_path = os.path.join(
                xyz_out_dir,
                f"{inchi_key}_conf{conf_num:03d}.xyz",
            )

            try:
                with open(xyz_path, "w") as out:
                    out.write(f"{natoms}\n")
                    out.write(f"{energy_line}\n")
                    for atom_line in atom_block:
                        out.write(atom_line + "\n")
            except Exception as e:
                self.log(f"[WARNING] Failed to write CREST XYZ for {inchi_key} conf {conf_num}: {e}")
                break

            # Energy from crest.energies if available
            energy = energies[conf_num] if conf_num < len(energies) else None

            record = ConformerRecord(
                inchi_key=inchi_key,
                conf_num=conf_num,
                xyz_path=xyz_path,
                energy=energy,
                smiles=smiles,
                provenance={
                    "backend": "crest",
                    "crest_version": crest_version,
                    "gfn": gfn_mode,
                    "generation_timestamp": timestamp,
                    "source_row": int(idx),
                    "charge": charge,
                    "spin": spin,
                },
            )

            record.optimisation_history.append({
                "stage": "generation",
                "engine": "crest",
                "energy": energy,
                "xyz_path": xyz_path,
                "timestamp": timestamp,
            })

            records.append(record)
            conf_num += 1

            # Advance to next block
            i = end

        self.log(f"[CREST] Parsed {len(records)} conformers for {inchi_key}")
        return records


    # ------------------------------------------------------------
    # Backend: OpenBabel (pure generation)
    # ------------------------------------------------------------
    def _backend_openbabel(
        self,
        valid_rows,
        num_confs,
        seed,
        threads,
        xyz_out_dir,
        metadata_index,
        rdkit_version,
        crest_version,
        openbabel_version,
    ):
        if not HAVE_OPENBABEL:
            self.fail("OpenBabel backend requested but openbabel-wheel is not installed/importable.")

        from openbabel import openbabel as ob
        import numpy as np

        conformer_set = ConformerSet()
        timestamp = datetime.utcnow().isoformat() + "Z"

        for idx, inchi_key, smiles, row in valid_rows:
            self.log_section(f"Generating conformers (OpenBabel) for {inchi_key} ({smiles})")

            # RDKit: initial 3D seed
            base_mol = Chem.MolFromSmiles(smiles)
            if base_mol is None:
                msg = f"Row {idx}: RDKit failed to parse SMILES '{smiles}'"
                self.log(f"[WARNING] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.update_progress(idx, success=False)
                continue

            mol_seed = Chem.AddHs(base_mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            params.numThreads = threads
            params.useSmallRingTorsions = True
            params.useMacrocycleTorsions = True
            params.enforceChirality = True

            try:
                conf_ids = AllChem.EmbedMultipleConfs(mol_seed, 1, params)
            except Exception as e:
                msg = f"RDKit embedding failed for OpenBabel seed {inchi_key}: {e}"
                self.log(f"[WARNING] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.update_progress(idx, success=False)
                continue

            if len(conf_ids) == 0:
                msg = f"No initial conformer generated for OpenBabel seed {inchi_key}"
                self.log(f"[WARNING] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.update_progress(idx, success=False)
                continue

            workdir = os.path.join(self.outputs_dir, "openbabel_work", inchi_key)
            os.makedirs(workdir, exist_ok=True)

            seed_xyz = os.path.join(workdir, "input.xyz")
            self._write_xyz(mol_seed, conf_ids[0], seed_xyz)

            # OpenBabel: load seed
            obmol = ob.OBMol()
            conv = ob.OBConversion()
            conv.SetInAndOutFormats("xyz", "xyz")
            conv.ReadFile(obmol, seed_xyz)

            # OpenBabel GA search (no FD tricks)
            search = ob.OBConformerSearch()
            search.Setup(obmol, int(num_confs))
            search.SetScore(ob.OBRMSDConformerScore())
            search.SetFilter(ob.OBStericConformerFilter())

            search.Search()
            search.GetConformers(obmol)

            n_generated = obmol.NumConformers()
            self.log(f"[OpenBabel] Generated {n_generated} conformers for {inchi_key}")

            if n_generated == 0:
                msg = f"OpenBabel generated zero conformers for {inchi_key}"
                self.log(f"[WARNING] {msg}")
                if self.strict_mode:
                    self.fail(msg)
                self.update_progress(idx, success=False)
                continue

            # Write OBabel conformers to temp XYZs
            ob_xyz_dir = os.path.join(workdir, "xyz")
            os.makedirs(ob_xyz_dir, exist_ok=True)

            for cidx in range(n_generated):
                obmol.SetConformer(cidx)
                conv.WriteFile(obmol, os.path.join(ob_xyz_dir, f"{inchi_key}_c{cidx:03d}.xyz"))

            # RDKit container + MMFF energies (no pruning)
            rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            for cidx in range(n_generated):
                xyz_path = os.path.join(ob_xyz_dir, f"{inchi_key}_c{cidx:03d}.xyz")
                m2 = Chem.MolFromXYZFile(xyz_path)
                if m2 is not None:
                    rdkit_mol.AddConformer(m2.GetConformer(), assignId=True)

            mmff_results = AllChem.MMFFOptimizeMoleculeConfs(rdkit_mol, numThreads=threads)
            energies = np.array([res[1] for res in mmff_results], dtype=float)

            # Final XYZs + records, IDs stable
            for conf_num in range(n_generated):
                xyz_path = os.path.join(xyz_out_dir, f"{inchi_key}_conf{conf_num:03d}.xyz")
                self._write_xyz(rdkit_mol, conf_num, xyz_path)

                charge, spin = self._get_charge_spin(inchi_key, self.parameters, metadata_index)

                record = ConformerRecord(
                    inchi_key=inchi_key,
                    conf_num=conf_num,
                    xyz_path=xyz_path,
                    energy=float(energies[conf_num]),
                    smiles=smiles,
                    provenance={
                        "backend": "openbabel",
                        "openbabel_version": openbabel_version,
                        "rdkit_version": rdkit_version,
                        "seed": seed,
                        "threads": threads,
                        "generation_timestamp": timestamp,
                        "source_row": int(idx),
                        "charge": charge,
                        "spin": spin,
                    },
                )

                record.optimisation_history.append({
                    "stage": "generation",
                    "engine": "openbabel",
                    "energy": float(energies[conf_num]),
                    "xyz_path": xyz_path,
                    "timestamp": timestamp,
                })

                conformer_set.add(record)

            self.update_progress(idx)

        return conformer_set


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

        self.log(f"Summary written to: {summary_path}")
        self.log(f"Energies written to: {energies_path}")
        self.log(f"Generation metadata written to: {meta_path}")

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
                f" - Energy evaluation failures: {self.warnings['optimisation_failures']} conformer(s)."
            )
        if self.warnings["zero_conformers"] > 0:
            self.log(
                f" - Zero conformers generated for "
                f"{self.warnings['zero_conformers']} molecule(s)."
            )
        if self.warnings["crest_failures"] > 0:
            self.log(
                f" - CREST failures: {self.warnings['crest_failures']} molecule(s)."
            )
        if self.warnings["openbabel_failures"] > 0:
            self.log(
                f" - OpenBabel failures: {self.warnings['openbabel_failures']} molecule(s)."
            )
