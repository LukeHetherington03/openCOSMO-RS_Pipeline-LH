import os
import shutil
import subprocess
import json

from rdkit import Chem
from rdkit.Chem import rdmolops

from modules.stages.base_stage import BaseStage
from modules.parsers.cosmo.orca_log_parser import OrcaLogParser
from modules.utils.atomic_write import AtomicWriter


class OrcacosmoStage(BaseStage):
    """
    ORCA COSMO stage (two-stage version):
      - Stage 1: CPCM BP86 def2-TZVP PAL8 SP
      - Stage 2: CPCM BP86 def2-TZVPD SP + ELPROP
      - Loads CPCM radii from external file
      - Parses ONLY energy + dipole from log
      - Computes adjacency via RDKit
      - Produces .orcacosmo by concatenating raw ORCA outputs
      - Generates orcacosmo_summary.json
    """

    # ------------------------------------------------------------
    # Entry point
    # ------------------------------------------------------------
    def execute(self):
        self.log_header("Starting ORCA COSMO Stage")

        self.orca_command = self.parameters.get("orca_command", "orca")
        self.cpcm_radii_file = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/CONSTANT_FILES/cpcm/cpcm_radii.inp"

        self._prepare_directories()
        self._load_xyzs_from_summary()

        lookup_ids = self._discover_lookup_ids()
        self.set_items(lookup_ids)

        self.successful_outputs = []

        for lookup_id in list(self.job.pending_items):
            xyz_path = os.path.join(self.inputs_dir, f"{lookup_id}.xyz")

            if not os.path.exists(xyz_path):
                self.log(f"[WARNING] Missing XYZ for {lookup_id}, skipping.")
                self.update_progress(lookup_id, success=False)
                continue

            try:
                self._process_single(lookup_id, xyz_path)
                self.log(f"Completed {lookup_id}")
                self.update_progress(lookup_id)
            except Exception as e:
                self.log(f"[ERROR] {lookup_id}: {e}")
                self.update_progress(lookup_id, success=False)

        self._write_summary()
        self.job.mark_complete()
        self.log_header("ORCA COSMO Stage Complete")

    # ------------------------------------------------------------
    # Load XYZs from summary file
    # ------------------------------------------------------------
    def _load_xyzs_from_summary(self):
        summary_file = self.parameters.get("summary_file")
        if summary_file is None:
            self.fail("ORCA COSMO stage requires summary_file")

        if not os.path.exists(summary_file):
            self.fail(f"summary_file does not exist: {summary_file}")

        with open(summary_file) as f:
            entries = json.load(f)

        for entry in entries:
            lookup_id = entry.get("lookup_id")
            xyz_src = entry.get("xyz_path")

            if not lookup_id or not xyz_src:
                continue

            if not os.path.exists(xyz_src):
                self.log(f"[WARNING] XYZ missing for {lookup_id}: {xyz_src}")
                continue

            dst = os.path.join(self.inputs_dir, f"{lookup_id}.xyz")
            shutil.copy(xyz_src, dst)

    # ------------------------------------------------------------
    # Directory setup
    # ------------------------------------------------------------
    def _prepare_directories(self):
        self.tmp_exec = os.path.join(self.outputs_dir, "tmp_exec")
        self.raw_dir = os.path.join(self.outputs_dir, "raw_outputs")
        self.parsed_dir = os.path.join(self.outputs_dir, "parsed_outputs")

        os.makedirs(self.tmp_exec, exist_ok=True)
        os.makedirs(self.raw_dir, exist_ok=True)
        os.makedirs(self.parsed_dir, exist_ok=True)

    # ------------------------------------------------------------
    # Discover lookup IDs
    # ------------------------------------------------------------
    def _discover_lookup_ids(self):
        ids = []
        for fname in os.listdir(self.inputs_dir):
            if fname.endswith(".xyz"):
                ids.append(os.path.splitext(fname)[0])
        ids.sort()
        self.log(f"Discovered {len(ids)} XYZ structures")
        return ids

    # ------------------------------------------------------------
    # Process a single structure
    # ------------------------------------------------------------
    def _process_single(self, lookup_id, xyz_path):
        workdir = os.path.join(self.tmp_exec, lookup_id)
        os.makedirs(workdir, exist_ok=True)

        # Copy XYZ
        local_xyz = os.path.join(workdir, f"{lookup_id}.xyz")
        shutil.copy(xyz_path, local_xyz)

        # ORCA input + log paths
        inp = os.path.join(workdir, f"{lookup_id}.inp")
        log = os.path.join(workdir, f"{lookup_id}.log")

        # Write two-stage ORCA input
        self._write_orca_input(inp, f"{lookup_id}.xyz", lookup_id)

        # Run ORCA
        self._run_orca(workdir, lookup_id, log)

        # Parse log for energy + dipole
        log_data = self._parse_log(log, lookup_id)

        # Compute adjacency
        adjacency = self._adjacency_from_xyz(xyz_path)

        # Build merged data
        merged = {
            "lookup_id": lookup_id,
            "energy": log_data.get("energy"),
            "dipole": log_data.get("dipole", [0.0, 0.0, 0.0]),
            "adjacency_matrix": adjacency,
            "atomic_polarizabilities": [],  # to be added later
        }

        # Write .orcacosmo (raw concatenation)
        orcacosmo_path = os.path.join(self.parsed_dir, f"{lookup_id}.orcacosmo")
        self._write_orcacosmo(lookup_id, xyz_path, merged, orcacosmo_path)

        # Write JSON summary
        merged_path = os.path.join(self.parsed_dir, f"{lookup_id}_cosmo_full.json")
        with AtomicWriter(merged_path) as f:
            json.dump(merged, f, indent=2)

        # Track success
        self.successful_outputs.append({
            "lookup_id": lookup_id,
            "orcacosmo_path": orcacosmo_path,
            "json_path": merged_path,
            "xyz_path": xyz_path,
            "energy": merged.get("energy"),
        })

        # Clean scratch
        shutil.rmtree(workdir, ignore_errors=True)

    # ------------------------------------------------------------
    # ORCA input writer (two-stage CPCM)
    # ------------------------------------------------------------
    def _write_orca_input(self, path, xyz_name, base):

        # Load CPCM radii block
        with open(self.cpcm_radii_file) as r:
            radii_block = r.read()

        with open(path, "w") as f:
            f.write("%MaxCore 2000\n")
            f.write("%pal nprocs 8 end\n\n")

            # -----------------------------
            # Stage 1: TZVP SP
            # -----------------------------
            f.write(f'%base "{base}_tzvp"\n\n')

            f.write("%cpcm\n")
            f.write(radii_block)
            f.write("end\n\n")

            f.write("! CPCM BP86 def2-TZVP PAL8 SP\n\n")
            f.write(f"* xyzfile 0 1 {xyz_name}\n\n")

            # -----------------------------
            # Stage 2: TZVPD SP + ELPROP
            # -----------------------------
            f.write("$new_job\n\n")
            f.write(f'%base "{base}_tzvpd"\n\n')

            f.write("! CPCM BP86 def2-TZVPD SP\n\n")

            f.write("%elprop\n")
            f.write("  Polar 1\n")
            f.write("  Polaratom 1\n")
            f.write("end\n\n")

            f.write("%cpcm\n")
            f.write(radii_block)
            f.write("end\n\n")

            f.write("* xyzfile 0 1 geo_opt_tzvp.xyz\n")

    # ------------------------------------------------------------
    # Run ORCA
    # ------------------------------------------------------------
    def _run_orca(self, workdir, job_name, log_path):
        with open(log_path, "w") as f:
            subprocess.run(
                [self.orca_command, f"{job_name}.inp"],
                cwd=workdir,
                stdout=f,
                stderr=f,
                check=True,
            )

        # Fix .tmp files
        for suffix in ["tzvp", "tzvpd"]:
            tmp_cpcm = os.path.join(workdir, f"{job_name}_{suffix}.cpcm.tmp")
            final_cpcm = os.path.join(workdir, f"{job_name}_{suffix}.cpcm")
            if os.path.exists(tmp_cpcm) and not os.path.exists(final_cpcm):
                shutil.move(tmp_cpcm, final_cpcm)

            tmp_corr = os.path.join(workdir, f"{job_name}_{suffix}.cpcm_corr.tmp")
            final_corr = os.path.join(workdir, f"{job_name}_{suffix}.cpcm_corr")
            if os.path.exists(tmp_corr) and not os.path.exists(final_corr):
                shutil.move(tmp_corr, final_corr)

        # Copy raw outputs
        outdir = os.path.join(self.raw_dir, "two_stage_jobs")
        os.makedirs(outdir, exist_ok=True)

        for ext in ["inp", "log", "xyz", "gbw", "prop",
                    "tzvp.cpcm", "tzvp.cpcm_corr",
                    "tzvpd.cpcm", "tzvpd.cpcm_corr"]:
            src = os.path.join(workdir, f"{job_name}.{ext}")
            if os.path.exists(src):
                shutil.copy(src, os.path.join(outdir, f"{job_name}.{ext}"))

    # ------------------------------------------------------------
    # Log parser (energy + dipole only)
    # ------------------------------------------------------------
    def _parse_log(self, log_file, lookup_id):
        if not os.path.exists(log_file):
            raise FileNotFoundError(f"Missing ORCA log file: {log_file}")

        try:
            log_data = OrcaLogParser(log_file).parse()
            return {
                "energy": log_data.get("energy"),
                "dipole": log_data.get("dipole", [0.0, 0.0, 0.0]),
            }
        except Exception as e:
            self.log(f"[WARNING] Log parser failed for {lookup_id}: {e}")
            return {"energy": None, "dipole": [0.0, 0.0, 0.0]}

    # ------------------------------------------------------------
    # Adjacency
    # ------------------------------------------------------------
    def _adjacency_from_xyz(self, xyz_path):
        mol = Chem.MolFromXYZFile(xyz_path)
        if mol is None:
            return []
        adj = rdmolops.GetAdjacencyMatrix(mol, useBO=True)
        return adj.tolist()

    # ------------------------------------------------------------
    # Raw concatenation writer
    # ------------------------------------------------------------
    def _write_orcacosmo(self, lookup_id, xyz_path, data, out_path):
        atoms = self._load_xyz(xyz_path)

        # Use TZVPD files
        base = os.path.join(self.raw_dir, "two_stage_jobs", lookup_id)
        cpcm_file = f"{base}.tzvpd.cpcm"
        cpcm_corr_file = f"{base}.tzvpd.cpcm_corr"

        with open(out_path, "w") as f:

            # HEADER
            f.write(f"{lookup_id} : DFT_CPCM_BP86_def2-TZVP+def2-TZVPD_SP\n")

            # ENERGY
            f.write("\n##################################################\n#ENERGY\n")
            f.write(f"FINAL SINGLE POINT ENERGY      {data['energy']}\n")

            # DIPOLE
            f.write("\n##################################################\n#DIPOLE MOMENT (Debye)\n")
            f.write(" ".join(str(x) for x in data["dipole"]) + "\n")

            # XYZ
            f.write("\n##################################################\n#XYZ_FILE\n")
            f.write(f"{len(atoms)}\n")
            f.write(f"Coordinates from ORCA-job geo_opt_tzvp E {data['energy']}\n")
            for el, x, y, z in atoms:
                f.write(f"  {el:<2}  {x:20.14f}  {y:20.14f}  {z:20.14f}\n")

            # COSMO (CPCM)
            f.write("\n##################################################\n#COSMO\n")
            if os.path.exists(cpcm_file):
                with open(cpcm_file) as r:
                    f.write(r.read())
            else:
                f.write("# Missing CPCM file\n")

            # COSMO_corrected (CPCM_CORR)
            f.write("\n##################################################\n#COSMO_corrected\n")
            if os.path.exists(cpcm_corr_file):
                with open(cpcm_corr_file) as r:
                    f.write(r.read())
            else:
                f.write("# Missing CPCM_CORR file\n")

            # ADJACENCY
            f.write("\n##################################################\n#ADJACENCY_MATRIX\n")
            for row in data["adjacency_matrix"]:
                f.write("".join(f"{int(v):4d}" for v in row) + "\n")

            # POLARIZABILITIES (empty for now)
            f.write("\n##################################################\n#ATOMIC_POLARIZABILITIES\n")
            for line in data["atomic_polarizabilities"]:
                f.write(line + "\n")

    # ------------------------------------------------------------
    # XYZ loader
    # ------------------------------------------------------------
    def _load_xyz(self, path):
        atoms = []
        with open(path) as f:
            lines = f.readlines()

        start = 2 if len(lines) > 2 and lines[0].strip().isdigit() else 0

        for line in lines[start:]:
            parts = line.split()
            if len(parts) == 4:
                atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))

        return atoms

    # ------------------------------------------------------------
    # Summary writer
    # ------------------------------------------------------------
    def _write_summary(self):
        summary_path = os.path.join(self.outputs_dir, "orcacosmo_summary.json")

        with AtomicWriter(summary_path) as f:
            json.dump(self.successful_outputs, f, indent=2)

        self.log(f"Wrote summary with {len(self.successful_outputs)} entries: {summary_path}")
