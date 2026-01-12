#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
from dataclasses import dataclass
from modules.utils import molecule_utils


# ============================================================
#  Result Dataclass
# ============================================================

@dataclass
class OrcaCosmoBuildResult:
    lookup_id: str
    input_xyz: str
    input_json: str
    output_orcacosmo: str


# ============================================================
#  Exception
# ============================================================

class OrcaCosmoBuilderError(Exception):
    pass


# ============================================================
#  ORCA COSMO Builder Stage
# ============================================================

class OrcaCosmoStage:
    """
    Builds ORCA-compatible .orcacosmo files from:
      - <lookup_id>.xyz
      - <lookup_id>_cosmo_surface.json

    Job model:
      - job.items = ["ABC_conf0", "ABC_conf1", ...]
      - job.pending_items = same list initially
      - job.completed_items = []

    This stage:
      - Resolves files dynamically from job.parameters["input_folder"]
      - Extracts InChIKey from JSON
      - Uses molecule_utils for adjacency
      - Skips missing pairs
      - Writes .orcacosmo into output_parsed/
    """

    def __init__(self, job):
        from modules.provenance.job_manager import Job
        if not isinstance(job, Job):
            raise TypeError("OrcaCosmoBuilder must be instantiated with a Job object")

        self.job = job
        self.xyz_suffix = ".xyz"
        self.cosmo_suffix = "_cosmo_surface.json"

    # ------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------
    def run(self):
        job = self.job

        try:
            job.write_pipeline_log("Starting ORCA COSMO builder stage")

            input_folder = job.parameters["input_folder"]
            results = []

            # --------------------------------------------------------
            # Loop over pending lookup_ids
            # --------------------------------------------------------
            for lookup_id in job.pending_items.copy():

                xyz_path = os.path.join(input_folder, f"{lookup_id}{self.xyz_suffix}")
                json_path = os.path.join(input_folder, f"{lookup_id}{self.cosmo_suffix}")

                # Skip missing pairs
                if not os.path.exists(xyz_path) or not os.path.exists(json_path):
                    job.write_pipeline_log(
                        f"Skipping {lookup_id}: missing required files "
                        f"(xyz={os.path.exists(xyz_path)}, json={os.path.exists(json_path)})"
                    )
                    job.update_progress(lookup_id)
                    continue

                job.write_pipeline_log(f"Building ORCA COSMO file for {lookup_id}")

                result = self._build_single(lookup_id, xyz_path, json_path)
                results.append(result)

                job.update_progress(lookup_id)

            job.mark_complete()
            job.write_pipeline_log("ORCA COSMO builder completed successfully.")

            return results

        except Exception as e:
            job.fail(str(e))
            raise

    # ------------------------------------------------------------
    # Build a single ORCA COSMO file
    # ------------------------------------------------------------
    def _build_single(self, lookup_id: str, xyz_path: str, json_path: str) -> OrcaCosmoBuildResult:
        job = self.job

        # Load COSMO JSON
        with open(json_path) as f:
            cosmo_data = json.load(f)

        # Extract InChIKey from JSON (correct)
        inchi_key = cosmo_data["inchi_key"]

        # Load molecule info (RDKit mol + adjacency)
        mol_info = molecule_utils.get_molecule_info_from_inchi_key(inchi_key)
        adjacency = mol_info["adjacency"]

        # Load XYZ atoms
        xyz_atoms = self._load_xyz(xyz_path)

        # Prepare output path
        output_dir = os.path.join(job.output_dir, "output_parsed")
        os.makedirs(output_dir, exist_ok=True)

        out_path = os.path.join(output_dir, f"{lookup_id}.orcacosmo")

        # Write ORCA COSMO file
        with open(out_path, "w") as f:
            self._write_header(f, lookup_id, cosmo_data)
            self._write_energy_section(f, cosmo_data)
            self._write_dipole_section(f, cosmo_data)
            self._write_xyz_section(f, xyz_atoms)
            self._write_cosmo_section(f, cosmo_data)
            self._write_surface_points(f, cosmo_data)
            self._write_cosmo_corrected(f, cosmo_data)
            self._write_adjacency(f, adjacency)
            self._write_polarizabilities(f, cosmo_data)

        return OrcaCosmoBuildResult(
            lookup_id=lookup_id,
            input_xyz=xyz_path,
            input_json=json_path,
            output_orcacosmo=out_path,
        )

    # ============================================================
    # Section writers
    # ============================================================

    def _write_header(self, f, lookup_id: str, data: dict):
        method_str = data.get("method", "UNKNOWN_METHOD")
        f.write(f"{lookup_id} : {method_str}\n\n")

    def _write_energy_section(self, f, data: dict):
        f.write("##################################################\n")
        f.write("#ENERGY\n")
        energy = data.get("total_energy", 0.0)
        f.write(f"FINAL SINGLE POINT ENERGY      {energy:20.12f}\n\n")

    def _write_dipole_section(self, f, data: dict):
        f.write("##################################################\n")
        f.write("#DIPOLE MOMENT (Debye)\n")
        dx = dy = dz = 0.0
        f.write(f"{dx:.6f} {dy:.6f} {dz:.6f}\n\n")

    def _write_xyz_section(self, f, xyz_atoms):
        f.write("##################################################\n")
        f.write("#XYZ_FILE\n")
        f.write(f"{len(xyz_atoms)}\n")
        f.write("Coordinates from pipeline XYZ\n")

        for el, x, y, z in xyz_atoms:
            f.write(f"  {el:<2}  {x:20.14f} {y:20.14f} {z:20.14f}\n")

        f.write("\n")

    def _write_cosmo_section(self, f, data: dict):
        f.write("##################################################\n")
        f.write("#COSMO\n")
        f.write("\n")

    def _write_surface_points(self, f, data: dict):
        f.write("#------------------------------------------------------------\n")
        f.write("# SURFACE POINTS (A.U.)    (Hint - charge NOT scaled by FEps)\n")
        f.write("#------------------------------------------------------------\n")
        f.write("\n")

    def _write_cosmo_corrected(self, f, data: dict):
        f.write("##################################################\n")
        f.write("#COSMO_corrected\n")
        f.write("\n")

    def _write_adjacency(self, f, adjacency):
        f.write("##################################################\n")
        f.write("#ADJACENCY_MATRIX\n")
        for row in adjacency:
            f.write("".join(f"{v:4d}" for v in row) + "\n")
        f.write("\n")

    def _write_polarizabilities(self, f, data: dict):
        if "polarizabilities" in data:
            f.write("##################################################\n")
            f.write("#ATOMIC POLARIZABILITIES\n")
            for line in data["polarizabilities"]:
                f.write(line)
            f.write("\n")

    # ============================================================
    # Helpers
    # ============================================================

    def _load_xyz(self, path: str):
        atoms = []
        with open(path) as f:
            lines = f.readlines()

        # Skip natoms + comment if present
        start_idx = 0
        if len(lines) >= 2:
            try:
                int(lines[0].strip())
                start_idx = 2
            except ValueError:
                pass

        for line in lines[start_idx:]:
            parts = line.split()
            if len(parts) != 4:
                continue
            el = parts[0]
            x, y, z = map(float, parts[1:])
            atoms.append((el, x, y, z))

        return atoms
