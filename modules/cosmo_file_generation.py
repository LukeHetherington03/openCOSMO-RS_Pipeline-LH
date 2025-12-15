import os
import glob
import subprocess
import pandas as pd
import re
import csv
from rdkit import Chem
from .molecule_utils import MoleculeUtils

class CosmoFileGenerator:
    def __init__(self, cosmo_root="pipeline_data/6_cosmo_files"):
        self.cosmo_root = cosmo_root
        os.makedirs(self.cosmo_root, exist_ok=True)

    def _run_orca(self, inp_file, cwd=None):
        """
        Run ORCA on a given input file.
        By default, executes in pipeline_data/tmp_exec so scratch files
        are isolated and ignored by Git.

        Parameters
        ----------
        inp_file : str
            Path to the ORCA input file (.inp).
        cwd : str, optional
            Working directory for the ORCA process. If None, defaults to
            /home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/pipeline_data/tmp_exec.

        Returns
        -------
        success : bool
            True if ORCA terminated normally, False otherwise.

        Raises
        ------
        FileNotFoundError
            If the cwd provided does not exist or is not a directory.
        """

        # Default to tmp_exec if no cwd provided
        if cwd is None:
            cwd = os.path.join(
                "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/pipeline_data", "tmp_exec"
            )
            os.makedirs(cwd, exist_ok=True)

        # Validate cwd
        if not os.path.exists(cwd):
            raise FileNotFoundError(f"Working directory does not exist: {cwd}")
        if not os.path.isdir(cwd):
            raise FileNotFoundError(f"Working directory is not a directory: {cwd}")

        out_file = os.path.splitext(inp_file)[0] + ".out"
        try:
            with open(out_file, "w") as fout:
                subprocess.run(
                    ["orca", inp_file],
                    cwd=cwd,
                    stdout=fout,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=True
                )
            return True
        except subprocess.CalledProcessError as e:
            print(f"ORCA job failed for {os.path.basename(inp_file)} (returncode={e.returncode}). "
                  f"See {out_file} for details.")
            return False

    def _build_orcacosmo_file(self,
                              structname,
                              method,
                              filename_final_log,
                              filename_final_xyz,
                              filename_final_cpcm,
                              filename_final_cpcm_corr=None,
                              mol=None):
        """
        Build a consolidated .orcacosmo file from ORCA outputs.

        Notes
        -----
        This method was previously named `_concatenate_output`. It has been renamed
        to `_build_orcacosmo_file` to better reflect its purpose: assembling and
        structuring the outputs from a single conformer run into one .orcacosmo file.
        """
        out_path = f"{structname}.orcacosmo"

        with open(out_path, 'w') as file:
            file.write(f"{structname} : {method}\n")

            # ---------------- ENERGY ----------------
            file.write('\n' + '#' * 50 + '\n')
            file.write('#ENERGY\n')
            line_final_energy = ''
            dipole_moment = None

            with open(filename_final_log, 'r') as log_file:
                terminated_normally = False
                did_not_converge = False
                for line in log_file:
                    if '****ORCA TERMINATED NORMALLY****' in line:
                        terminated_normally = True
                    if 'The optimization did not converge but reached' in line:
                        did_not_converge = True
                    re_match = re.match(r'.*FINAL\s+SINGLE\s+POINT\s+ENERGY.+', line)
                    if re_match:
                        line_final_energy = line
                    if line.strip().startswith('x,y,z [Debye]:'):
                        dipole_moment = ' '.join(line.strip().split()[-3:])
                if not terminated_normally or did_not_converge:
                    raise RuntimeError('The calculation did not converge or did not terminate normally.')

            file.write(line_final_energy)

            if dipole_moment:
                file.write('\n' + '#' * 50 + '\n')
                file.write('#DIPOLE MOMENT (Debye)\n')
                file.write(f'{dipole_moment}\n')

            # ---------------- XYZ GEOMETRY ----------------
            file.write('\n' + '#' * 50 + '\n')
            file.write('#XYZ_FILE\n')
            with open(filename_final_xyz, 'r') as xyz_file:
                for line in xyz_file:
                    file.write(line)

            # ---------------- COSMO DATA ----------------
            if filename_final_cpcm is not None:
                file.write('\n' + '#' * 50 + '\n')
                file.write('#COSMO\n')
                with open(filename_final_cpcm, 'r') as cpcm_file:
                    for line in cpcm_file:
                        file.write(line)

            # ---------------- CORRECTED COSMO DATA ----------------
            if filename_final_cpcm_corr is not None and os.path.exists(filename_final_cpcm_corr):
                file.write('\n' + '#' * 50 + '\n')
                file.write('#COSMO_corrected\n')
                with open(filename_final_cpcm_corr, 'r') as cpcm_file:
                    for line in cpcm_file:
                        file.write(line)

            # ---------------- ADJACENCY MATRIX ----------------
            if mol:
                file.write('\n' + '#' * 50 + '\n')
                file.write('#ADJACENCY_MATRIX\n')
                adjacency_matrix = Chem.rdmolops.GetAdjacencyMatrix(mol, useBO=True)
                for i_row in range(adjacency_matrix.shape[0]):
                    line = ''.join(['{:4d}'.format(int(v)) for v in adjacency_matrix[i_row, :]]) + '\n'
                    file.write(line)

        return out_path

    def _build_output_dir(self, dataset, generation_engine, pruning_method, optimisation_engine):
        """Build output directory path based on metadata."""
        root = os.path.join(self.cosmo_root, dataset,
                            generation_engine, pruning_method,
                            optimisation_engine)
        os.makedirs(root, exist_ok=True)
        return root

    def run_cosmo_for_xyz(self,
                        xyz_file,
                        method="B3LYP",
                        basis="def2-SVP",
                        solvent="Water",
                        charge=0,
                        multiplicity=1):
        """
        Wrapper: take a single XYZ file, run ORCA in tmp_exec, and build the .orcacosmo file.
        Only .orcacosmo is stored in cosmo_root; .inp/.out/.cpcm stay in tmp_exec.
        """
        lookup_id = os.path.splitext(os.path.basename(xyz_file))[0]
        inchi = lookup_id.split("_conf")[0]

        mol_dir = os.path.join(self.cosmo_root, inchi)
        os.makedirs(mol_dir, exist_ok=True)

        tmp_exec = os.path.join(os.path.dirname(self.cosmo_root), "tmp_exec")
        os.makedirs(tmp_exec, exist_ok=True)


        # Write .inp into tmp_exec
        inp_file = os.path.join(tmp_exec, f"{lookup_id}.inp")
        MoleculeUtils.xyz_to_orca_inp(
            xyz_file, inp_file,
            method=method,
            basis=basis,
            solvent=solvent,
            charge=charge,
            multiplicity=multiplicity
        )

        # Run ORCA in tmp_exec
        ok = self._run_orca(inp_file, cwd=tmp_exec)
        if not ok:
            raise RuntimeError(f"ORCA job failed for {lookup_id}")

        # Build .orcacosmo in mol_dir using tmp_exec outputs
        log_file = os.path.splitext(inp_file)[0] + ".out"
        cpcm_file = os.path.splitext(inp_file)[0] + ".cpcm"
        cpcm_corr_file = os.path.splitext(inp_file)[0] + ".cpcm_corr"

        orcacosmo_path = self._build_orcacosmo_file(
            structname=os.path.join(mol_dir, lookup_id),
            method=f"{method}_{basis}_CPCM({solvent})",
            filename_final_log=log_file,
            filename_final_xyz=xyz_file,
            filename_final_cpcm=cpcm_file,
            filename_final_cpcm_corr=cpcm_corr_file,
            mol=None
        )

        return orcacosmo_path


