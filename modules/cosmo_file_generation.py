import os
import glob
import subprocess
import pandas as pd
import re
from rdkit import Chem
from .molecule_utils import MoleculeUtils

class CosmoFileGenerator:
    def __init__(self, cosmo_root="pipeline_data/6_cosmo_files"):
        self.cosmo_root = cosmo_root
        os.makedirs(self.cosmo_root, exist_ok=True)

    def _run_orca(self, inp_file, cwd=None):
        """Run ORCA on a given input file, write output to .out, suppress terminal spam."""
        out_file = os.path.splitext(inp_file)[0] + ".out"
        try:
            with open(out_file, "w") as fout:
                # Redirect both stdout and stderr into the .out file
                result = subprocess.run(
                    ["orca", inp_file],
                    cwd=cwd,
                    stdout=fout,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=True
                )
            return True
        except subprocess.CalledProcessError as e:
            # Only print a concise error message, not the full stdout/stderr
            print(f"ORCA job failed for {os.path.basename(inp_file)} (returncode={e.returncode}). "
                f"See {out_file} for details.")
            return False


    def _concatenate_output(self, structname, method,
                            filename_final_log,
                            filename_final_xyz,
                            filename_final_cpcm,
                            filename_final_cpcm_corr=None,
                            mol=None):
        """Build a .orcacosmo file from ORCA outputs."""
        out_path = f"{structname}.orcacosmo"
        with open(out_path, 'w') as file:
            file.write(f"{structname} : {method}\n")

            file.write('\n'+'#'*50+'\n')
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
                file.write('\n'+'#'*50+'\n')
                file.write('#DIPOLE MOMENT (Debye)\n')
                file.write(f'{dipole_moment}\n')

            file.write('\n'+'#'*50+'\n')
            file.write('#XYZ_FILE\n')
            with open(filename_final_xyz, 'r') as xyz_file:
                for line in xyz_file:
                    file.write(line)

            if filename_final_cpcm is not None:
                file.write('\n'+'#'*50+'\n')
                file.write('#COSMO\n')
                with open(filename_final_cpcm, 'r') as cpcm_file:
                    for line in cpcm_file:
                        file.write(line)

            if filename_final_cpcm_corr is not None and os.path.exists(filename_final_cpcm_corr):
                file.write('\n'+'#'*50+'\n')
                file.write('#COSMO_corrected\n')
                with open(filename_final_cpcm_corr, 'r') as cpcm_file:
                    for line in cpcm_file:
                        file.write(line)

            if mol:
                file.write('\n'+'#'*50+'\n')
                file.write('#ADJACENCY_MATRIX\n')
                adjacency_matrix = Chem.rdmolops.GetAdjacencyMatrix(mol, useBO=True)
                for i_row in range(adjacency_matrix.shape[0]):
                    line = ''.join(['{:4d}'.format(int(v)) for v in adjacency_matrix[i_row, :]]) + '\n'
                    file.write(line)
        return out_path

    def _build_output_dir(self, df):
        """Build output directory path based on dataset metadata."""
        dataset = df["dataset"].iloc[0]
        generation_engine = df["generation_engine"].iloc[0]
        pruning_method = df["pruning_method"].iloc[0]
        optimisation_engine = df["optimisation_engine"].iloc[0]

        root = os.path.join(self.cosmo_root, dataset,
                            generation_engine, pruning_method,
                            optimisation_engine)
        os.makedirs(root, exist_ok=True)
        return root

    def generate_orca_cosmo_from_folder(self, optimisation_dir,
                                        method="B3LYP",
                                        basis="def2-SVP",
                                        solvent="Water",
                                        charge=0,
                                        multiplicity=1):
        """
        Generate ORCA CPCM jobs for all conformers in a folder.
        Groups conformers by InChIKey and writes .orcacosmo files.
        Returns dict { inchi_key: { 'dir': path, 'n_confs': int } }
        """
        summary_csv = os.path.join(optimisation_dir, "_optimisation_summary.csv")
        if not os.path.exists(summary_csv):
            raise FileNotFoundError(f"Missing summary file: {summary_csv}")

        df = pd.read_csv(summary_csv)
        cosmo_root_for_run = self._build_output_dir(df)

        df["inchi_key"] = df["lookup_id"].str.split("_conf").str[0]

        results = {}
        for inchi, group in df.groupby("inchi_key"):
            mol_dir = os.path.join(cosmo_root_for_run, inchi)
            os.makedirs(mol_dir, exist_ok=True)

            conf_count = 0
            for _, row in group.iterrows():
                lookup_id = row["lookup_id"]
                xyz_file = row["xyz_file"]

                inp_file = os.path.join(mol_dir, f"{lookup_id}.inp")
                MoleculeUtils.xyz_to_orca_inp(
                    xyz_file, inp_file,
                    method=method,
                    basis=basis,
                    solvent=solvent,
                    charge=charge,
                    multiplicity=multiplicity
                )

                ok = self._run_orca(inp_file, cwd=mol_dir)
                if not ok:
                    print(f"ORCA job failed for {lookup_id}. See diagnostic output above.")
                    continue

                # Build .orcacosmo file
                log_file = os.path.splitext(inp_file)[0] + ".out"
                cpcm_file = os.path.splitext(inp_file)[0] + ".cpcm"
                cpcm_corr_file = os.path.splitext(inp_file)[0] + ".cpcm_corr"
                orcacosmo_path = self._concatenate_output(
                    structname=os.path.join(mol_dir, lookup_id),
                    method=f"{method}_{basis}_CPCM({solvent})",
                    filename_final_log=log_file,
                    filename_final_xyz=xyz_file,
                    filename_final_cpcm=cpcm_file,
                    filename_final_cpcm_corr=cpcm_corr_file,
                    mol=None  # optionally pass RDKit mol if available
                )
                print(f"Generated {orcacosmo_path}")
                conf_count += 1

            results[inchi] = {"dir": mol_dir, "n_confs": conf_count}

        return results
