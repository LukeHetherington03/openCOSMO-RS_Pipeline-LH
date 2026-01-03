import os          # for filesystem paths and directory creation
import subprocess  # for running external programs (xtb, gxtb, orca)
import pandas as pd  # for reading/writing CSVs (lookup and master tables)
import shutil 
import time 

# At module level (conformer_optimisation.py)
OPTIMISATION_METHODS = {
    "xtb": {"runner": "_run_xtb"},
    "gxtb": {"runner": "_run_gxtb"},
    "forcefield": {"runner": "_run_forcefield"},
    "dft": {"runner": "_run_dft", "mode": "standard"},
    "dft_fast": {"runner": "_run_dft", "mode": "fast"},
    "dft_final": {"runner": "_run_dft", "mode": "final"},
    "dft_cpcm": {"runner": "_run_dft", "mode": "cpcm"},
}


class ConformerOptimiser:
    def __init__(self,
                 input_xyz_dir,
                 lookup_csv,
                 output_dir,
                 tmp_exec,
                 engine="gxtb",
                 mute=False,
                 master_csv=None,
                 **params):
        """
        Parameters
        ----------
        input_xyz_dir : str
            Directory containing the input .xyz conformer files.
        lookup_csv : str
            Path to the lookup CSV produced by pruning stage.
            This file contains the conformer IDs to be optimised.
        output_dir : str
            Directory where optimisation results will be written.
            Subfolders 'xyz/' and 'log/' will be created inside this directory.
        tmp_exec : str
            Scratch directory for temporary execution files.
            External binaries (xtb, gxtb, orca) will run here.
        mute : Bool
            Job counter verbosity
        engine : str, default="gxtb"
            Optimisation backend to use. Supported values:
                - "xtb"  : classic GFN-xTB program
                - "gxtb" : new g-xtb binary
                - "orca" : ORCA DFT engine
        master_csv : str, optional
            Path to pipeline-specific master conformer table.
            This will be updated with optimisation energies/status.
        params : dict
            Engine-specific parameters, e.g.:
                - level="GFN2-xTB" for xtb
                - functional="B3LYP", basis="def2-SVP" for orca
                - max_iter=250 for iteration limits
        """
        self.input_xyz_dir = input_xyz_dir
        self.lookup_csv = lookup_csv
        self.output_dir = output_dir
        self.tmp_exec = tmp_exec
        self.engine = engine.lower()
        self.mute = mute
        self.master_csv = master_csv
        self.params = params

        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.tmp_exec, exist_ok=True)

    def optimise(self):
        df_lookup = pd.read_csv(self.lookup_csv, comment="#")
        df_master = pd.read_csv(self.master_csv)
        total_jobs = len(df_lookup)

        results = []

        for idx, row in enumerate(df_lookup.itertuples(index=False), start=1):

            if not self.mute:
                self._print_progress(idx, total_jobs, row.lookup_id)

            xyz_file = self._find_xyz_file(row.lookup_id)
            if not xyz_file:
                if not self.mute:
                    print(f"Missing XYZ for {row.lookup_id}, skipping.")
                continue

            start = time.perf_counter()
            result = self._run_engine(xyz_file)
            elapsed = time.perf_counter() - start


            # Update master CSV
            mask = df_master["lookup_id"] == row.lookup_id
            for col, val in result["updates"].items():
                df_master.loc[mask, col] = val

            results.append({
                "lookup_id": row.lookup_id,
                "energy": result["updates"].get("energy_xtb")
                        or result["updates"].get("energy_gxtb")
                        or result["updates"].get("energy_dft"),
                "status": result["updates"].get("status_xtb")
                        or result["updates"].get("status_gxtb")
                        or result["updates"].get("status_dft"),
                "xyz_file": result["out_file"],
                "log_file": result["log_file"],
                "elapsed_seconds": elapsed
            })

        df_master.to_csv(self.master_csv, index=False)

        summary_path = os.path.join(self.output_dir, "_optimisation_summary.csv")
        pd.DataFrame(results).to_csv(summary_path, index=False)
        return summary_path

    def _run_engine(self, xyz_file):
        """Dispatch to the correct engine runner. Only takes an .xyz file."""
        if self.engine not in OPTIMISATION_METHODS:
            raise ValueError(f"Unknown optimisation engine: {self.engine}")

        method = OPTIMISATION_METHODS[self.engine]
        runner = getattr(self, method["runner"])
        if "mode" in method:
            return runner(xyz_file, mode=method["mode"])
        else:
            return runner(xyz_file)


    def _run_gxtb(self, xyz_file, quiet=False, ncores=None):
        import os, shutil, subprocess, math

        base = os.path.splitext(os.path.basename(xyz_file))[0]

        # Clean tmp_exec scratch before run
        tmp_xyz = os.path.join(self.tmp_exec, f"{base}.xyz")
        if os.path.exists(tmp_xyz):
            os.remove(tmp_xyz)
        opt_tmp = os.path.join(self.tmp_exec, "xtbopt.xyz")
        if os.path.exists(opt_tmp):
            os.remove(opt_tmp)

        # Copy input xyz into tmp_exec
        shutil.copy(xyz_file, tmp_xyz)

        # Prepare output and log paths
        log_file = os.path.join(self.output_dir, "log", f"{base}_gxtb.log")
        out_file = os.path.join(self.output_dir, "xyz", f"{base}_opt.xyz")
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        os.makedirs(os.path.dirname(log_file), exist_ok=True)

        # Decide number of cores
        max_cores = os.cpu_count() or 1
        if ncores is None:
            ncores = max(1, math.floor(max_cores * 0.8))  # default: 80% of max cores

        # Environment with OMP_NUM_THREADS
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(ncores)

        # Run xtb with gxtb as driver
        cmd = [
            "xtb", tmp_xyz,
            "--driver", "gxtb -grad -c xtbdriver.xyz",
            "--opt",
            "--iterations", str(self.params.get("max_iter", 250))
        ]

        energy, status = None, "unknown"
        try:
            if quiet:
                subprocess.run(cmd, cwd=self.tmp_exec,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT,
                            check=True, env=env)
            else:
                with open(log_file, "w") as log:
                    subprocess.run(cmd, cwd=self.tmp_exec,
                                stdout=log, stderr=log,
                                check=True, env=env)

            # xtb writes xtbopt.xyz in cwd
            if os.path.exists(opt_tmp):
                shutil.move(opt_tmp, out_file)
                status = "converged"

            # Parse energy from log
            if not quiet and os.path.exists(log_file):
                with open(log_file) as log:
                    for line in log:
                        if "TOTAL ENERGY" in line.upper() or line.strip().lower().startswith("total"):
                            try:
                                energy = float(line.split()[-1])
                                status = "converged"
                            except Exception:
                                pass
                            break
        except subprocess.CalledProcessError:
            status = "failed"

        # Clean tmp_exec scratch
        for fn in [tmp_xyz]:
            if os.path.exists(fn):
                os.remove(fn)

        return {
            "out_file": out_file if os.path.exists(out_file) else None,
            "log_file": None if quiet else log_file,
            "updates": {
                "energy_gxtb": energy,
                "status_gxtb": status,
                "ncores_used": ncores
            }
        }




    def _run_xtb(self, xyz_file):
        base = os.path.splitext(os.path.basename(xyz_file))[0]
        tmp_xyz = os.path.join(self.tmp_exec, f"{base}.xyz")
        shutil.copy(xyz_file, tmp_xyz)

        log_file = os.path.join(self.tmp_exec, f"{base}_xtb.log")
        out_file = os.path.join(self.output_dir, "xyz", f"{base}_opt.xyz")
        os.makedirs(os.path.dirname(out_file), exist_ok=True)

        cmd = ["xtb", tmp_xyz, "--opt", "--gfn", self.params.get("level","2"),
            "--iterations", str(self.params.get("max_iter",250))]

        try:
            with open(log_file, "w") as log:
                subprocess.run(cmd, cwd=self.tmp_exec, stdout=log, stderr=log, check=True)

            opt_tmp = os.path.join(self.tmp_exec, "xtbopt.xyz")
            if os.path.exists(opt_tmp):
                shutil.move(opt_tmp, out_file)

            energy, status = None, "unknown"
            with open(log_file) as log:
                for line in log:
                    if "TOTAL ENERGY" in line.upper():
                        try:
                            energy = float(line.split()[-2])
                            status = "converged"
                        except Exception:
                            pass
        except subprocess.CalledProcessError:
            energy, status = None, "failed"

        return {
            "out_file": out_file,
            "log_file": log_file,
            "updates": {
                "energy_xtb": energy,
                "status_xtb": status
            }
        }
    
    def _generate_dft_input(self, xyz_file, mode="standard"):
        """Generate an ORCA input file for DFT depending on mode and leave it in tmp_exec for inspection."""
        import os, shutil

        # Defaults
        functional = self.params.get("functional", "B3LYP")
        basis = self.params.get("basis", "def2-SVP")
        max_iter = self.params.get("max_iter", 250)
        charge = self.params.get("charge", 0)
        multiplicity = self.params.get("multiplicity", 1)

        base = os.path.splitext(os.path.basename(xyz_file))[0]
        tmp_xyz = os.path.join(self.tmp_exec, f"{base}.xyz")
        shutil.copy(xyz_file, tmp_xyz)

        inp_file = os.path.join(self.tmp_exec, f"{base}.inp")
        lines = []

        # Mode branching
        if mode == "fast":
            lines.append("! BP86 def2-SVP Opt")
        elif mode == "final":
            lines.append("! BP86 def2-TZVP Opt")
        elif mode == "cpcm":
            lines.append("! BP86 def2-TZVP Opt CPCM")
            lines.append("%cpcm")
            if "cpcm_radii" in self.params:
                for atom_num, radius in self.params["cpcm_radii"].items():
                    lines.append(f"  radius[{atom_num}] {radius}")
            if "cut_area" in self.params:
                lines.append(f"  cut_area {self.params['cut_area']}")
            lines.append("end")
        else:
            lines.append(f"! {functional} {basis} Opt")

        # Geometry optimisation block
        lines.append("%geom")
        lines.append(f"  MaxIter {max_iter}")
        lines.append("end")

        # Output verbosity block (correct syntax)
        lines.append("%output")
        lines.append("  PrintLevel Mini")
        lines.append("end")

        # Coordinates
        lines.append(f"* xyz {charge} {multiplicity}")
        with open(tmp_xyz) as xyz:
            lines.extend(xyz.readlines()[2:])  # skip atom count + comment
        lines.append("*")

        # Write input file
        with open(inp_file, "w") as f:
            f.write("\n".join(lines))

        # Debug print so you know where to look
        print(f"ORCA input written to: {inp_file}")

        return inp_file, base


    def _run_dft(self, xyz_file, mode="standard", quiet=False):
        import os, subprocess, shutil

        inp_file, base = self._generate_dft_input(xyz_file, mode=mode)

        # Use mode in filenames
        log_file = os.path.join(self.output_dir, "log", f"{base}_{mode}.log")
        out_file = os.path.join(self.output_dir, "xyz", f"{base}_{mode}_opt.xyz")
        prop_file = os.path.join(self.output_dir, "log", f"{base}_{mode}.property.txt")

        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        os.makedirs(os.path.dirname(log_file), exist_ok=True)

        energy, status = None, "unknown"
        try:
            with open(log_file, "w") as log:
                subprocess.run(["orca", inp_file],
                            cwd=self.tmp_exec,
                            stdout=log,
                            stderr=log,
                            check=True)

            # Copy final geometry
            final_xyz = os.path.join(self.tmp_exec, f"{base}.xyz")
            if os.path.exists(final_xyz):
                shutil.copy(final_xyz, out_file)
                status = "converged"

            # Copy property file if present
            tmp_prop = os.path.join(self.tmp_exec, f"{base}.property.txt")
            if os.path.exists(tmp_prop):
                shutil.copy(tmp_prop, prop_file)

            # Parse energy from property file if available
            if os.path.exists(prop_file):
                with open(prop_file) as pf:
                    for line in pf:
                        if line.strip().startswith("FINAL SINGLE POINT ENERGY") or line.strip().startswith("Total Energy"):
                            try:
                                energy = float(line.split()[-1])
                            except Exception:
                                pass
                            break
        except subprocess.CalledProcessError:
            status = "failed"

        return {
            "out_file": out_file if os.path.exists(out_file) else None,
            "log_file": log_file if os.path.exists(log_file) else None,
            "property_file": prop_file if os.path.exists(prop_file) else None,
            "updates": {
                "energy_dft": energy,
                "status_dft": status,
                "mode_dft": mode,
                "functional_dft": self.params.get("functional"),
                "basis_dft": self.params.get("basis"),
                "charge_dft": self.params.get("charge"),
                "multiplicity_dft": self.params.get("multiplicity"),
                "solvent_dft": self.params.get("solvent"),
            }
        }


    def _find_xyz_file(self, lookup_id):
        """
        Return the path to the XYZ file matching this lookup_id.
        Accepts either exact filename match or prefix match.
        """
        base = str(lookup_id)

        # Exact match
        exact = os.path.join(self.input_xyz_dir, f"{base}.xyz")
        if os.path.exists(exact):
            return exact

        # Prefix match (e.g. lookup_id = "mol_001" matches "mol_001_conf3.xyz")
        for fn in os.listdir(self.input_xyz_dir):
            if fn.startswith(base) and fn.endswith(".xyz"):
                return os.path.join(self.input_xyz_dir, fn)

        return None

    def _print_progress(self, idx, total, lookup_id):
        """Print a clean progress line showing job count and CPU allocation."""
        max_cores = os.cpu_count() or 1
        requested_cores = self.params.get("ncores")

        if requested_cores is None:
            requested_cores = max(1, int(max_cores * 0.8))

        cpu_pct = (requested_cores / max_cores) * 100

        print(
            f"[{idx:4d} / {total}]  "
            f"{self.engine:<6}  |  using {cpu_pct:4.0f}% CPUs "
            f"({requested_cores} cores)  |  lookup_id = {lookup_id}"
        )
