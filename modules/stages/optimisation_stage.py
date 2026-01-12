import os
import json
import time
import shutil
import subprocess
import pandas as pd

from modules.utils.atomic_write import AtomicWriter
from modules.utils.output_parsers import XtbCosmoParser, XtbTopologyParser


OPTIMISATION_METHODS = {
    "xtb": {"runner": "_run_xtb"},
    "gxtb": {"runner": "_run_gxtb"},
    "dft": {"runner": "_run_dft", "mode": "standard"},
    "dft_fast": {"runner": "_run_dft", "mode": "fast"},
    "dft_final": {"runner": "_run_dft", "mode": "final"},
    "dft_cpcm": {"runner": "_run_dft", "mode": "cpcm"},
}


class OptimisationRecord:
    def __init__(self, lookup_id, energy, xyz_path, log_path, method, status, metadata=None):
        self.lookup_id = lookup_id
        self.energy = energy
        self.xyz_path = xyz_path
        self.log_path = log_path
        self.method = method
        self.status = status
        self.metadata = metadata or {}

    def to_summary_row(self):
        return {
            "lookup_id": self.lookup_id,
            "energy": self.energy,
            "method": self.method,
            "status": self.status,
            "xyz_path": self.xyz_path,
            "log_path": self.log_path,
        }

    def to_energy_entry(self):
        return {
            "lookup_id": self.lookup_id,
            "energy": self.energy,
            "xyz_path": self.xyz_path,
            "log_path": self.log_path,
            "metadata": {
                "method": self.method,
                "status": self.status,
                **self.metadata
            }
        }


class OptimisationStage:

    # ------------------------------------------------------------
    # Entry point (now clean and readable)
    # ------------------------------------------------------------
    def run(self, job):
        params = job.parameters
        engine = params.get("engine", "gxtb").lower()
        max_iter = params.get("max_iter", 250)

        job.log_header("Starting optimisation stage")
        job.log(f"Engine: {engine}, max_iter={max_iter}")

        inputs, entries = self._prepare_inputs(job, params)
        molecules = self._group_by_molecule(entries)

        job.log(f"Input: {len(entries)} conformers across {len(molecules)} molecules")

        records = self._optimise_all(job, molecules, engine, max_iter)

        self._write_outputs(job, records, params, inputs["missing"])
        job.mark_complete()

    # ------------------------------------------------------------
    # Input preparation
    # ------------------------------------------------------------
    def _prepare_inputs(self, job, params):
        summary_file = params["summary_file"]

        inputs_dir = job.inputs_dir
        outputs_dir = job.outputs_dir

        inputs_xyz_dir = os.path.join(inputs_dir, "xyz")
        outputs_xyz_dir = os.path.join(outputs_dir, "xyz")
        outputs_log_dir = os.path.join(outputs_dir, "log")

        os.makedirs(inputs_xyz_dir, exist_ok=True)
        os.makedirs(outputs_xyz_dir, exist_ok=True)
        os.makedirs(outputs_log_dir, exist_ok=True)

        # Copy energies.json
        local_energies_path = os.path.join(inputs_dir, "energies.json")
        shutil.copy(summary_file, local_energies_path)

        with open(local_energies_path) as f:
            entries = json.load(f)

        # Copy XYZs
        missing = []
        for entry in entries:
            src = entry["xyz_path"]
            dst = os.path.join(inputs_xyz_dir, os.path.basename(src))

            if os.path.exists(src):
                shutil.copy(src, dst)
                entry["xyz_path"] = dst
            else:
                missing.append(entry["lookup_id"])

        if missing:
            job.log(f"[WARNING] Missing XYZs: {len(missing)}")
            for m in missing:
                job.log(f"  - {m}", indent=2)

        entries = [e for e in entries if e["lookup_id"] not in missing]

        if not entries:
            raise ValueError("No valid conformers available for optimisation.")

        return {"inputs_dir": inputs_dir, "outputs_dir": outputs_dir, "missing": missing}, entries

    # ------------------------------------------------------------
    # Group by molecule
    # ------------------------------------------------------------
    def _group_by_molecule(self, entries):
        molecules = {}
        for e in entries:
            mol = e["lookup_id"].split("_")[0]
            molecules.setdefault(mol, []).append(e)
        return molecules

    # ------------------------------------------------------------
    # Optimise all molecules
    # ------------------------------------------------------------
    def _optimise_all(self, job, molecules, engine, max_iter):
        records = []
        global_start = time.perf_counter()
        total_confs = sum(len(v) for v in molecules.values())
        global_counter = 0

        for mol_idx, (mol_id, mol_entries) in enumerate(molecules.items(), start=1):
            mol_records, global_counter = self._optimise_molecule(
                job, mol_id, mol_entries, mol_idx, len(molecules),
                engine, max_iter, global_counter, total_confs
            )
            records.extend(mol_records)

        stage_elapsed = time.perf_counter() - global_start
        job.log_header("Optimisation complete")
        job.log(f"Total conformers attempted: {total_confs}")
        job.log(f"Total conformers optimised: {len(records)}")
        job.log(f"Stage time: {stage_elapsed:.2f} seconds")

        return records

    # ------------------------------------------------------------
    # Optimise a single molecule
    # ------------------------------------------------------------
    def _optimise_molecule(self, job, mol_id, mol_entries, mol_idx, total_mols,
                           engine, max_iter, global_counter, total_confs):

        job.log_section(f"Molecule {mol_idx}/{total_mols}: {mol_id}")

        mol_start = time.perf_counter()
        mol_times = []
        mol_records = []

        for entry in mol_entries:
            global_counter += 1
            record, elapsed = self._optimise_single(
                job, entry, engine, max_iter, global_counter, total_confs
            )
            mol_records.append(record)
            mol_times.append(elapsed)

        # Molecule summary
        mol_elapsed = time.perf_counter() - mol_start
        avg_time = sum(mol_times) / len(mol_times)

        job.log(f"Completed molecule {mol_id}:", indent=1)
        job.log(f"{len(mol_entries)} conformers processed", indent=2)
        job.log(f"Average time per conformer: {avg_time:.2f} seconds", indent=2)
        job.log(f"Total time: {mol_elapsed:.2f} seconds", indent=2)

        return mol_records, global_counter

    # ------------------------------------------------------------
    # Optimise a single conformer
    # ------------------------------------------------------------
    def _optimise_single(self, job, entry, engine, max_iter, global_counter, total_confs):
        lookup_id = entry["lookup_id"]
        xyz_path = entry["xyz_path"]

        job.log(f"[{global_counter}/{total_confs}] Optimising {lookup_id}", indent=2)

        start = time.perf_counter()
        result = self._run_engine(job, xyz_path, engine, max_iter)
        elapsed = time.perf_counter() - start

        job.log(f"Completed {lookup_id} in {elapsed:.2f} seconds ({result['status']})", indent=3)

        record = OptimisationRecord(
            lookup_id=lookup_id,
            energy=result["energy"],
            xyz_path=result["xyz_out"],
            log_path=result["log_out"],
            method=engine,
            status=result["status"],
            metadata={
                "elapsed_seconds": elapsed,
                "backend_meta": result.get("backend_meta", {})
            }
        )

        return record, elapsed

    # ------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------
    def _write_outputs(self, job, records, params, missing):
        outputs_dir = job.outputs_dir

        # summary.csv
        summary_path = os.path.join(outputs_dir, "summary.csv")
        rows = [r.to_summary_row() for r in records]
        pd.DataFrame(rows).to_csv(summary_path, index=False)

        # energies.json
        energies_path = os.path.join(outputs_dir, "energies.json")
        energy_entries = [r.to_energy_entry() for r in records]

        with AtomicWriter(energies_path) as f:
            json.dump(energy_entries, f, indent=2)

        # job_state.json
        job_state_path = os.path.join(job.job_dir, "job_state.json")
        with AtomicWriter(job_state_path) as f:
            json.dump(
                {
                    "stage": "optimisation",
                    "engine": params.get("engine"),
                    "max_iter": params.get("max_iter"),
                    "num_input": len(records) + len(missing),
                    "num_output": len(records),
                    "missing_xyz": missing,
                },
                f,
                indent=2,
            )

        job.log(f"Outputs written to: {outputs_dir}")

    # ------------------------------------------------------------
    # Engine dispatch (unchanged)
    # ------------------------------------------------------------
    def _run_engine(self, job, xyz_path, engine, max_iter):
        if engine not in OPTIMISATION_METHODS:
            raise ValueError(f"Unknown optimisation engine: {engine}")

        method = OPTIMISATION_METHODS[engine]
        runner = getattr(self, method["runner"])

        if "mode" in method:
            return runner(job, xyz_path, max_iter=max_iter, mode=method["mode"])
        else:
            return runner(job, xyz_path, max_iter=max_iter)

    # ------------------------------------------------------------
    # Backends (unchanged)
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # gXTB backend
    # ------------------------------------------------------------
    def _run_gxtb(self, job, xyz_path, max_iter, mode=None):
        """
        Proper gXTB optimisation using:
            xtb <xyz> --driver "gxtb -grad -c xtbdriver.xyz" --opt
        """

        import os, shutil, subprocess, math

        base = os.path.splitext(os.path.basename(xyz_path))[0]

        # ------------------------------------------------------------
        # tmp_exec scratch inside the job directory
        # ------------------------------------------------------------
        tmp_exec = os.path.join(job.job_dir, "tmp_exec")
        os.makedirs(tmp_exec, exist_ok=True)

        tmp_xyz = os.path.join(tmp_exec, f"{base}.xyz")
        opt_tmp = os.path.join(tmp_exec, "xtbopt.xyz")

        # Clean scratch
        for f in (tmp_xyz, opt_tmp):
            if os.path.exists(f):
                os.remove(f)

        shutil.copy(xyz_path, tmp_xyz)

        # ------------------------------------------------------------
        # Output paths
        # ------------------------------------------------------------
        log_out = os.path.join(job.outputs_dir, "log", f"{base}_gxtb.log")
        xyz_out = os.path.join(job.outputs_dir, "xyz", f"{base}_opt.xyz")

        os.makedirs(os.path.dirname(log_out), exist_ok=True)
        os.makedirs(os.path.dirname(xyz_out), exist_ok=True)

        # ------------------------------------------------------------
        # CPU allocation
        # ------------------------------------------------------------
        max_cores = os.cpu_count() or 1
        ncores = max(1, math.floor(max_cores * 0.8))

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(ncores)

        # ------------------------------------------------------------
        # Correct gXTB command
        # ------------------------------------------------------------
        cmd = [
            "xtb", tmp_xyz,
            "--driver", "gxtb -grad -c xtbdriver.xyz",
            "--opt",
            "--iterations", str(max_iter)
        ]

        energy = None
        status = "unknown"

        try:
            # Run gXTB
            with open(log_out, "w") as log:
                subprocess.run(
                    cmd,
                    cwd=tmp_exec,
                    stdout=log,
                    stderr=log,
                    check=True,
                    env=env
                )

            # Move optimised geometry
            if os.path.exists(opt_tmp):
                shutil.move(opt_tmp, xyz_out)
                status = "converged"

            # Parse energy from log
            with open(log_out) as log:
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

        # Clean scratch
        if os.path.exists(tmp_xyz):
            os.remove(tmp_xyz)

        return {
            "energy": energy,
            "status": status,
            "xyz_out": xyz_out if os.path.exists(xyz_out) else None,
            "log_out": log_out,
            "backend_meta": {
                "mode": mode,
                "ncores_used": ncores
            }
        }

    # ------------------------------------------------------------
    # Classic XTB backend
    # ------------------------------------------------------------
    def _run_xtb(self, job, xyz_path, max_iter):
        base = os.path.splitext(os.path.basename(xyz_path))[0]

        tmp_xyz = os.path.join(job.inputs_dir, "xyz", f"{base}_tmp.xyz")
        shutil.copy(xyz_path, tmp_xyz)

        log_out = os.path.join(job.outputs_dir, "log", f"{base}_xtb.log")
        xyz_out = os.path.join(job.outputs_dir, "xyz", f"{base}_opt.xyz")

        cmd = [
            "xtb", str(tmp_xyz),
            "--opt",
            "--gfn", job.parameters.get("level", "2"),
            "--iterations", str(max_iter)
        ]

        energy = None
        status = "unknown"

        try:
            with open(log_out, "w") as log:
                subprocess.run(cmd, stdout=log, stderr=log, check=True)

            opt_tmp = "xtbopt.xyz"
            if os.path.exists(opt_tmp):
                shutil.move(opt_tmp, xyz_out)
                status = "converged"

            with open(log_out) as log:
                for line in log:
                    if "TOTAL ENERGY" in line.upper():
                        try:
                            energy = float(line.split()[-2])
                        except Exception:
                            pass
                        break

        except subprocess.CalledProcessError:
            status = "failed"

        return {
            "energy": energy,
            "status": status,
            "xyz_out": xyz_out if os.path.exists(xyz_out) else None,
            "log_out": log_out,
            "backend_meta": {}
        }

    # ------------------------------------------------------------
    # ORCA DFT backend
    # ------------------------------------------------------------
    def _run_dft(self, job, xyz_path, max_iter, mode="standard"):
        base = os.path.splitext(os.path.basename(xyz_path))[0]

        tmp_xyz = os.path.join(job.inputs_dir, "xyz", f"{base}.xyz")
        shutil.copy(xyz_path, tmp_xyz)

        inp_file = os.path.join(job.inputs_dir, f"{base}.inp")
        log_out = os.path.join(job.outputs_dir, "log", f"{base}_{mode}.log")
        xyz_out = os.path.join(job.outputs_dir, "xyz", f"{base}_{mode}_opt.xyz")

        # Write ORCA input
        with open(inp_file, "w") as f:
            f.write(f"! B3LYP def2-SVP Opt\n")
            f.write(f"* xyz 0 1\n")
            with open(tmp_xyz) as xyz:
                f.writelines(xyz.readlines()[2:])
            f.write("*\n")

        energy = None
        status = "unknown"

        try:
            with open(log_out, "w") as log:
                subprocess.run(["orca", str(inp_file)],
                               stdout=log, stderr=log, check=True)

            final_xyz = f"{base}.xyz"
            if os.path.exists(final_xyz):
                shutil.copy(final_xyz, xyz_out)
                status = "converged"

            # Parse energy
            with open(log_out) as log:
                for line in log:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        try:
                            energy = float(line.split()[-1])
                        except Exception:
                            pass
                        break

        except subprocess.CalledProcessError:
            status = "failed"

        return {
            "energy": energy,
            "status": status,
            "xyz_out": xyz_out if os.path.exists(xyz_out) else None,
            "log_out": log_out,
            "backend_meta": {"mode": mode}
        }
