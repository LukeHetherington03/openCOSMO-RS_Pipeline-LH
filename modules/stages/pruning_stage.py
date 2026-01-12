import os
import json
import time
import shutil
import pandas as pd

from modules.utils.atomic_write import AtomicWriter


class PruningStage:
    """
    Pruning stage:
      - Takes an energies.json from a previous stage
      - Copies it into inputs/
      - Copies all referenced XYZs into inputs/xyz/
      - Applies pruning strategy
      - Writes pruned energies.json, summary.csv, xyz/ into outputs/
    """

    # ------------------------------------------------------------
    # Entry point (clean + readable)
    # ------------------------------------------------------------
    def run(self, job):
        params = job.parameters

        summary_file = params.get("summary_file")
        strategy = params.get("strategy")
        strategy_params = params.get("strategy_params", {})

        if summary_file is None:
            raise ValueError("PruningStage requires 'summary_file' parameter.")

        job.log_header("Starting pruning stage")
        job.log(f"Strategy: {strategy} {strategy_params}")
        job.log(f"Input energies: {summary_file}")

        # Prepare directories + load entries
        inputs, entries = self._prepare_inputs(job, summary_file)

        # Group by molecule
        groups = self._group_by_molecule(entries)
        total_confs = len(entries)
        total_mols = len(groups)

        job.log(f"Input: {total_confs} conformers across {total_mols} molecules")

        # Apply pruning
        pruned = self._prune_all(job, groups, strategy, strategy_params)

        # Write outputs
        self._write_outputs(job, pruned, params, inputs["missing"])

        job.log_header("Pruning complete")
        job.log(f"{len(pruned)} conformers retained")
        job.log(f"Outputs written to: {job.outputs_dir}")

        job.mark_complete()

    # ------------------------------------------------------------
    # Input preparation
    # ------------------------------------------------------------
    def _prepare_inputs(self, job, summary_file):
        inputs_dir = job.inputs_dir
        outputs_dir = job.outputs_dir

        inputs_xyz_dir = os.path.join(inputs_dir, "xyz")
        outputs_xyz_dir = os.path.join(outputs_dir, "xyz")

        os.makedirs(inputs_xyz_dir, exist_ok=True)
        os.makedirs(outputs_xyz_dir, exist_ok=True)

        # Copy energies.json
        local_energies_path = os.path.join(inputs_dir, "energies.json")
        shutil.copy(summary_file, local_energies_path)

        # Load entries
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

        # Filter missing
        entries = [e for e in entries if e["lookup_id"] not in missing]

        if not entries:
            raise ValueError("No valid conformers available for pruning.")

        return {
            "inputs_dir": inputs_dir,
            "outputs_dir": outputs_dir,
            "missing": missing,
            "outputs_xyz_dir": outputs_xyz_dir,
        }, entries

    # ------------------------------------------------------------
    # Group by molecule
    # ------------------------------------------------------------
    def _group_by_molecule(self, entries):
        groups = {}
        for e in entries:
            mol_id = e["lookup_id"].split("_conf")[0]
            groups.setdefault(mol_id, []).append(e)
        return groups

    # ------------------------------------------------------------
    # Apply pruning to all molecules
    # ------------------------------------------------------------
    def _prune_all(self, job, groups, strategy, strategy_params):
        pruned = []
        total_mols = len(groups)

        for idx, (mol_id, conformers) in enumerate(groups.items(), start=1):
            job.log_section(f"Molecule {idx}/{total_mols}: {mol_id}")

            mol_start = time.perf_counter()
            kept = self._apply_strategy(conformers, strategy, strategy_params)
            elapsed = time.perf_counter() - mol_start

            job.log(f"{len(conformers)} → {len(kept)} retained", indent=1)
            job.log(f"Completed in {elapsed:.2f} seconds", indent=1)

            pruned.extend(kept)

        return pruned

    # ------------------------------------------------------------
    # Strategy dispatcher
    # ------------------------------------------------------------
    def _apply_strategy(self, conformers, strategy, params):
        if strategy == "top_n":
            return self._top_n(conformers, params.get("n"))

        if strategy == "energy_window":
            return self._energy_window(conformers, params.get("threshold"))

        if strategy == "rmsd":
            return conformers  # placeholder

        raise ValueError(f"Unknown pruning strategy: {strategy}")

    # ------------------------------------------------------------
    # Strategies
    # ------------------------------------------------------------
    def _top_n(self, conformers, n):
        if n is None:
            raise ValueError("top_n strategy requires parameter 'n'")
        return sorted(conformers, key=lambda e: e["energy"])[:n]

    def _energy_window(self, conformers, threshold):
        if threshold is None:
            raise ValueError("energy_window strategy requires 'threshold'")
        sorted_group = sorted(conformers, key=lambda e: e["energy"])
        min_energy = sorted_group[0]["energy"]
        return [e for e in sorted_group if e["energy"] <= min_energy + threshold]

    # ------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------
    def _write_outputs(self, job, pruned, params, missing):
        outputs_dir = job.outputs_dir
        outputs_xyz_dir = os.path.join(outputs_dir, "xyz")

        # Copy XYZs
        for entry in pruned:
            src = entry["xyz_path"]
            dst = os.path.join(outputs_xyz_dir, os.path.basename(src))
            shutil.copy(src, dst)
            entry["xyz_path"] = dst

        # summary.csv
        summary_path = os.path.join(outputs_dir, "summary.csv")
        df = pd.DataFrame(
            [
                {
                    "lookup_id": e["lookup_id"],
                    "energy": e["energy"],
                    "xyz_path": e["xyz_path"],
                    "log_path": e["log_path"],
                }
                for e in pruned
            ]
        )
        df.to_csv(summary_path, index=False)

        # energies.json
        pruned_energies_path = os.path.join(outputs_dir, "energies.json")
        with AtomicWriter(pruned_energies_path) as f:
            json.dump(pruned, f, indent=2)

        # job_state.json
        job_state_path = os.path.join(job.job_dir, "job_state.json")
        with AtomicWriter(job_state_path) as f:
            json.dump(
                {
                    "stage": "pruning",
                    "strategy": params.get("strategy"),
                    "strategy_params": params.get("strategy_params"),
                    "input_energies": params.get("summary_file"),
                    "num_input": len(pruned) + len(missing),
                    "num_output": len(pruned),
                    "missing_xyz": missing,
                },
                f,
                indent=2,
            )
