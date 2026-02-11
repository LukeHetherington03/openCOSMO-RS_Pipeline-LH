#!/usr/bin/env python3

"""
run_visualisations.py

Generates visualisations for multiple optimisation pipelines:
- gXTB only
- gXTB → ORCA (DFT)
- ORCA only

SSH-safe: saves PNGs instead of showing plots.
Elapsed time is converted from seconds → hours.
"""

import json
import pathlib
from collections import defaultdict
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt


# ------------------------------------------------------------
# Metadata loader
# ------------------------------------------------------------

def load_metadata(metadata_dir: str) -> Dict[str, dict]:
    meta_index = {}
    meta_path = pathlib.Path(metadata_dir)

    for file in meta_path.glob("*.json"):
        with open(file) as f:
            data = json.load(f)
            ikey = data.get("inchi_key")
            if ikey:
                meta_index[ikey] = data

    return meta_index


# ------------------------------------------------------------
# Extract times for gXTB, ORCA, and combined
# ------------------------------------------------------------

def extract_times(
    energies: List[dict],
    metadata: Dict[str, dict]
) -> List[Tuple[int, float, float, float]]:
    """
    Returns a list of tuples:
    (rotatable_bonds, gxtb_hours, orca_hours, combined_hours)
    Missing stages are recorded as 0 hours.
    """
    records = []

    for entry in energies:
        ikey = entry.get("inchi_key")
        meta = metadata.get(ikey)
        if not meta:
            continue

        rot_bonds = meta.get("rotatable_bonds")
        if rot_bonds is None:
            continue

        gxtb_seconds = 0.0
        orca_seconds = 0.0

        for stage in entry.get("optimisation_history", []):
            engine = stage.get("engine")

            if engine == "gxtb":
                gxtb_seconds = stage.get("elapsed_seconds", 0.0)

            if engine == "orca":
                orca_seconds = stage.get("elapsed_seconds", 0.0)

        gxtb_hours = gxtb_seconds / 3600.0
        orca_hours = orca_seconds / 3600.0
        combined_hours = gxtb_hours + orca_hours

        records.append((rot_bonds, gxtb_hours, orca_hours, combined_hours))

    return records


# ------------------------------------------------------------
# Plotting function
# ------------------------------------------------------------

def plot_multi_pipeline_times(
    energies_paths: Dict[str, str],
    metadata_dir: str,
    output_dir: str
):
    """
    energies_paths: dict with keys:
        "gxtb", "gxtb_dft", "dft"
        each value is a path to an energies.json file
    """

    output_path = pathlib.Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    metadata = load_metadata(metadata_dir)

    # Storage for plotting
    pipeline_data = {}

    for label, path in energies_paths.items():
        with open(path) as f:
            energies = json.load(f)

        records = extract_times(energies, metadata)
        if not records:
            raise ValueError(f"No optimisation records found for {label}")

        pipeline_data[label] = records

    # Plot
    plt.figure(figsize=(10, 7))

    colours = {
        "gxtb": "blue",
        "gxtb_dft": "green",
        "dft": "red"
    }

    labels = {
        "gxtb": "gXTB only",
        "gxtb_dft": "gXTB → ORCA",
        "dft": "ORCA only"
    }

    for key, records in pipeline_data.items():
        rot = [r[0] for r in records]

        if key == "gxtb":
            y = [r[1] for r in records]  # gxtb_hours
        elif key == "dft":
            y = [r[2] for r in records]  # orca_hours
        else:
            y = [r[3] for r in records]  # combined_hours

        # Aggregate mean per rotatable bond count
        groups = defaultdict(list)
        for rb, val in zip(rot, y):
            groups[rb].append(val)

        sorted_rb = sorted(groups.keys())
        mean_vals = [sum(groups[rb]) / len(groups[rb]) for rb in sorted_rb]

        plt.plot(
            sorted_rb,
            mean_vals,
            marker="o",
            linewidth=2,
            color=colours[key],
            label=labels[key]
        )

    plt.xlabel("Rotatable Bonds")
    plt.ylabel("Elapsed Time (hours)")
    plt.title("Optimisation Time vs Rotatable Bonds (gXTB, ORCA, Combined)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()

    out_file = output_path / "multi_pipeline_times.png"
    plt.savefig(out_file, dpi=300)
    plt.close()

    print(f"Saved plot to: {out_file}")


# ------------------------------------------------------------
# INPUTS — edit these to run your analysis
# ------------------------------------------------------------

if __name__ == "__main__":
    ENERGIES = {
        "gxtb": "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/post_analysis_results/inputs/energies_gxtb.json",
        "gxtb_dft": "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/post_analysis_results/inputs/energies_gxtb_dft.json",
        "dft": "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/post_analysis_results/inputs/energies_dft.json"
    }

    METADATA_DIR = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/CONSTANT_FILES/molecule_metadata"
    OUTPUT_DIR = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/post_analysis_results/visualisations"

    plot_multi_pipeline_times(
        energies_paths=ENERGIES,
        metadata_dir=METADATA_DIR,
        output_dir=OUTPUT_DIR
    )
