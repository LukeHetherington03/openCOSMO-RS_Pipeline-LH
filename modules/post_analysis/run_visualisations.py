#!/usr/bin/env python3

"""
run_visualisations.py

Self-contained script for generating post-analysis visualisations
for the openCOSMO-RS pipeline.

Edit the INPUTS section at the bottom to point to your energies.json,
metadata directory, and output directory.

This version is SSH-safe: it saves PNGs instead of showing plots.
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
# Energies → (rotatable_bonds, elapsed_hours) extractor
# ------------------------------------------------------------

def extract_gxtb_records(
    energies: List[dict],
    metadata: Dict[str, dict]
) -> List[Tuple[int, float]]:
    """
    Extract (rotatable_bonds, elapsed_hours) pairs for each conformer
    that has a gXTB optimisation stage.
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

        gxtb_stage = next(
            (s for s in entry.get("optimisation_history", [])
             if s.get("engine") == "gxtb"),
            None
        )
        if not gxtb_stage:
            continue

        elapsed_seconds = gxtb_stage.get("elapsed_seconds")
        if elapsed_seconds is None:
            continue

        elapsed_hours = elapsed_seconds / 3600.0

        records.append((rot_bonds, elapsed_hours))

    return records


# ------------------------------------------------------------
# Plotting function (SSH-safe)
# ------------------------------------------------------------

def plot_gxtb_time_vs_rotatable_bonds(
    energies_path: str,
    metadata_dir: str,
    output_dir: str,
    aggregate: bool = True
):
    """
    Generate a PNG plot of gXTB optimisation time (hours) vs rotatable bonds.
    Saves the figure to output_dir.
    """

    output_path = pathlib.Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    metadata = load_metadata(metadata_dir)

    with open(energies_path) as f:
        energies = json.load(f)

    records = extract_gxtb_records(energies, metadata)
    if not records:
        raise ValueError("No gXTB optimisation records found.")

    rot, hours = zip(*records)

    plt.figure(figsize=(9, 6))
    plt.scatter(rot, hours, alpha=0.6, label="Conformers")

    if aggregate:
        groups = defaultdict(list)
        for rb, h in records:
            groups[rb].append(h)

        sorted_rb = sorted(groups.keys())
        mean_hours = [sum(groups[rb]) / len(groups[rb]) for rb in sorted_rb]

        plt.plot(
            sorted_rb,
            mean_hours,
            marker="o",
            color="red",
            linewidth=2,
            label="Mean per rotatable bond count"
        )

    plt.xlabel("Rotatable Bonds")
    plt.ylabel("gXTB Elapsed Time (hours)")
    plt.title("gXTB Optimisation Time vs Rotatable Bonds")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()

    out_file = output_path / "gxtb_time_vs_rotatable_bonds_hours.png"
    plt.savefig(out_file, dpi=300)
    plt.close()

    print(f"Saved plot to: {out_file}")


# ------------------------------------------------------------
# INPUTS — edit these to run your analysis
# ------------------------------------------------------------

if __name__ == "__main__":
    ENERGIES_PATH = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/post_analysis_results/inputs/energies.json"
    METADATA_DIR = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/CONSTANT_FILES/molecule_metadata"
    OUTPUT_DIR = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/post_analysis_results/visualisations"

    plot_gxtb_time_vs_rotatable_bonds(
        energies_path=ENERGIES_PATH,
        metadata_dir=METADATA_DIR,
        output_dir=OUTPUT_DIR,
        aggregate=True
    )
