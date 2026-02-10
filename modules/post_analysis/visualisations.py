"""
visualisations.py

Plotting utilities for post-analysis of the openCOSMO-RS pipeline.
All functions are explicit, minimal, and reviewer-friendly.
"""

import json
import pathlib
from collections import defaultdict
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt


def _load_metadata(metadata_dir: str) -> Dict[str, dict]:
    """
    Load all molecule metadata JSON files from a directory.
    Returns a dict keyed by InChIKey.
    """
    meta_index = {}
    meta_path = pathlib.Path(metadata_dir)

    for file in meta_path.glob("*.json"):
        with open(file) as f:
            data = json.load(f)
            ikey = data.get("inchi_key")
            if ikey:
                meta_index[ikey] = data

    return meta_index


def _extract_gxtb_records(
    energies: List[dict],
    metadata: Dict[str, dict]
) -> List[Tuple[int, float]]:
    """
    Extract (rotatable_bonds, elapsed_seconds) pairs for each conformer
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

        # Find the gXTB stage
        gxtb_stage = next(
            (s for s in entry.get("optimisation_history", [])
             if s.get("engine") == "gxtb"),
            None
        )
        if not gxtb_stage:
            continue

        elapsed = gxtb_stage.get("elapsed_seconds")
        if elapsed is None:
            continue

        records.append((rot_bonds, elapsed))

    return records


def plot_gxtb_time_vs_rotatable_bonds(
    energies_path: str,
    metadata_dir: str,
    aggregate: bool = True
):
    """
    Plot gXTB optimisation elapsed time as a function of rotatable bonds.

    Parameters
    ----------
    energies_path : str
        Path to the energies.json file.
    metadata_dir : str
        Directory containing molecule metadata JSON files.
    aggregate : bool
        If True, overlays mean elapsed time per rotatable bond count.
    """
    # Load metadata
    metadata = _load_metadata(metadata_dir)

    # Load energies
    with open(energies_path) as f:
        energies = json.load(f)

    # Extract paired data
    records = _extract_gxtb_records(energies, metadata)
    if not records:
        raise ValueError("No gXTB optimisation records found.")

    rot, time = zip(*records)

    # Raw scatter
    plt.figure(figsize=(9, 6))
    plt.scatter(rot, time, alpha=0.6, label="Conformers")

    # Aggregated means
    if aggregate:
        groups = defaultdict(list)
        for rb, t in records:
            groups[rb].append(t)

        sorted_rb = sorted(groups.keys())
        mean_time = [sum(groups[rb]) / len(groups[rb]) for rb in sorted_rb]

        plt.plot(
            sorted_rb,
            mean_time,
            marker="o",
            color="red",
            linewidth=2,
            label="Mean per rotatable bond count"
        )

    plt.xlabel("Rotatable Bonds")
    plt.ylabel("gXTB Elapsed Time (s)")
    plt.title("gXTB Optimisation Time vs Rotatable Bonds")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()
