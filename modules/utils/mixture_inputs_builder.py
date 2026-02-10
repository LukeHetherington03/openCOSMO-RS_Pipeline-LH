#!/usr/bin/env python3
from pathlib import Path


def build_mixture_inputs(
    solute_meta: dict,
    solute_dir: Path,
    solvent_dir: Path,
    n_solute_confs: int,
    n_solvent_confs: int,
    defaults: dict,
    temperature: float,
) -> str:
    """
    Build a valid mixture_inputs.txt string for the legacy COSMO‑RS script.
    """

    # Metadata
    inchi = solute_meta["inchi_key"]
    smiles = solute_meta["smiles"]
    Tm = solute_meta["melting_temp"]

    # Defaults
    calc_type = defaults["calculation_type"]
    initial_solute_x = defaults["initial_solute_x"]
    solvent_x = 1.0 - initial_solute_x

    Gfus_mode = defaults["saturation"]["Gfus_mode"]
    Hfus = defaults["saturation"]["Hfus"]
    w = defaults["saturation"]["SORcf"]

    # Multiplicities (assume 1 for all conformers)
    solute_mults = " ".join(["1"] * n_solute_confs)
    solvent_mults = " ".join(["1"] * n_solvent_confs)

    text = f"""calc_type {calc_type}
T {temperature}
sat {inchi} {smiles}
Tm {Tm}
Gfus {Gfus_mode}
Hfus {Hfus}
w {w}

{inchi} {initial_solute_x} {solute_dir} {n_solute_confs} {solute_mults}
{defaults["default_solvent"]} {solvent_x} {solvent_dir} {n_solvent_confs} {solvent_mults}
"""

    return text
