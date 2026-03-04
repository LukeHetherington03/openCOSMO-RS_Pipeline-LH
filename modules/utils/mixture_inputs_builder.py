#!/usr/bin/env python3
"""
mixture_inputs_builder.py
─────────────────────────
Builds a valid mixture_inputs.txt string for the openCOSMO-RS engine.

All solvent paths passed to this builder must point to COPIED conformer
directories inside the job's results tree, never to the source bank in
CONSTANT_FILES.  The SolubilityStage is responsible for copying before
calling this builder.

Melting temperature handling
────────────────────────────
    "liquid" string or numeric Tm <= temperature  ->  meltingtemp liquid
    numeric Tm > temperature                      ->  meltingtemp <value>
    "N/A" string                                  ->  meltingtemp N/A
"""

import re
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# Parameter parsing helpers
# ─────────────────────────────────────────────────────────────────────────────

def parse_solvents_from_args(parameters: dict) -> list:
    """
    Parse solvent specification from a stage parameters dict.

    Returns [{"name": str, "ratio": float}, ...].
    Call normalise_solvent_ratios() to convert ratios to absolute mol fracs.

    Format A — numbered kwargs  [recommended for pipeline_spec]
        solvent1: "water",   solvent1_mol_frac: 0.6
        solvent2: "ethanol", solvent2_mol_frac: 0.4

    Format B — explicit list
        solvents: [{"name": "water", "mol_frac": 0.6}, ...]

    Format C — legacy single solvent
        solvent_name: "water"
    """
    if "solvents" in parameters:
        return [
            {"name": str(s["name"]), "ratio": float(s.get("mol_frac", 1.0))}
            for s in parameters["solvents"]
        ]

    numbered: dict = {}
    for key, val in parameters.items():
        m = re.match(r"^solvent(\d+)$", key)
        if m:
            numbered.setdefault(int(m.group(1)), {})["name"] = str(val)
            continue
        m = re.match(r"^solvent(\d+)_mol_frac$", key)
        if m:
            numbered.setdefault(int(m.group(1)), {})["ratio"] = float(val)

    if numbered:
        result = []
        for n in sorted(numbered):
            entry = numbered[n]
            if "name" not in entry:
                raise ValueError(
                    f"solvent{n}_mol_frac specified but solvent{n} name is missing"
                )
            result.append({"name": entry["name"], "ratio": entry.get("ratio", 1.0)})
        return result

    return [{"name": parameters.get("solvent_name", "water"), "ratio": 1.0}]


def normalise_solvent_ratios(solvents: list, solute_x: float) -> list:
    """
    Scale solvent ratios so they sum to (1.0 - solute_x).

    Returns [{"name": str, "mol_frac": float}, ...].

    Examples
    --------
    [{name:water, ratio:0.6}, {name:ethanol, ratio:0.4}], solute_x=0.01
        -> water=0.594, ethanol=0.396  (total with solute = 1.000)

    [{name:water, ratio:3}, {name:ethanol, ratio:2}], solute_x=0.01
        -> same result — ratios are scaled proportionally
    """
    if not solvents:
        raise ValueError("At least one solvent must be specified")

    total_ratio = sum(s["ratio"] for s in solvents)
    if total_ratio <= 0:
        raise ValueError(f"Solvent ratios sum to {total_ratio} — must be positive")

    solvent_budget = 1.0 - float(solute_x)
    if solvent_budget <= 0:
        raise ValueError(f"initial_solute_x={solute_x} leaves no room for solvents")

    return [
        {
            "name":     s["name"],
            "mol_frac": round(s["ratio"] / total_ratio * solvent_budget, 8),
        }
        for s in solvents
    ]


# ─────────────────────────────────────────────────────────────────────────────
# Builder
# ─────────────────────────────────────────────────────────────────────────────

def build_mixture_inputs(
    solute_meta:     dict,
    solute_dir:      Path,
    solvent_dirs:    dict,
    n_solute_confs:  int,
    n_solvent_confs: dict,
    defaults:        dict,
    temperature:     float,
    solvents:        list,
) -> str:
    """
    Build a valid mixture_inputs.txt string for the openCOSMO-RS engine.

    Parameters
    ----------
    solute_meta : dict
        Must contain: inchi_key, smiles, melting_temp, Gfus_mode, Hfus.

    solute_dir : Path
        Copied solute conformer directory (inside job results tree).

    solvent_dirs : dict[str, Path]
        name -> copied conformer directory (inside job results tree).
        Must NOT point to CONSTANT_FILES bank.

    n_solute_confs : int

    n_solvent_confs : dict[str, int]

    defaults : dict
        Pipeline defaults JSON (calculation_type, initial_solute_x,
        saturation.Gfus_mode, saturation.Hfus, saturation.SORcf).

    temperature : float
        Calculation temperature in Kelvin.

    solvents : list[dict]
        Pre-normalised list from normalise_solvent_ratios():
        [{"name": str, "mol_frac": float}, ...].
    """
    inchi  = solute_meta["inchi_key"]
    smiles = solute_meta["smiles"]
    Tm     = solute_meta["melting_temp"]

    calc_type        = defaults["calculation_type"]
    initial_solute_x = defaults["initial_solute_x"]

    Gfus_mode = solute_meta.get("Gfus_mode") or defaults["saturation"]["Gfus_mode"]
    Hfus      = solute_meta.get("Hfus")      or defaults["saturation"]["Hfus"]
    w         = defaults["saturation"]["SORcf"]

    # Resolve meltingtemp keyword
    if isinstance(Tm, str) and Tm.strip().lower() == "liquid":
        meltingtemp_str = "liquid"
    elif isinstance(Tm, str) and Tm.strip().upper() == "N/A":
        meltingtemp_str = "N/A"
    else:
        try:
            Tm_float = float(Tm)
            meltingtemp_str = "liquid" if Tm_float <= temperature else str(Tm_float)
        except (TypeError, ValueError):
            meltingtemp_str = str(Tm)

    lines = [
        f"calc_type {calc_type}",
        f"temperature {temperature}",
        f"saturation {inchi} {smiles}",
        f"meltingtemp {meltingtemp_str}",
        f"Gfus {Gfus_mode}",
        f"Hfus {Hfus}",
        f"SORcf {w}",
        "",
    ]

    solute_mults = " ".join(["1"] * n_solute_confs)
    lines.append(
        f"{inchi}\t{initial_solute_x}\t{solute_dir}\t{n_solute_confs}\t{solute_mults}"
    )

    for s in solvents:
        name          = s["name"]
        mol_frac      = s["mol_frac"]
        solvent_dir   = solvent_dirs[name]
        n_confs       = n_solvent_confs[name]
        solvent_mults = " ".join(["1"] * n_confs)
        lines.append(
            f"{name}\t{mol_frac}\t{solvent_dir}\t{n_confs}\t{solvent_mults}"
        )

    return "\n".join(lines) + "\n"