import json
from pathlib import Path

BASE = Path("/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/CONSTANT_FILES/molecule_metadata")

def load_metadata(inchi_key: str):
    path = BASE / f"{inchi_key}.json"

    with open(path, "r") as f:
        meta = json.load(f)

    # Normalise melting temp
    Tm = meta.get("melting_temp", "N/A")
    try:
        Tm = float(Tm)
    except Exception:
        pass  # "liquid", "N/A", etc.

    return {
        "inchi_key": meta["inchi_key"],
        "smiles": meta.get("smiles", ""),
        "Tm": Tm,
        "Gfus_mode": meta.get("Gfus_mode", "MyrdalYalkowsky"),
        "Hfus": meta.get("Hfus", "N/A"),
        "experimental_solubility_mol_frac": meta.get("experimental_solubility_mol_frac", None),
        "mol_name": meta.get("mol_name", meta["inchi_key"]),
    }
