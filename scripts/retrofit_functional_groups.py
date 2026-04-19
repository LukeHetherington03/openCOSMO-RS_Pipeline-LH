"""
retrofit_functional_groups.py

Adds functional_groups counts to existing molecule_metadata JSON files
without re-running the full pipeline cleaning stage.

Works on any molecule_metadata/ directory — pass one or more paths as
arguments, or run with no arguments to use the two standard locations:
  1. CONSTANT_FILES/molecule_metadata/
  2. All final_outputs/molecule_metadata/ dirs found under pipeline_data/

Usage:
    python3 scripts/retrofit_functional_groups.py
    python3 scripts/retrofit_functional_groups.py path/to/molecule_metadata/

Behaviour:
  - Reads SMILES from each JSON, computes functional group counts via RDKit
  - Adds a "functional_groups" key and bumps "metadata_version" to 4
  - Skips files that already have metadata_version >= 4 (already retrofitted)
  - Writes back to the same file (in-place)
  - Prints a summary of files updated / skipped / failed
"""

import json
import sys
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Import the canonical SMARTS dict and version from the cleaning stage
from modules.stages.cleaning_stage import FUNCTIONAL_GROUP_SMARTS, METADATA_VERSION

# Pre-compile patterns once
_PATTERNS = {}
for name, smarts in FUNCTIONAL_GROUP_SMARTS.items():
    pat = Chem.MolFromSmarts(smarts)
    if pat is None:
        print(f"WARNING: could not compile SMARTS for '{name}': {smarts}")
    _PATTERNS[name] = pat


def compute_fg(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {k: None for k in FUNCTIONAL_GROUP_SMARTS}
    return {
        name: len(mol.GetSubstructMatches(pat)) if pat is not None else 0
        for name, pat in _PATTERNS.items()
    }


def retrofit_directory(meta_dir: Path) -> dict:
    counts = {"updated": 0, "already_current": 0, "no_smiles": 0, "failed": 0}
    json_files = sorted(meta_dir.glob("*.json"))
    if not json_files:
        print(f"  No JSON files found in {meta_dir}")
        return counts

    for path in json_files:
        try:
            with open(path) as f:
                meta = json.load(f)

            version = meta.get("metadata_version", 0)
            if version >= METADATA_VERSION and "functional_groups" in meta:
                counts["already_current"] += 1
                continue

            smiles = meta.get("smiles")
            if not smiles:
                counts["no_smiles"] += 1
                continue

            meta["functional_groups"] = compute_fg(smiles)
            meta["metadata_version"]  = METADATA_VERSION

            with open(path, "w") as f:
                json.dump(meta, f, indent=2)

            counts["updated"] += 1

        except Exception as e:
            print(f"  ERROR: {path.name}: {e}")
            counts["failed"] += 1

    return counts


def find_all_metadata_dirs() -> list[Path]:
    dirs = []

    # 1. CONSTANT_FILES
    const = ROOT / "CONSTANT_FILES" / "molecule_metadata"
    if const.is_dir():
        dirs.append(const)

    # 2. All final_outputs under pipeline_data
    for p in sorted((ROOT / "pipeline_data").rglob("final_outputs/molecule_metadata")):
        if p.is_dir():
            dirs.append(p)

    return dirs


def main():
    if len(sys.argv) > 1:
        target_dirs = [Path(a) for a in sys.argv[1:]]
    else:
        target_dirs = find_all_metadata_dirs()

    if not target_dirs:
        print("No molecule_metadata directories found.")
        sys.exit(1)

    print(f"Retrofitting functional groups (METADATA_VERSION → {METADATA_VERSION})")
    print(f"Patterns: {list(FUNCTIONAL_GROUP_SMARTS.keys())}\n")

    total = {"updated": 0, "already_current": 0, "no_smiles": 0, "failed": 0}

    for d in target_dirs:
        print(f"Directory: {d}")
        result = retrofit_directory(d)
        for k in total:
            total[k] += result[k]
        print(f"  updated={result['updated']}  already_current={result['already_current']}"
              f"  no_smiles={result['no_smiles']}  failed={result['failed']}")

    print(f"\nTotal — updated: {total['updated']}  skipped: {total['already_current']}"
          f"  no_smiles: {total['no_smiles']}  failed: {total['failed']}")


if __name__ == "__main__":
    main()
