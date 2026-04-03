"""
modules/utils/rdkit_helper.py

Shared RDKit helpers for XYZ file I/O and mol-from-XYZ with bond perception.

Used by:
  - modules/stages/optimisation_stage.py  (_backend_forcefield)
  - modules/utils/molecule_utils.py       (mol_from_xyz)
"""

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Geometry import Point3D


_COV_RADIUS_FACTOR = 1.3  # multiplier on covalent radii for distance-based fallback


def load_conformer_from_xyz(filepath: str):
    """
    Read an XYZ file and return (Conformer, atom_symbols, comment).

    Does not attempt bond perception — use load_mol_from_xyz for that.
    """
    with open(filepath, "r") as f:
        n_atoms      = int(f.readline().strip())
        comment      = f.readline().strip()
        atom_symbols = []
        conformer    = Chem.Conformer(n_atoms)
        for i in range(n_atoms):
            parts = f.readline().strip().split()
            atom_symbols.append(parts[0])
            conformer.SetAtomPosition(
                i, Point3D(float(parts[1]), float(parts[2]), float(parts[3]))
            )
    return conformer, atom_symbols, comment


def load_mol_from_xyz(xyz_file: str, charge: int = 0) -> Chem.Mol:
    """
    Build an RDKit Mol from an XYZ file with automatic bond perception.

    Strategy:
      1. Try Chem.MolFromXYZFile + rdDetermineBonds.DetermineBonds.
      2. On failure, fall back to distance-based bonding using covalent radii
         scaled by _COV_RADIUS_FACTOR (1.3).  All bonds assigned as SINGLE;
         adequate for forcefield setup which will re-assign orders internally.
    """
    mol = None
    try:
        mol = Chem.MolFromXYZFile(xyz_file)
        rdDetermineBonds.DetermineBonds(mol, charge=charge)
    except Exception:
        mol = None

    if mol is None:
        conformer, atom_symbols, _ = load_conformer_from_xyz(xyz_file)
        rw = Chem.RWMol()
        for sym in atom_symbols:
            rw.AddAtom(Chem.Atom(sym))
        rw.AddConformer(conformer)

        pt    = Chem.GetPeriodicTable()
        d_mat = Chem.Get3DDistanceMatrix(rw)
        n     = rw.GetNumAtoms()
        for i in range(n):
            rc_i = pt.GetRcovalent(rw.GetAtomWithIdx(i).GetAtomicNum()) * _COV_RADIUS_FACTOR
            for j in range(i + 1, n):
                rc_j = pt.GetRcovalent(rw.GetAtomWithIdx(j).GetAtomicNum()) * _COV_RADIUS_FACTOR
                if d_mat[i, j] <= rc_i + rc_j:
                    rw.AddBond(i, j, Chem.BondType.SINGLE)

        mol = rw.GetMol()

    return mol


def write_xyz(mol: Chem.Mol, path: str, comment: str = "optimised") -> None:
    """Write an RDKit mol's conformer to an XYZ file."""
    conf = mol.GetConformer()
    with open(path, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n")
        f.write(f"{comment}\n")
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            f.write(f"{atom.GetSymbol()}  {pos.x:.16f}  {pos.y:.16f}  {pos.z:.16f}\n")
