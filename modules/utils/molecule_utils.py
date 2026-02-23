from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops
from rdkit.Chem import rdchem
from rdkit.Chem import rdDetermineBonds
import requests
import re

class MoleculeUtils:
    """
    Utility functions for RDKit molecule handling.
    These operate directly on RDKit Mol objects and do NOT parse
    external engine output (that belongs in output_parsers).
    """

    @staticmethod
    def rotatable_bonds(mol):
        """
        Count rotatable bonds using RDKit's definition.
        Returns an integer. If RDKit fails, returns 0.
        """
        try:
            return rdMolDescriptors.CalcNumRotatableBonds(mol)
        except Exception:
            return 0

    @staticmethod
    def mol_to_xyz_block(mol, conf_id=0):
        """
        Convert an RDKit molecule + conformer into an XYZ block (no header).
        Returns a string like:

            C   0.000000   1.234567   2.345678
            H   0.123456   1.567890   2.678901
            ...

        This is used by the generation stage to write XYZ files.
        """
        try:
            conf = mol.GetConformer(conf_id)
        except Exception:
            raise ValueError(f"Conformer ID {conf_id} not found in molecule.")

        lines = []
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)
            symbol = atom.GetSymbol()
            lines.append(f"{symbol:2s}  {pos.x: .6f}  {pos.y: .6f}  {pos.z: .6f}")

        return "\n".join(lines) + "\n"




    def get_molecule_info_from_inchi(inchi: str):
        """
        Returns RDKit molecule and adjacency matrix from an InChI string.
        """

        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            raise ValueError(f"RDKit failed to parse InChI: {inchi}")

        # Sanitize (perceive bonds, valence, aromaticity)
        Chem.SanitizeMol(mol)

        # Compute adjacency with bond orders
        adj = rdmolops.GetAdjacencyMatrix(mol, useBO=True)

        return {
            "mol": mol,
            "adjacency": adj.astype(int).tolist(),
            "atomic_numbers": [a.GetAtomicNum() for a in mol.GetAtoms()],
            "inchi": inchi,
            "smiles": Chem.MolToSmiles(mol),
    }


    @staticmethod
    def adjacency_from_smiles(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        n = mol.GetNumAtoms()
        adj = [[0]*n for _ in range(n)]

        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            adj[i][j] = 1
            adj[j][i] = 1

        return adj



    @staticmethod
    def mol_from_xyz(xyz_atoms):
        """
        Build an RDKit molecule from XYZ atoms.
        xyz_atoms = [(element, x, y, z), ...]
        """
        mol = rdchem.RWMol()
        conf = Chem.Conformer(len(xyz_atoms))

        atom_indices = []

        for i, (el, x, y, z) in enumerate(xyz_atoms):
            atom = Chem.Atom(el)
            idx = mol.AddAtom(atom)
            atom_indices.append(idx)
            conf.SetAtomPosition(idx, (x, y, z))

        mol.AddConformer(conf)

        return mol



    @staticmethod
    def mol_from_xyz(xyz_atoms):
        mol = rdchem.RWMol()
        conf = Chem.Conformer(len(xyz_atoms))

        for i, (el, x, y, z) in enumerate(xyz_atoms):
            idx = mol.AddAtom(Chem.Atom(el))
            conf.SetAtomPosition(idx, (x, y, z))

        mol.AddConformer(conf)

        # Let RDKit guess bonds from geometry
        rdDetermineBonds.DetermineBonds(mol)

        return mol


    @staticmethod
    def adjacency_from_xyz(xyz_atoms):
        mol = MoleculeUtils.mol_from_xyz(xyz_atoms)
        adj = rdmolops.GetAdjacencyMatrix(mol, useBO=True)
        return adj


    # ======================================================================
    # NEW SECTION: Melting point retrieval (PubChem PUG-View)
    # ======================================================================

    PUBCHEM_CID_URL = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}/cids/JSON"
    )
    PUBCHEM_RECORD_URL = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON"
    )

    @staticmethod
    def _pubchem_cid_from_inchikey(inchikey: str):
        """Return PubChem CID for an InChIKey, or None."""
        try:
            url = MoleculeUtils.PUBCHEM_CID_URL.format(inchikey)
            r = requests.get(url, timeout=10)
            data = r.json()
            cids = data.get("IdentifierList", {}).get("CID", [])
            return cids[0] if cids else None
        except Exception:
            return None

    @staticmethod
    def _pubchem_full_record(cid: int):
        """Return full PUG-View JSON record for a CID, or None."""
        try:
            url = MoleculeUtils.PUBCHEM_RECORD_URL.format(cid)
            r = requests.get(url, timeout=10)
            return r.json()
        except Exception:
            return None

    @staticmethod
    def _recursive_find_melting_point(node):
        """
        Recursively search the PUG-View JSON tree for melting point.
        Returns float or None.
        """
        if isinstance(node, dict):
            # Check if this is the melting point section
            if node.get("TOCHeading") == "Melting Point":
                info_list = node.get("Information", [])
                for info in info_list:
                    value = info.get("Value", {})
                    strings = value.get("StringWithMarkup", [])
                    for s in strings:
                        text = s.get("String", "")
                        m = re.search(r"-?\d+(\.\d+)?", text)
                        if m:
                            return float(m.group(0))

            # Recurse into dict values
            for v in node.values():
                mp = MoleculeUtils._recursive_find_melting_point(v)
                if mp is not None:
                    return mp

        elif isinstance(node, list):
            for item in node:
                mp = MoleculeUtils._recursive_find_melting_point(item)
                if mp is not None:
                    return mp

        return None
    
    @staticmethod
    def get_melting_point(inchikey: str):
        """
        Returns melting_temp and melting_temp_source in the format expected
        by the CleaningStage metadata writer.

        Output:
            {
                "melting_temp": float or "N/A",
                "melting_temp_source": "pubchem" or "missing"
            }
        """
        if not inchikey:
            return {"melting_temp": "N/A", "melting_temp_source": "missing"}

        cid = MoleculeUtils._pubchem_cid_from_inchikey(inchikey)
        if cid is None:
            return {"melting_temp": "N/A", "melting_temp_source": "missing"}

        record = MoleculeUtils._pubchem_full_record(cid)
        if not record:
            return {"melting_temp": "N/A", "melting_temp_source": "missing"}

        mp_celsius = MoleculeUtils._recursive_find_melting_point(record)
        if mp_celsius is None:
            return {"melting_temp": "N/A", "melting_temp_source": "missing"}

        # Convert °C → K
        mp_kelvin = mp_celsius + 273.15

        return {"melting_temp": mp_kelvin, "melting_temp_source": "pubchem"}

