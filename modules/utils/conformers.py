import json
from typing import List, Dict, Tuple, Optional

from modules.utils.atomic_write import AtomicWriter

class ConformerRecord:
    """
    Represents a single conformer with a stable schema used across
    Generation, Pruning, and Optimisation stages.
    """

    def __init__(
        self,
        inchi_key: str,
        conf_num: int,
        xyz_path: str,
        energy: Optional[float],
        smiles: str,
        provenance: Optional[Dict] = None,
        optimisation_history: Optional[List[Dict]] = None,
    ):
        self.inchi_key = inchi_key
        self.conf_num = int(conf_num)
        self.lookup_id = f"{inchi_key}_conf{conf_num:03d}"
        self.xyz_path = xyz_path
        self.energy = energy
        self.smiles = smiles

        # Misc metadata
        self.provenance = provenance or {}

        # Chronological list of optimisation steps
        self.optimisation_history = optimisation_history or []

    # ------------------------------------------------------------
    # Serialisation helpers
    # ------------------------------------------------------------
    def to_dict(self) -> Dict:
        return {
            "lookup_id": self.lookup_id,
            "inchi_key": self.inchi_key,
            "conf_num": self.conf_num,
            "xyz_path": self.xyz_path,
            "energy": self.energy,
            "smiles": self.smiles,
            "provenance": self.provenance,
            "optimisation_history": self.optimisation_history,
        }

    @classmethod
    def from_dict(cls, data: Dict):
        return cls(
            inchi_key=data["inchi_key"],
            conf_num=data["conf_num"],
            xyz_path=data["xyz_path"],
            energy=data.get("energy"),
            smiles=data.get("smiles", ""),
            provenance=data.get("provenance", {}),
            optimisation_history=data.get("optimisation_history", []),
        )

class ConformerSet:
    """
    A container for a list of ConformerRecord objects.
    Provides loading/saving and convenience operations.
    """

    def __init__(self, records: Optional[List[ConformerRecord]] = None):
        self.records: List[ConformerRecord] = records or []

    # ------------------------------------------------------------
    # Basic operations
    # ------------------------------------------------------------
    def add(self, record: ConformerRecord):
        self.records.append(record)

    def __len__(self):
        return len(self.records)

    def __iter__(self):
        return iter(self.records)

    # ------------------------------------------------------------
    # Grouping
    # ------------------------------------------------------------
    def group_by_molecule(self) -> Dict[str, List[ConformerRecord]]:
        grouped = {}
        for rec in self.records:
            grouped.setdefault(rec.inchi_key, []).append(rec)
        return grouped

    # ------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------
    def to_list(self) -> List[Dict]:
        return [rec.to_dict() for rec in self.records]

    @classmethod
    def from_list(cls, data: List[Dict]):
        return cls(records=[ConformerRecord.from_dict(d) for d in data])

    # ------------------------------------------------------------
    # File IO
    # ------------------------------------------------------------
    @classmethod
    def load(cls, path: str):
        with open(path, "r") as f:
            data = json.load(f)
        return cls.from_list(data)

    def save(self, path: str):
        with AtomicWriter(path) as f:
            json.dump(self.to_list(), f, indent=2)
