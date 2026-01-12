# modules/utils/output_parsers.py

import os
from dataclasses import dataclass, field
from typing import List, Tuple


# ============================================================
#  RDKit Parser
# ============================================================

class RdkitParser:
    def __init__(self, path):
        self.path = path

    def get_xyz(self):
        with open(self.path) as f:
            lines = f.readlines()
        return "".join(lines[2:])  # skip atom count + comment

    def get_energy(self):
        with open(self.path) as f:
            second_line = f.readlines()[1]
        if "energy=" in second_line:
            return float(second_line.split("energy=")[1])
        return None

    def get_metadata(self):
        return {"engine": "rdkit"}


# ============================================================
#  XTB LOG Parser (not COSMO)
# ============================================================

class XtbLogParser:
    def __init__(self, path):
        self.path = path

    def get_energy(self):
        with open(self.path) as f:
            for line in f:
                if "TOTAL ENERGY" in line:
                    return float(line.split()[-1])
        raise RuntimeError("Energy not found in XTB log")

    def get_xyz(self):
        coords = []
        capture = False
        with open(self.path) as f:
            for line in f:
                if line.strip().startswith("$coord"):
                    capture = True
                    continue
                if capture:
                    if line.strip().startswith("$end"):
                        break
                    coords.append(line)
        return "".join(coords)

    def get_metadata(self):
        return {"engine": "xtb"}


# ============================================================
#  ORCA Parser
# ============================================================

class OrcaParser:
    def __init__(self, path):
        self.path = path

    def get_energy(self):
        with open(self.path) as f:
            for line in f:
                if "FINAL SINGLE POINT ENERGY" in line:
                    return float(line.split()[-1])
        raise RuntimeError("Energy not found in ORCA log")

    def get_xyz(self):
        coords = []
        capture = False
        with open(self.path) as f:
            for line in f:
                if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                    capture = True
                    next(f)
                    next(f)
                    continue
                if capture:
                    if line.strip() == "":
                        break
                    coords.append(line)
        return "".join(coords)

    def get_metadata(self):
        return {"engine": "orca"}


# ============================================================
#  XTB COSMO Parser
# ============================================================

@dataclass
class CosmoAtom:
    index: int
    element: str
    radius: float
    position: Tuple[float, float, float]


@dataclass
class CosmoSegment:
    index: int
    atom_index: int
    position: Tuple[float, float, float]
    charge: float
    area: float
    sigma_raw: float
    potential: float


@dataclass
class XtbCosmoData:
    atoms: List[CosmoAtom] = field(default_factory=list)
    segments: List[CosmoSegment] = field(default_factory=list)
    area: float | None = None
    fepsi: float | None = None
    dielectric_energy: float | None = None
    total_energy: float | None = None
    raw_text: str | None = None


class XtbCosmoParser:
    """
    Reads a TURBOMOLE-format COSMO file produced by XTB.
    """

    def __init__(self, path):
        self.path = path

    def get_metadata(self):
        return {"engine": "xtb-cosmo"}

    def get_energy(self):
        return self.get_cosmo_data().total_energy

    def get_xyz(self):
        raise RuntimeError("COSMO files do not contain XYZ geometry.")

    def get_cosmo_data(self) -> XtbCosmoData:
        with open(self.path, "r") as f:
            lines = f.readlines()

        data = XtbCosmoData(raw_text="".join(lines))

        i = 0
        n = len(lines)

        while i < n:
            line = lines[i].strip()

            # --------------------------
            # $cosmo_data
            # --------------------------
            if line.startswith("$cosmo_data"):
                i += 1
                while i < n and not lines[i].startswith("$"):
                    if "area=" in lines[i]:
                        data.area = float(lines[i].split("=")[1])
                    if "fepsi=" in lines[i]:
                        data.fepsi = float(lines[i].split("=")[1])
                    i += 1
                continue

            # --------------------------
            # $coord_rad
            # --------------------------
            if line.startswith("$coord_rad"):
                i += 2
                while i < n and lines[i].strip() and not lines[i].startswith("$"):
                    parts = lines[i].split()
                    if len(parts) >= 6:
                        idx = int(parts[0])
                        x, y, z = map(float, parts[1:4])
                        element = parts[4]
                        radius = float(parts[5])
                        data.atoms.append(
                            CosmoAtom(idx, element, radius, (x, y, z))
                        )
                    i += 1
                continue

            # --------------------------
            # $segment_information
            # --------------------------
            if line.startswith("$segment_information"):
                while i < n and not lines[i].lstrip().startswith("1"):
                    i += 1
                while i < n and lines[i].strip() and not lines[i].startswith("$"):
                    parts = lines[i].split()
                    if len(parts) >= 8:
                        seg_idx = int(parts[0])
                        atom_idx = int(parts[1])
                        x, y, z = map(float, parts[2:5])
                        charge = float(parts[5])
                        area = float(parts[6])
                        sigma_raw = float(parts[7])
                        potential = float(parts[8]) if len(parts) > 8 else 0.0
                        data.segments.append(
                            CosmoSegment(
                                seg_idx, atom_idx, (x, y, z),
                                charge, area, sigma_raw, potential
                            )
                        )
                    i += 1
                continue

            # --------------------------
            # $cosmo_energy
            # --------------------------
            if line.startswith("$cosmo_energy"):
                i += 1
                while i < n and not lines[i].startswith("$"):
                    if "Total energy corrected" in lines[i]:
                        data.total_energy = float(lines[i].split("=")[1])
                    if "Dielectric energy" in lines[i]:
                        data.dielectric_energy = float(lines[i].split("=")[1])
                    i += 1
                continue

            i += 1

        return data


# ============================================================
#  XTB Topology Parser (xtbtopo.mol)
# ============================================================

@dataclass
class XtbTopology:
    atoms: List[str]
    adjacency_list: List[Tuple[int, int, int]]
    adjacency_matrix: List[List[int]]


class XtbTopologyParser:
    """
    Parses xtbtopo.mol (V2000 MOL format) to extract adjacency information.
    """

    def __init__(self, path):
        self.path = path

    def get_topology_data(self) -> XtbTopology:
        atoms = []
        adjacency_list = []

        with open(self.path, "r") as f:
            lines = f.readlines()

        # MOL format: atom/bond counts on line 4
        counts_line = lines[3]
        natoms = int(counts_line[0:3])
        nbonds = int(counts_line[3:6])

        # Atom block
        atom_start = 4
        atom_end = atom_start + natoms

        for line in lines[atom_start:atom_end]:
            parts = line.split()
            element = parts[3]
            atoms.append(element)

        # Bond block
        bond_start = atom_end
        bond_end = bond_start + nbonds

        for line in lines[bond_start:bond_end]:
            parts = line.split()
            i = int(parts[0])
            j = int(parts[1])
            order = int(parts[2])
            adjacency_list.append((i, j, order))

        # Build adjacency matrix
        adj = [[0 for _ in range(natoms)] for _ in range(natoms)]
        for i, j, order in adjacency_list:
            adj[i-1][j-1] = order
            adj[j-1][i-1] = order

        return XtbTopology(
            atoms=atoms,
            adjacency_list=adjacency_list,
            adjacency_matrix=adj
        )
