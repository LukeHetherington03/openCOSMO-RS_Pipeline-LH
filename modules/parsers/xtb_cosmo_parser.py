# modules/parsers/xtb_cosmo_parser.py

import os
import re
import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional

# Unit conversions
kJ_per_kcal = 4.184
angstrom_per_bohr = 0.52917721092
kJdivmol_per_hartree = 2625.499639479


# ============================================================
# Dataclasses
# ============================================================

@dataclass
class CosmoAtom:
    index: int
    element: str
    position: np.ndarray
    radius: Optional[float] = None


@dataclass
class CosmoSegment:
    index: int
    atom_index: int
    position: np.ndarray
    charge: float
    area: float
    sigma_raw: float
    potential: float


@dataclass
class XtbCosmoData:
    filepath: str
    method: str
    area: float
    volume: float
    energy_tot: float
    energy_dielectric: float
    atoms: List[CosmoAtom] = field(default_factory=list)
    segments: List[CosmoSegment] = field(default_factory=list)

    # Optional sigma-profile fields
    sigma_averaged: Optional[np.ndarray] = None
    sigma_moments: Optional[np.ndarray] = None

    def to_dict(self):
        return {
            "filepath": self.filepath,
            "method": self.method,
            "area": self.area,
            "volume": self.volume,
            "energy_tot": self.energy_tot,
            "energy_dielectric": self.energy_dielectric,
            "n_atoms": len(self.atoms),
            "n_segments": len(self.segments),
            "atoms": [
                {
                    "index": a.index,
                    "element": a.element,
                    "position": a.position.tolist(),
                    "radius": a.radius,
                }
                for a in self.atoms
            ],
            "segments": [
                {
                    "index": s.index,
                    "atom_index": s.atom_index,
                    "position": s.position.tolist(),
                    "charge": s.charge,
                    "area": s.area,
                    "sigma_raw": s.sigma_raw,
                    "potential": s.potential,
                }
                for s in self.segments
            ],
            "sigma_averaged": (
                self.sigma_averaged.tolist()
                if isinstance(self.sigma_averaged, np.ndarray)
                else None
            ),
            "sigma_moments": (
                self.sigma_moments.tolist()
                if isinstance(self.sigma_moments, np.ndarray)
                else None
            ),
        }


# ============================================================
# Parser Class
# ============================================================

class XtbCosmoParser:
    """
    Parser for XTB/Turbomole COSMO (.cosmo) files.
    Includes optional sigma-profile analysis methods.
    """

    def __init__(self, filepath: str):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"COSMO file not found: {filepath}")

        self.filepath = filepath

        # Internal storage (mirrors SigmaProfileParser structure)
        self._data = {
            "method": "",
            "area": None,
            "volume": None,
            "energy_tot": None,
            "energy_dielectric": None,
            "atm_nr": [],
            "atm_pos": [],
            "atm_elmnt": [],
            "atm_rad": [],
            "seg_nr": [],
            "seg_atm_nr": [],
            "seg_pos": [],
            "seg_charge": [],
            "seg_area": [],
            "seg_sigma_raw": [],
            "seg_potential": [],
        }

        self._read_turbomolesp()

    # ============================================================
    # Public API
    # ============================================================

    def get_cosmo_data(self) -> XtbCosmoData:
        """Return parsed COSMO data as a structured dataclass."""

        atoms = [
            CosmoAtom(
                index=int(i),
                element=str(el),
                position=np.array(pos, dtype=float),
                radius=float(rad) if rad is not None else None
            )
            for i, el, pos, rad in zip(
                self._data["atm_nr"],
                self._data["atm_elmnt"],
                self._data["atm_pos"],
                self._data["atm_rad"],
            )
        ]

        segments = [
            CosmoSegment(
                index=int(i),
                atom_index=int(a),
                position=np.array(pos, dtype=float),
                charge=float(q),
                area=float(area),
                sigma_raw=float(sigma),
                potential=float(pot),
            )
            for i, a, pos, q, area, sigma, pot in zip(
                self._data["seg_nr"],
                self._data["seg_atm_nr"],
                self._data["seg_pos"],
                self._data["seg_charge"],
                self._data["seg_area"],
                self._data["seg_sigma_raw"],
                self._data["seg_potential"],
            )
        ]

        return XtbCosmoData(
            filepath=self.filepath,
            method=self._data["method"],
            area=self._data["area"],
            volume=self._data["volume"],
            energy_tot=self._data["energy_tot"],
            energy_dielectric=self._data["energy_dielectric"],
            atoms=atoms,
            segments=segments,
        )

    # ============================================================
    # Internal Parsing Logic (from SigmaProfileParser)
    # ============================================================

    def _read_single_float(self, line, variable, regex, scaling_factor):
        match = re.match(regex, line)
        if match:
            self._data[variable] = float(match.groups()[0]) * scaling_factor

    def _read_turbomole_atom_section(self, f):
        line = True
        mode = None

        while line:
            line = next(f).strip()
            parts = line.split()

            if len(parts) not in (4, 6):
                if mode:
                    break
                continue

            if not mode:
                mode = len(parts)

            if len(parts) == 6:
                atm_nr = int(parts[0]) - 1
                atm_pos = [float(v) for v in parts[1:4]]
                atm_elmnt = parts[4].title()
                atm_rad = float(parts[5])
            else:
                atm_nr = len(self._data["atm_nr"])
                atm_pos = [float(v) for v in parts[1:4]]
                atm_elmnt = parts[0].title()
                atm_rad = None

            self._data["atm_nr"].append(atm_nr)
            self._data["atm_pos"].append(atm_pos)
            self._data["atm_elmnt"].append(atm_elmnt)
            self._data["atm_rad"].append(atm_rad)

        multiplier = 1 if mode == 4 else angstrom_per_bohr

        self._data["atm_nr"] = np.array(self._data["atm_nr"], dtype=int)
        self._data["atm_pos"] = np.array(self._data["atm_pos"], dtype=float) * multiplier
        self._data["atm_rad"] = np.array(self._data["atm_rad"], dtype=float)

    def _read_turbomole_seg_section(self, f):
        line = next(f)

        while line:
            try:
                line = next(f).strip()
            except StopIteration:
                break

            parts = line.split()
            if len(parts) != 9:
                if self._data["seg_nr"]:
                    break
                continue

            self._data["seg_nr"].append(int(parts[0]) - 1)
            self._data["seg_atm_nr"].append(int(parts[1]) - 1)
            self._data["seg_pos"].append([float(v) for v in parts[2:5]])
            self._data["seg_charge"].append(float(parts[5]))
            self._data["seg_area"].append(float(parts[6]))
            self._data["seg_sigma_raw"].append(float(parts[7]))
            self._data["seg_potential"].append(float(parts[8]))

        self._data["seg_nr"] = np.array(self._data["seg_nr"], dtype=int)
        self._data["seg_atm_nr"] = np.array(self._data["seg_atm_nr"], dtype=int)
        self._data["seg_pos"] = np.array(self._data["seg_pos"], dtype=float) * angstrom_per_bohr
        self._data["seg_charge"] = np.array(self._data["seg_charge"], dtype=float)
        self._data["seg_area"] = np.array(self._data["seg_area"], dtype=float)
        self._data["seg_sigma_raw"] = np.array(self._data["seg_sigma_raw"], dtype=float)
        self._data["seg_potential"] = (
            np.array(self._data["seg_potential"], dtype=float)
            * kJdivmol_per_hartree
            * angstrom_per_bohr
        )

    def _read_turbomolesp(self):
        with open(self.filepath, "r") as f:
            for i, line in enumerate(f):
                line = line.strip()

                if line == "$info":
                    line = next(f).strip()
                    self._data["method"] = (
                        f'{line.split(";")[-2]}_{line.split(";")[-1]}'.lower()
                    )

                self._read_single_float(
                    line, "area", r"area\s*=\s*([0-9+-.eE]+)", angstrom_per_bohr**2
                )

                self._read_single_float(
                    line, "volume", r"volume\s*=\s*([0-9+-.eE]+)", angstrom_per_bohr**3
                )

                self._read_single_float(
                    line,
                    "energy_tot",
                    r"Total\s+energy\s+corrected.*=\s*([0-9+-.eE]+)",
                    kJdivmol_per_hartree,
                )

                self._read_single_float(
                    line,
                    "energy_dielectric",
                    r"Dielectric\s+energy\s+\[a\.u\.\]\s*=\s*([0-9+-.eE]+)",
                    kJdivmol_per_hartree,
                )

                if line == "$coord_rad" or (i == 0 and line == "$coord_car"):
                    self._read_turbomole_atom_section(f)

                if line == "$segment_information":
                    self._read_turbomole_seg_section(f)

                # After parsing completes
                if len(self._data["seg_nr"]) == 0:
                    raise ValueError(
                        f"COSMO file '{self.filepath}' contains zero segments — "
                        "XTB may have failed or produced an incomplete COSMO surface."
                    )
                
                if self._data["energy_tot"] is None:
                    raise ValueError(
                        f"COSMO file '{self.filepath}' is missing total energy — "
                        "XTB may have failed silently."
                    )



    # ============================================================
    # Sigma-profile analysis (optional)
    # ============================================================

    def calculate_averaged_sigmas(self, averaging_radius=0.5):
        seg_pos = self._data["seg_pos"]
        seg_area = self._data["seg_area"]
        sigmas_raw = self._data["seg_sigma_raw"]

        seg_radii_sq = seg_area / np.pi
        r_av_sq = averaging_radius**2

        sigmas_avg = np.zeros_like(sigmas_raw)

        for i in range(len(sigmas_raw)):
            d2 = np.sum((seg_pos - seg_pos[i])**2, axis=1)
            denom = seg_radii_sq + r_av_sq
            weights = (seg_radii_sq * r_av_sq / denom) * np.exp(-d2 / denom)
            sigmas_avg[i] = np.sum(sigmas_raw * weights) / np.sum(weights)

        self._data["seg_sigma_averaged"] = sigmas_avg
        return sigmas_avg

    def calculate_sigma_moments(self, sigmas=None):
        if sigmas is None:
            sigmas = self._data.get("seg_sigma_averaged")
            if sigmas is None:
                sigmas = self.calculate_averaged_sigmas()

        areas = self._data["seg_area"]
        n_moments = 7
        moments = np.zeros(n_moments)

        for i in range(n_moments):
            moments[i] = np.sum((sigmas**i) * areas)
            if i > 1:
                moments[i] *= 100**i

        self._data["sigma_moments"] = moments
        return moments

    def cluster_segments_into_segmenttypes(self, descriptors, descriptor_ranges):
        # Direct port of SPP logic (unchanged)
        # Included for completeness
        areas = self._data["seg_area"]
        descriptors = np.transpose(np.array(descriptors))
        n_desc = descriptors.shape[1]

        def cluster_one(seg_desc, seg_area, idx):
            seg_types = []
            seg_type_areas = []
            drange = descriptor_ranges[idx]

            for d, area in zip(seg_desc, seg_area):
                val = d[idx]
                ind_left = np.flatnonzero(drange <= val)[-1]
                left_val = drange[ind_left]

                left_desc = d.copy()
                left_desc[idx] = left_val
                seg_types.append(left_desc)

                if left_val == val:
                    seg_type_areas.append(area)
                else:
                    right_val = drange[ind_left + 1]
                    delta = right_val - left_val

                    right_desc = d.copy()
                    right_desc[idx] = right_val

                    seg_types.append(right_desc)
                    seg_type_areas.append(area * ((right_val - val) / delta))
                    seg_type_areas.append(area * ((val - left_val) / delta))

            return seg_types, seg_type_areas

        all_types = []
        all_areas = []

        for d, area in zip(descriptors, areas):
            types, type_areas = [d.tolist()], [area]

            for i in range(n_desc):
                types, type_areas = cluster_one(types, type_areas, i)

            for t, a in zip(types, type_areas):
                if t not in all_types:
                    all_types.append(t)
                    all_areas.append(0)
                idx = all_types.index(t)
                all_areas[idx] += a

        all_types = np.array(all_types)
        all_areas = np.array(all_areas)

        # Sort
        descs = [all_types[:, i].tolist() for i in range(n_desc)]
        descs.reverse()
        idx = np.lexsort(descs)

        return all_types[idx], all_areas[idx]

    def cluster_and_create_sigma_profile(self, sigmas="seg_sigma_averaged",
                                         sigmas_range=np.arange(-0.03, 0.03, 0.001)):

        if sigmas == "seg_sigma_averaged":
            if "seg_sigma_averaged" not in self._data:
                self.calculate_averaged_sigmas()
            sigmas = self._data["seg_sigma_averaged"]

        descriptors = [sigmas]
        descriptor_ranges = [sigmas_range]

        clustered_sigmas, clustered_areas = self.cluster_segments_into_segmenttypes(
            descriptors, descriptor_ranges
        )

        clustered_sigmas = clustered_sigmas.reshape(-1)

        sp_areas = []
        sp_sigmas = []

        for sigma in sigmas_range:
            area = 0.0
            if sigma in clustered_sigmas:
                area = clustered_areas[clustered_sigmas == sigma][0]
            sp_areas.append(area)
            sp_sigmas.append(sigma)

        return np.array(sp_sigmas), np.array(sp_areas)
