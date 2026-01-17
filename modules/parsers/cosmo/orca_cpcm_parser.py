# modules/parsers/cosmo/orca_cpcm_parser.py

import os

class OrcaCpcmParserError(Exception):
    pass


class OrcaCpcmParser:
    """
    Parses ORCA .cpcm file:
      - COSMO metadata
      - cavity volume, area
      - dielectric energies
      - Lebedev grid info
      - isodensity thresholds
      - cartesian coordinates + radii + Z
      - surface points
    """

    def __init__(self, filepath: str):
        if not os.path.exists(filepath):
            raise OrcaCpcmParserError(f"CPCM file not found: {filepath}")
        self.filepath = filepath

    def parse(self):
        metadata = {}
        cartesian = []
        surface_points = []

        mode = None  # "cartesian" or "surface"

        with open(self.filepath, "r") as f:
            lines = [l.rstrip() for l in f]

        i = 0
        while i < len(lines):
            line = lines[i].strip()

            # Metadata (fixed order in ORCA)
            if line.isdigit():
                # number of atoms
                metadata["n_atoms"] = int(line)
                metadata["n_surface_points"] = int(lines[i+1])
                metadata["surface_type"] = int(lines[i+2])
                metadata["epsilon_function_type"] = int(lines[i+3])
                metadata["print_level"] = int(lines[i+4])
                metadata["feps_x_flag"] = int(lines[i+5])
                metadata["n_lebedev_points"] = int(lines[i+6])
                metadata["isodensity_scheme"] = int(lines[i+7])
                metadata["threshold_h"] = float(lines[i+8])
                metadata["threshold_non_h"] = float(lines[i+9])
                metadata["feps_x_parameter"] = float(lines[i+10])
                metadata["cutoff_segment_area"] = float(lines[i+11])
                metadata["cutoff_switching_function"] = float(lines[i+12])
                metadata["volume"] = float(lines[i+13])
                metadata["area"] = float(lines[i+14])
                metadata["dielectric_energy"] = float(lines[i+15])
                metadata["one_electron_energy"] = float(lines[i+16])
                i += 17
                continue

            # Cartesian block
            if "CARTESIAN COORDINATES" in line:
                mode = "cartesian"
                i += 1
                continue

            if mode == "cartesian" and line and not line.startswith("#"):
                parts = line.split()
                if len(parts) == 5:
                    x, y, z, r, Z = parts
                    cartesian.append([float(x), float(y), float(z), float(r), int(Z)])
                else:
                    mode = None

            # Surface points
            if "SURFACE POINTS" in line:
                mode = "surface"
                i += 1
                continue

            if mode == "surface" and line and not line.startswith("#"):
                parts = line.split()
                if len(parts) >= 10:
                    surface_points.append([
                        float(parts[0]), float(parts[1]), float(parts[2]),
                        float(parts[3]), float(parts[4]), float(parts[5]),
                        float(parts[6]), float(parts[7]), float(parts[8]),
                        int(parts[9])
                    ])
                else:
                    mode = None

            i += 1

        return {
            "metadata": metadata,
            "cartesian": cartesian,
            "surface_points": surface_points,
        }
