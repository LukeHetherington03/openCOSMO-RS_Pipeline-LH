import os
import re

class OrcaCpcmParserError(Exception):
    pass


class OrcaCpcmParser:
    """
    Robust parser for ORCA .cpcm files.
    Handles:
      - metadata block (all numeric header values)
      - cartesian coordinates + radii + Z
      - surface points block
      - ignores comments (# ...)
      - tolerant to ORCA version differences
    """

    def __init__(self, filepath: str):
        if not os.path.exists(filepath):
            raise OrcaCpcmParserError(f"CPCM file not found: {filepath}")
        self.filepath = filepath

    def _strip_comment(self, line: str) -> str:
        """Remove trailing ORCA comments (# ...)."""
        return line.split("#")[0].strip()

    def parse(self):
        metadata = {}
        cartesian = []
        surface_points = []

        mode = None  # None, "cartesian", "surface"

        with open(self.filepath, "r") as f:
            raw_lines = f.readlines()

        # Preprocess: strip comments, keep raw for debugging
        lines = [self._strip_comment(l) for l in raw_lines]

        i = 0
        n = len(lines)

        # ------------------------------------------------------------
        # 1. Parse metadata block (first long numeric block)
        # ------------------------------------------------------------
        numeric_values = []

        while i < n:
            line = lines[i].strip()

            # Stop when we hit the Cartesian header
            if "CARTESIAN COORDINATES" in line:
                break

            # Collect numeric lines
            if re.match(r"^[\d\.\-Ee+]+$", line):
                numeric_values.append(line)

            i += 1

        # Now interpret numeric_values in order
        # ORCA CPCM metadata always appears in the same order
        try:
            idx = 0
            metadata["n_atoms"] = int(numeric_values[idx]); idx += 1
            metadata["n_surface_points"] = int(numeric_values[idx]); idx += 1
            metadata["surface_type"] = int(numeric_values[idx]); idx += 1
            metadata["epsilon_function_type"] = int(numeric_values[idx]); idx += 1
            metadata["print_level"] = int(numeric_values[idx]); idx += 1
            metadata["feps_x_flag"] = int(numeric_values[idx]); idx += 1
            metadata["draco_scheme"] = int(numeric_values[idx]); idx += 1

            metadata["n_lebedev_points"] = int(numeric_values[idx]); idx += 1
            metadata["isodensity_scheme"] = int(numeric_values[idx]); idx += 1

            metadata["threshold_h"] = float(numeric_values[idx]); idx += 1
            metadata["threshold_non_h"] = float(numeric_values[idx]); idx += 1

            metadata["feps_x_parameter"] = float(numeric_values[idx]); idx += 1

            metadata["cutoff_segment_area"] = float(numeric_values[idx]); idx += 1
            metadata["cutoff_switching_function"] = float(numeric_values[idx]); idx += 1

            metadata["volume"] = float(numeric_values[idx]); idx += 1
            metadata["area"] = float(numeric_values[idx]); idx += 1

            metadata["dielectric_energy"] = float(numeric_values[idx]); idx += 1
            metadata["one_electron_energy"] = float(numeric_values[idx]); idx += 1

        except Exception as e:
            raise OrcaCpcmParserError(f"Failed to parse CPCM metadata: {e}")

        # ------------------------------------------------------------
        # 2. Parse Cartesian block
        # ------------------------------------------------------------
        while i < n:
            line = lines[i].strip()

            if "CARTESIAN COORDINATES" in line:
                mode = "cartesian"
                i += 1
                continue

            if mode == "cartesian":
                if not line or line.startswith("#"):
                    i += 1
                    continue

                parts = line.split()
                if len(parts) == 5:
                    x, y, z, r, Z = parts
                    cartesian.append([
                        float(x), float(y), float(z),
                        float(r), int(Z)
                    ])
                    i += 1
                    continue
                else:
                    mode = None

            if "SURFACE POINTS" in line:
                break

            i += 1

        # ------------------------------------------------------------
        # 3. Parse Surface Points block
        # ------------------------------------------------------------
        while i < n:
            line = lines[i].strip()

            if "SURFACE POINTS" in line:
                mode = "surface"
                i += 1
                continue

            if mode == "surface":
                if not line or line.startswith("#"):
                    i += 1
                    continue

                parts = line.split()
                if len(parts) >= 10:
                    surface_points.append([
                        float(parts[0]), float(parts[1]), float(parts[2]),
                        float(parts[3]), float(parts[4]), float(parts[5]),
                        float(parts[6]), float(parts[7]), float(parts[8]),
                        int(parts[9])
                    ])
                    i += 1
                    continue
                else:
                    mode = None

            i += 1

        return {
            "metadata": metadata,
            "cartesian": cartesian,
            "surface_points": surface_points,
        }
