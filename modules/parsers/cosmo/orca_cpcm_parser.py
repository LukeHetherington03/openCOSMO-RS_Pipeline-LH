#!/usr/bin/env python3
"""
Robust ORCA CPCM Parser

Parses:
  - Metadata block
  - Cartesian coordinates block
  - Surface points block
"""

import os
import re


class OrcaCpcmParserError(Exception):
    pass


class OrcaCpcmParser:
    """
    Fully structured CPCM parser with explicit block handlers:
      - metadata
      - cartesian coordinates
      - surface points
    """

    def __init__(self, filepath: str):
        if not os.path.exists(filepath):
            raise OrcaCpcmParserError(f"CPCM file not found: {filepath}")
        self.filepath = filepath

        # Regex for Cartesian lines: X Y Z radius Z_atomic
        self.cartesian_re = re.compile(
            r"""
            ^\s*
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # X
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # Y
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # Z
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # radius
            (\d+)                                  # atomic number
            \s*$
            """,
            re.VERBOSE,
        )

        # Regex for surface points (handles wrapped floats, variable spacing)
        self.surface_point_re = re.compile(
            r"""
            ^\s*
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # X
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # Y
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # Z
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # area
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # potential
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # charge
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # w_leb
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # Switch_F
            ([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+   # G_width
            (\d+)                                  # atom index
            \s*$
            """,
            re.VERBOSE,
        )

    # ------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------
    def parse(self):
        lines = self._read_lines()

        metadata = self._parse_metadata_block(lines)
        cartesian = self._parse_cartesian_block(lines)
        surface_points = self._parse_surface_block(lines)

        return {
            "metadata": metadata,
            "cartesian": cartesian,
            "surface_points": surface_points,
        }

    # ------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------
    def _read_lines(self):
        with open(self.filepath, "r") as f:
            return [line.rstrip("\n") for line in f]

    # ------------------------------------------------------------
    # METADATA BLOCK
    # ------------------------------------------------------------
    def _parse_metadata_block(self, lines):
        metadata = {}

        def extract_number(line):
            parts = line.strip().split()
            for p in parts:
                try:
                    if any(ch in p for ch in (".", "E", "e")):
                        return float(p)
                    return int(p)
                except ValueError:
                    continue
            return None

        for line in lines:
            s = line.strip()

            if "Number of atoms" in s:
                metadata["n_atoms"] = extract_number(s)

            elif "Number of surface points" in s:
                metadata["n_surface_points"] = extract_number(s)

            elif "Surface type" in s:
                metadata["surface_type"] = extract_number(s)

            elif "Epsilon function type" in s:
                metadata["epsilon_function_type"] = extract_number(s)

            elif "Print level" in s:
                metadata["print_level"] = extract_number(s)

            elif "FEps X flag" in s:
                metadata["feps_x_flag"] = extract_number(s)

            elif "Number of Leb" in s:
                metadata["n_lebedev_points"] = extract_number(s)

            elif "Isodensity discr" in s:
                metadata["isodensity_scheme"] = extract_number(s)

            elif "Threshold for H atoms" in s:
                metadata["threshold_h"] = extract_number(s)

            elif "Threshold for non-H atoms" in s:
                metadata["threshold_non_h"] = extract_number(s)

            elif "FEps X parameter" in s:
                metadata["feps_x_parameter"] = extract_number(s)

            elif "cutoff segment area" in s:
                metadata["cutoff_segment_area"] = extract_number(s)

            elif "cutoff switching function" in s:
                metadata["cutoff_switching_function"] = extract_number(s)

            elif s.endswith("# Volume"):
                metadata["volume"] = extract_number(s)

            elif s.endswith("# Area"):
                metadata["area"] = extract_number(s)

            elif "CPCM dielectric energy" in s:
                metadata["dielectric_energy"] = extract_number(s)

            elif "One-electron operator energy" in s:
                metadata["one_electron_energy"] = extract_number(s)

        return metadata

    # ------------------------------------------------------------
    # CARTESIAN BLOCK (section-specific)
    # ------------------------------------------------------------
    def _parse_cartesian_block(self, lines):
        cartesian = []
        in_block = False

        for line in lines:
            s = line.rstrip("\n")

            # Start of block
            if "CARTESIAN COORDINATES" in s:
                in_block = True
                continue

            if not in_block:
                continue

            # End of block: blank line or next header
            if not s.strip():
                if cartesian:
                    break
                else:
                    continue
            if s.strip().startswith("#") and "CARTESIAN COORDINATES" not in s:
                # header or separator after we've started
                if cartesian:
                    break
                else:
                    continue

            # Try regex match
            m = self.cartesian_re.match(s)
            if not m:
                # If we already collected some lines and this doesn't match, assume block ended
                if cartesian:
                    break
                else:
                    continue

            x, y, z, radius, z_atomic = m.groups()
            cartesian.append([
                float(x),
                float(y),
                float(z),
                float(radius),
                int(z_atomic),
            ])

        return cartesian

    # ------------------------------------------------------------
    # SURFACE POINT BLOCK (section-specific)
    # ------------------------------------------------------------
    def _parse_surface_block(self, lines):
        surface_points = []
        in_block = False
        header_seen = False

        for line in lines:
            s = line.rstrip("\n")

            # Start of block
            if "SURFACE POINTS" in s:
                in_block = True
                header_seen = False
                continue

            if not in_block:
                continue

            # Skip header line with column names
            if not header_seen and "X" in s and "area" in s and "atom" in s:
                header_seen = True
                continue

            # End of block: blank line or next header
            if not s.strip():
                if surface_points:
                    break
                else:
                    continue
            if s.strip().startswith("#") and "SURFACE POINTS" not in s:
                if surface_points:
                    break
                else:
                    continue

            m = self.surface_point_re.match(s)
            if not m:
                # If we've already collected some points and this fails, assume block ended
                if surface_points:
                    break
                else:
                    continue

            vals = m.groups()
            surface_points.append([
                float(vals[0]),  # X
                float(vals[1]),  # Y
                float(vals[2]),  # Z
                float(vals[3]),  # area
                float(vals[4]),  # potential
                float(vals[5]),  # charge
                float(vals[6]),  # w_leb
                float(vals[7]),  # Switch_F
                float(vals[8]),  # G_width
                int(vals[9]),    # atom index
            ])

        return surface_points


# ------------------------------------------------------------
# Convenience API
# ------------------------------------------------------------
def parse(filepath: str):
    return OrcaCpcmParser(filepath).parse()


def parse_cpcm(filepath: str):
    return parse(filepath)


if __name__ == "__main__":
    import sys
    import json

    if len(sys.argv) != 2:
        print("Usage: python orca_cpcm_parser.py <cpcm_file>")
        sys.exit(1)

    fp = sys.argv[1]
    result = parse(fp)
    print(json.dumps(result, indent=2))
