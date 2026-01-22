#!/usr/bin/env python3
"""
Robust ORCA CPCM_CORR Parser

This parser handles:
  - Corrected dielectric energy
  - Total C-PCM charge
  - Corrected charges array
"""

import os
import re


class OrcaCpcmCorrParserError(Exception):
    pass


class OrcaCpcmCorrParser:
    """
    Fully structured CPCM_CORR parser with explicit block handlers:
      - corrected dielectric energy
      - total C-PCM charge
      - corrected charges list
    """

    def __init__(self, filepath: str):
        if not os.path.exists(filepath):
            raise OrcaCpcmCorrParserError(f"CPCM_CORR file not found: {filepath}")
        self.filepath = filepath

        # Regex for floating point numbers (robust)
        self.float_re = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")

    # ------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------
    def parse(self):
        lines = self._read_lines()

        corrected_energy = self._parse_corrected_energy_block(lines)
        total_charge = self._parse_total_charge_block(lines)
        corrected_charges = self._parse_corrected_charges_block(lines)

        return {
            "corrected_dielectric_energy": corrected_energy,
            "total_cpcm_charge": total_charge,
            "corrected_charges": corrected_charges,
        }

    # ------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------
    def _read_lines(self):
        with open(self.filepath, "r") as f:
            return [line.rstrip("\n") for line in f]

    # ------------------------------------------------------------
    # Corrected dielectric energy
    # ------------------------------------------------------------
    def _parse_corrected_energy_block(self, lines):
        for line in lines:
            if "Corrected dielectric energy" in line:
                m = self.float_re.search(line)
                if m:
                    return float(m.group(0))
        return None

    # ------------------------------------------------------------
    # Total C-PCM charge
    # ------------------------------------------------------------
    def _parse_total_charge_block(self, lines):
        for line in lines:
            if "Total C-PCM charge" in line:
                m = self.float_re.search(line)
                if m:
                    return float(m.group(0))
        return None

    # ------------------------------------------------------------
    # Corrected charges block
    # ------------------------------------------------------------
    def _parse_corrected_charges_block(self, lines):
        corrected_charges = []
        in_block = False

        for line in lines:
            s = line.strip()

            if "C-PCM corrected charges:" in s:
                in_block = True
                continue

            if in_block:
                # Stop if we hit a non-numeric line
                if not s or not self.float_re.fullmatch(s):
                    continue

                try:
                    corrected_charges.append(float(s))
                except ValueError:
                    continue

        return corrected_charges


# ------------------------------------------------------------
# Convenience API
# ------------------------------------------------------------
def parse(filepath: str):
    return OrcaCpcmCorrParser(filepath).parse()


def parse_cpcm_corr(filepath: str):
    return parse(filepath)


if __name__ == "__main__":
    import sys
    import json

    if len(sys.argv) != 2:
        print("Usage: python orca_cpcm_corr_parser.py <cpcm_corr_file>")
        sys.exit(1)

    fp = sys.argv[1]
    result = parse(fp)
    print(json.dumps(result, indent=2))
