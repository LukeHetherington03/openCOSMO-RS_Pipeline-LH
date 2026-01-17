# modules/parsers/cosmo/orca_cpcm_corr_parser.py

import os
import re

class OrcaCpcmCorrParserError(Exception):
    pass


class OrcaCpcmCorrParser:
    """
    Parses ORCA .cpcm_corr file:
      - corrected dielectric energy
      - total C-PCM charge
      - corrected charges list
    """

    def __init__(self, filepath: str):
        if not os.path.exists(filepath):
            raise OrcaCpcmCorrParserError(f"CPCM_CORR file not found: {filepath}")
        self.filepath = filepath

    def parse(self):
        corrected_charges = []
        corrected_dielectric_energy = None
        total_cpcm_charge = None

        with open(self.filepath, "r") as f:
            for line in f:
                stripped = line.strip()

                # Corrected dielectric energy
                if stripped.startswith("Corrected dielectric energy"):
                    m = re.search(r"=\s*([-\d\.Ee+]+)", stripped)
                    if m:
                        corrected_dielectric_energy = float(m.group(1))
                    continue

                # Total C-PCM charge
                if stripped.startswith("Total C-PCM charge"):
                    m = re.search(r"=\s*([-\d\.Ee+]+)", stripped)
                    if m:
                        total_cpcm_charge = float(m.group(1))
                    continue

                # Corrected charges block
                if re.match(r"^[\-\d\.Ee+]+$", stripped):
                    try:
                        corrected_charges.append(float(stripped))
                    except ValueError:
                        pass

        return {
            "corrected_dielectric_energy": corrected_dielectric_energy,
            "total_cpcm_charge": total_cpcm_charge,
            "corrected_charges": corrected_charges,
        }
