# modules/parsers/cosmo/orca_log_parser.py

import os
import re

class OrcaLogParserError(Exception):
    pass


class OrcaLogParser:
    """
    Extracts:
      - final single point energy
      - dipole moment (Debye)
      - atomic polarizabilities
      - SCF convergence info (future use)
      - method/basis metadata (future use)
      - timing information (future use)
    """

    def __init__(self, filepath: str):
        if not os.path.exists(filepath):
            raise OrcaLogParserError(f"Log file not found: {filepath}")
        self.filepath = filepath

    def parse(self):
        energy = None
        dipole = None
        polarizabilities = []

        method = None
        basis = None
        scf_converged = False
        total_time = None

        in_pol_block = False

        with open(self.filepath, "r") as f:
            for line in f:

                # Normal termination
                if "****ORCA TERMINATED NORMALLY****" in line:
                    scf_converged = True

                # Energy
                if "FINAL SINGLE POINT ENERGY" in line:
                    try:
                        energy = float(line.split()[-1])
                    except:
                        pass

                # Dipole
                if line.strip().startswith("x,y,z [Debye]:"):
                    parts = line.strip().split()
                    dipole = list(map(float, parts[-3:]))

                # Method / basis
                if "Method" in line and ":" in line:
                    method = line.split(":")[1].strip()
                if "Basis" in line and ":" in line:
                    basis = line.split(":")[1].strip()

                # Atomic polarizabilities block
                if re.match(r"\s*XX\s+YY\s+ZZ\s+XY\s+XZ\s+YZ", line):
                    in_pol_block = True
                    polarizabilities.append(line.rstrip())
                    continue

                if in_pol_block:
                    polarizabilities.append(line.rstrip())
                    if "Sum polar" in line:
                        in_pol_block = False

                # Timing
                if "TOTAL RUN TIME" in line:
                    total_time = line.strip()

        return {
            "energy": energy,
            "dipole": dipole,
            "polarizabilities": polarizabilities,
            "method": method,
            "basis": basis,
            "scf_converged": scf_converged,
            "total_time": total_time,
        }
