#!/usr/bin/env python3
"""
Robust ORCA Log File Parser

Structured parser with explicit block handlers:
  - Energy
  - Dipole
  - Atomic polarizabilities
  - Method / basis
  - Runtime
"""

import os
import re


class OrcaLogParserError(Exception):
    pass


class OrcaLogParser:
    """
    Fully structured ORCA log parser.
    """

    def __init__(self, filepath: str):
        if not os.path.exists(filepath):
            raise OrcaLogParserError(f"Log file not found: {filepath}")
        self.filepath = filepath

        # Regex patterns
        self.float_re = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
        self.pol_line_re = re.compile(
            r"""
            ^\s*(\d+)-([A-Za-z]+)\s*:\s*              # atom index + element
            ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+      # XX
            ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+      # YY
            ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+      # ZZ
            ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+      # XY
            ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+      # XZ
            ([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+      # YZ
            iso=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)  # iso
            """,
            re.VERBOSE
        )

    # ------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------
    def parse(self):
        lines = self._read_lines()

        return {
            "final_single_point_energy": self._parse_energy_block(lines),
            "dipole_moment_debye": self._parse_dipole_magnitude(lines),
            "dipole_components_debye": self._parse_dipole_components(lines),
            "scf_converged": self._parse_scf_convergence(lines),
            "terminated_normally": self._parse_normal_termination(lines),
            "atomic_polarizabilities": self._parse_atomic_polarizabilities(lines),
            "sum_polarizability": self._parse_sum_polarizability(lines),
            "method": self._parse_method(lines),
            "basis_set": self._parse_basis(lines),
            "total_run_time": self._parse_runtime(lines),
        }

    # ------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------
    def _read_lines(self):
        with open(self.filepath, "r") as f:
            return [line.rstrip("\n") for line in f]

    # ------------------------------------------------------------
    # Energy
    # ------------------------------------------------------------
    def _parse_energy_block(self, lines):
        for line in lines:
            if "FINAL SINGLE POINT ENERGY" in line:
                m = self.float_re.search(line)
                if m:
                    return float(m.group(0))
        return None

    # ------------------------------------------------------------
    # Dipole magnitude
    # ------------------------------------------------------------
    def _parse_dipole_magnitude(self, lines):
        for line in lines:
            if line.strip().startswith("Magnitude (Debye)"):
                m = self.float_re.search(line)
                if m:
                    return float(m.group(0))
        return None

    # ------------------------------------------------------------
    # Dipole components
    # ------------------------------------------------------------
    def _parse_dipole_components(self, lines):
        for line in lines:
            if line.strip().startswith("x,y,z [Debye]:"):
                nums = self.float_re.findall(line)
                if len(nums) >= 3:
                    return [float(nums[-3]), float(nums[-2]), float(nums[-1])]
        return None

    # ------------------------------------------------------------
    # SCF convergence / termination
    # ------------------------------------------------------------
    def _parse_scf_convergence(self, lines):
        return any("****ORCA TERMINATED NORMALLY****" in line for line in lines)

    def _parse_normal_termination(self, lines):
        return any("****ORCA TERMINATED NORMALLY****" in line for line in lines)

    # ------------------------------------------------------------
    # Atomic polarizabilities
    # ------------------------------------------------------------
    def _parse_atomic_polarizabilities(self, lines):
        pols = []
        in_block = False

        for line in lines:
            s = line.strip()

            if "ATOMIC POLARIZABILITIES" in s:
                in_block = True
                continue

            if in_block:
                # Skip separators, headers, notes
                if not s:
                    continue
                if s.startswith('-'):
                    continue
                if "NOTE" in s:
                    continue
                if "XX" in s and "YY" in s and "ZZ" in s:
                    continue

                # End when we reach the sum line
                if s.startswith("Sum polar"):
                    break

                m = self.pol_line_re.match(s)
                if m:
                    atom_index = int(m.group(1))
                    element = m.group(2)
                    XX = float(m.group(3))
                    YY = float(m.group(4))
                    ZZ = float(m.group(5))
                    XY = float(m.group(6))
                    XZ = float(m.group(7))
                    YZ = float(m.group(8))
                    iso = float(m.group(9))

                    pols.append({
                        "atom_index": atom_index,
                        "element": element,
                        "XX": XX,
                        "YY": YY,
                        "ZZ": ZZ,
                        "XY": XY,
                        "XZ": XZ,
                        "YZ": YZ,
                        "iso": iso,
                    })

        return pols

    # ------------------------------------------------------------
    # Sum polarizability
    # ------------------------------------------------------------
    def _parse_sum_polarizability(self, lines):
        for line in lines:
            s = line.strip()
            if s.startswith("Sum polar"):
                nums = self.float_re.findall(s)
                if len(nums) >= 7:
                    return {
                        "XX": float(nums[0]),
                        "YY": float(nums[1]),
                        "ZZ": float(nums[2]),
                        "XY": float(nums[3]),
                        "XZ": float(nums[4]),
                        "YZ": float(nums[5]),
                        "iso": float(nums[6]),
                    }
        return None

    # ------------------------------------------------------------
    # Method / basis / runtime
    # ------------------------------------------------------------
    def _parse_method(self, lines):
        for line in lines:
            if "Method" in line and ":" in line:
                return line.split(":", 1)[1].strip()
        return None

    def _parse_basis(self, lines):
        for line in lines:
            if "Basis" in line and ":" in line:
                return line.split(":", 1)[1].strip()
        return None

    def _parse_runtime(self, lines):
        for line in lines:
            if "TOTAL RUN TIME" in line.upper():
                m = self.float_re.search(line)
                if m:
                    return float(m.group(0))
        return None


def parse(filepath: str):
    return OrcaLogParser(filepath).parse()


def parse_log(filepath: str):
    return parse(filepath)


if __name__ == "__main__":
    import sys, json
    fp = sys.argv[1]
    print(json.dumps(parse(fp), indent=2))
