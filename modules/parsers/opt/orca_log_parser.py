# modules/parsers/opt/orca_log_parser.py

import re
from typing import Optional, Dict


class ORCALogParser:
    """
    Lean parser for ORCA optimisation log files.

    Extracts:
      - final energy (FINAL SINGLE POINT ENERGY)
      - number of optimisation steps
      - convergence status
      - elapsed time (seconds, from TOTAL RUN TIME)
    """

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    @staticmethod
    def parse(log_path: str) -> Dict[str, Optional[float]]:
        try:
            with open(log_path, "r", errors="ignore") as f:
                text = f.read()
        except Exception:
            return {
                "energy": None,
                "iterations": None,
                "converged": None,
                "elapsed_seconds": None,
            }

        return {
            "energy": ORCALogParser._parse_energy(text),
            "iterations": ORCALogParser._parse_iterations(text),
            "converged": ORCALogParser._parse_convergence(text),
            "elapsed_seconds": ORCALogParser._parse_elapsed(text),
        }

    # ------------------------------------------------------------------
    # ENERGY
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_energy(text: str) -> Optional[float]:
        """
        Extract FINAL SINGLE POINT ENERGY.
        """

        m = re.search(
            r"FINAL SINGLE POINT ENERGY\s+([-\d\.Ee]+)",
            text,
        )
        if m:
            try:
                return float(m.group(1))
            except Exception:
                pass

        return None

    # ------------------------------------------------------------------
    # ITERATIONS
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_iterations(text: str) -> Optional[int]:
        """
        Extract number of geometry optimisation steps.

        ORCA prints:
            GEOMETRY OPTIMIZATION CYCLE   12
        """

        cycles = re.findall(
            r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)",
            text,
        )
        if cycles:
            try:
                return int(cycles[-1])
            except Exception:
                pass

        return None

    # ------------------------------------------------------------------
    # CONVERGENCE
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_convergence(text: str) -> Optional[bool]:
        """
        Converged if ORCA prints the HURRAY banner.
        Failed if explicit ERROR markers appear.
        """

        if re.search(r"THE OPTIMIZATION HAS CONVERGED", text, re.IGNORECASE):
            return True

        if re.search(r"ERROR|FAILED|ABORT", text, re.IGNORECASE):
            return False

        return None

    # ------------------------------------------------------------------
    # ELAPSED TIME
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_elapsed(text: str) -> Optional[float]:
        """
        Extract TOTAL RUN TIME in seconds.

        ORCA prints:
            TOTAL RUN TIME: 0 days  0 hours  10 minutes  23 seconds
        """

        m = re.search(
            r"TOTAL RUN TIME:\s*(\d+)\s*days\s*(\d+)\s*hours\s*(\d+)\s*minutes\s*(\d+)\s*seconds",
            text,
        )
        if not m:
            return None

        try:
            days = int(m.group(1))
            hours = int(m.group(2))
            minutes = int(m.group(3))
            seconds = int(m.group(4))
            return (
                days * 86400
                + hours * 3600
                + minutes * 60
                + seconds
            )
        except Exception:
            return None
