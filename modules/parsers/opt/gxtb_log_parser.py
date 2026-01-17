# modules/parsers/opt/gxtb_log_parser.py

import re
from datetime import datetime
from typing import Optional, Dict


class GxTBLogParser:
    """
    Robust parser for gXTB log files.

    Extracts:
      - final energy (TOTAL ENERGY)
      - number of optimisation iterations
      - convergence status
      - elapsed time (seconds)
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
            "energy": GxTBLogParser._parse_energy(text),
            "iterations": GxTBLogParser._parse_iterations(text),
            "converged": GxTBLogParser._parse_convergence(text),
            "elapsed_seconds": GxTBLogParser._parse_elapsed(text),
        }

    # ------------------------------------------------------------------
    # ENERGY
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_energy(text: str) -> Optional[float]:
        """
        Extract ONLY the final TOTAL ENERGY (uppercase).
        Ignore earlier SCF iteration energies.
        """

        # Pattern: TOTAL ENERGY   -306.46305918
        m = re.search(r"TOTAL ENERGY\s+([-\d\.Ee]+)", text)
        if m:
            try:
                return float(m.group(1))
            except Exception:
                pass

        # Pattern: "total                         -306.46305918"
        m = re.search(r"\btotal\s+([-\d\.Ee]+)", text)
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
        Extract from:
            *** GEOMETRY OPTIMIZATION CONVERGED AFTER 6 ITERATIONS ***
        """

        m = re.search(
            r"GEOMETRY OPTIMIZATION CONVERGED AFTER\s+(\d+)\s+ITERATIONS",
            text,
            re.IGNORECASE,
        )
        if m:
            return int(m.group(1))

        return None

    # ------------------------------------------------------------------
    # CONVERGENCE
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_convergence(text: str) -> Optional[bool]:
        """
        Converged if iterations line exists.
        Failed if forrtl or severe errors appear.
        """

        # Converged
        if re.search(r"GEOMETRY OPTIMIZATION CONVERGED", text, re.IGNORECASE):
            return True

        # Hard failure patterns
        if re.search(r"forrtl:|severe|error", text, re.IGNORECASE):
            return False

        return None

    # ------------------------------------------------------------------
    # ELAPSED TIME
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_elapsed(text: str) -> Optional[float]:
        """
        Extract start and finish timestamps and compute elapsed seconds.

        Example:
            * started run on 2026/01/14 at 16:13:26.027
            * finished run on 2026/01/14 at 16:13:32.199
        """

        start = re.search(
            r"started run on\s+(\d{4}/\d{2}/\d{2})\s+at\s+([\d:.]+)",
            text,
        )
        end = re.search(
            r"finished run on\s+(\d{4}/\d{2}/\d{2})\s+at\s+([\d:.]+)",
            text,
        )

        if not start or not end:
            return None

        try:
            start_dt = datetime.strptime(
                start.group(1) + " " + start.group(2),
                "%Y/%m/%d %H:%M:%S.%f",
            )
            end_dt = datetime.strptime(
                end.group(1) + " " + end.group(2),
                "%Y/%m/%d %H:%M:%S.%f",
            )
            return (end_dt - start_dt).total_seconds()
        except Exception:
            return None
