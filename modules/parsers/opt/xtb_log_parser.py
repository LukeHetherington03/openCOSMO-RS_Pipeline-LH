import re
from typing import Optional, Dict


class XTBLogParser:
    """
    Placeholder parser for classic XTB log files.

    Extracts (best-effort):
      - final energy (TOTAL ENERGY)
      - convergence flag (heuristic)
      - iterations (not printed by XTB → None)
      - elapsed_seconds (not printed by XTB → None)
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
            "energy": XTBLogParser._parse_energy(text),
            "iterations": None,  # XTB does not report this
            "converged": XTBLogParser._parse_convergence(text),
            "elapsed_seconds": None,  # XTB does not report timestamps
        }

    # ------------------------------------------------------------------
    # ENERGY
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_energy(text: str) -> Optional[float]:
        """
        Extract TOTAL ENERGY from classic XTB logs.

        Typical line:
            TOTAL ENERGY       -123.456789 Eh
        """

        m = re.search(r"TOTAL ENERGY\s+([-\d\.Ee]+)", text)
        if m:
            try:
                return float(m.group(1))
            except Exception:
                pass

        return None

    # ------------------------------------------------------------------
    # CONVERGENCE
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_convergence(text: str) -> Optional[bool]:
        """
        XTB prints:
            *** GEOMETRY OPTIMIZATION CONVERGED ***
        """

        if re.search(r"GEOMETRY OPTIMIZATION CONVERGED", text, re.IGNORECASE):
            return True

        # Hard failure patterns
        if re.search(r"forrtl:|severe|error|failed", text, re.IGNORECASE):
            return False

        return None
