import re
from typing import Optional, Dict


class ForcefieldLogParser:
    """
    Parser for RDKit forcefield (MMFF94 / UFF) optimisation logs.

    Expected log format written by _backend_forcefield:

        FORCEFIELD_FAMILY  MMFF94
        OPTIMIZATION_STATUS  converged
        FORCEFIELD_ENERGY  -42.123456
        ATOMS  32

    Energy is in kcal/mol (RDKit native units).
    """

    @staticmethod
    def parse(log_path: str) -> Dict[str, Optional[float]]:
        try:
            with open(log_path, "r", errors="ignore") as f:
                text = f.read()
        except Exception:
            return {"energy": None, "converged": None, "iterations": None, "elapsed_seconds": None}

        return {
            "energy":          ForcefieldLogParser._parse_energy(text),
            "converged":       ForcefieldLogParser._parse_convergence(text),
            "iterations":      None,
            "elapsed_seconds": None,
        }

    @staticmethod
    def _parse_energy(text: str) -> Optional[float]:
        m = re.search(r"^FORCEFIELD_ENERGY\s+([-\d.Ee]+)", text, re.MULTILINE)
        if m:
            try:
                return float(m.group(1))
            except Exception:
                pass
        return None

    @staticmethod
    def _parse_convergence(text: str) -> Optional[bool]:
        m = re.search(r"^OPTIMIZATION_STATUS\s+(\S+)", text, re.MULTILINE)
        if m:
            status = m.group(1).lower()
            if status == "converged":
                return True
            if status == "failed":
                return False
        return None
