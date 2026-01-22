class OrcaLogParserRev:
    """
    Reconstructs ORCA-style log-derived blocks from parsed JSON.

    Input JSON structure:
    {
        "final_single_point_energy": float,
        "dipole_components_debye": [dx, dy, dz],
        "atomic_polarizabilities": [...],
        "sum_polarizability": {...}
    }
    """

    def __init__(self, log_json):
        self.data = log_json

    # ------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------
    def parse(self):
        """Return the parsed JSON (for orchestrator convenience)."""
        return self.data

    def reconstruct(self):
        """Return the full LOG block in canonical ORCA order."""
        out = []
        out.extend(self._energy_block())
        out.extend(self._dipole_block())
        out.extend(self._polarizability_block())
        return "\n".join(out)

    # ------------------------------------------------------------
    # Individual block APIs
    # ------------------------------------------------------------
    def reconstruct_energy(self):
        return "\n".join(self._energy_block()) + "\n"

    def reconstruct_dipole(self):
        return "\n".join(self._dipole_block()) + "\n"

    def reconstruct_polarizabilities(self):
        return "\n".join(self._polarizability_block()) + "\n"

    # ------------------------------------------------------------
    # ENERGY
    # ------------------------------------------------------------
    def _energy_block(self):
        e = self.data["final_single_point_energy"]
        return [
            "##################################################",
            "#ENERGY",
            f"FINAL SINGLE POINT ENERGY      {e}"
        ]

    # ------------------------------------------------------------
    # DIPOLE
    # ------------------------------------------------------------
    def _dipole_block(self):
        dx, dy, dz = self.data["dipole_components_debye"]
        return [
            "##################################################",
            "#DIPOLE MOMENT (Debye)",
            f"{dx} {dy} {dz}"
        ]

    # ------------------------------------------------------------
    # POLARIZABILITIES
    # ------------------------------------------------------------
    def _polarizability_block(self):
        out = []
        pols = self.data.get("atomic_polarizabilities", [])
        sum_pol = self.data.get("sum_polarizability")

        out.append("")
        out.append("##################################################")
        out.append("#ATOMIC POLARIZABILITIES")
        out.append("                                       XX                YY               ZZ              XY              XZ              YZ  ")
        out.append("                                  ------------------------------------------------------------------------------------------------")

        for p in pols:
            out.append(
                f"{p['atom_index']:3d}-{p['element']:<2}                     :"
                f"{self._fmt(p['XX'])}{self._fmt(p['YY'])}{self._fmt(p['ZZ'])}"
                f"{self._fmt(p['XY'])}{self._fmt(p['XZ'])}{self._fmt(p['YZ'])}"
                f" iso={self._fmt_iso(p['iso'])}"
            )

        out.append("                                  ------------------------------------------------------------------------------------------------")

        if sum_pol:
            out.append(
                "Sum polar. (a.u.)         :"
                f"{self._fmt(sum_pol['XX'])}{self._fmt(sum_pol['YY'])}{self._fmt(sum_pol['ZZ'])}"
                f"{self._fmt(sum_pol['XY'])}{self._fmt(sum_pol['XZ'])}{self._fmt(sum_pol['YZ'])}"
                f" iso={self._fmt_iso(sum_pol['iso'])}"
            )
        else:
            out.append("# No sum polarizability found")

        return out

    # ------------------------------------------------------------
    # Formatting helpers
    # ------------------------------------------------------------
    def _fmt(self, value):
        return f"{value:16.9f}"

    def _fmt_iso(self, value):
        return f"{value:14.9f}"
