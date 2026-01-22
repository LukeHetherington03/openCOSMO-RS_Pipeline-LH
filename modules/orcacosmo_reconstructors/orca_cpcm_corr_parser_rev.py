class OrcaCpcmCorrParserRev:
    """
    Reconstructs an ORCA-style CPCM_CORR block from parsed JSON.

    Expected JSON structure:
    {
        "corrected_dielectric_energy": float,
        "total_cpcm_charge": float,
        "corrected_charges": [float, float, ...]
    }
    """

    def __init__(self, cpcm_corr_json):
        self.data = cpcm_corr_json

    # ------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------
    def parse(self):
        """Return parsed CPCM_CORR JSON (for API symmetry)."""
        return self.data

    def reconstruct(self):
        """Return the full CPCM_CORR block exactly as ORCA prints it."""
        out = []

        # Header
        out.append("##################################################")
        out.append("#COSMO_corrected")

        # Corrected dielectric energy
        out.append(
            f"Corrected dielectric energy   = {self._fmt(self.data['corrected_dielectric_energy'])}"
        )

        # Total C-PCM charge
        out.append(
            f"Total C-PCM charge            = {self._fmt(self.data['total_cpcm_charge'])}"
        )

        # Charges header
        out.append("C-PCM corrected charges:")

        # Charges list (1 per line, indented 4 spaces)
        for charge in self.data.get("corrected_charges", []):
            out.append(f"    {self._fmt(charge)}")

        return "\n".join(out) + "\n"

    # ------------------------------------------------------------
    # Formatting helper
    # ------------------------------------------------------------
    def _fmt(self, value):
        """ORCA uses fixed-width 16.9f for CPCM_CORR floats."""
        return f"{value:16.9f}"
