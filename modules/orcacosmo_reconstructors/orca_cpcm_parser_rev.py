class OrcaCpcmParserRev:
    """
    Reconstructs an ORCA-style .cpcm file from parsed JSON.

    Expected JSON structure:
    {
        "metadata": {...},
        "cartesian": [[x,y,z,radius,Z], ...],
        "surface_points": [[x,y,z,area,potential,charge,w_leb,switch_f,g_width,atom], ...]
    }
    """

    def __init__(self, cpcm_json):
        self.data = cpcm_json
        self.meta = cpcm_json["metadata"]
        self.cart = cpcm_json["cartesian"]
        self.surf = cpcm_json["surface_points"]

    # ------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------
    def parse(self):
        return self.data

    def reconstruct(self):
        out = []
        out.extend(self._reconstruct_metadata())
        out.extend(self._reconstruct_cartesian())
        out.extend(self._reconstruct_surface_points())
        return "\n".join(out) + "\n"

    # ------------------------------------------------------------
    # Metadata
    # ------------------------------------------------------------
    def _reconstruct_metadata(self):
        m = self.meta
        out = []

        out.append(f"{m['n_atoms']}\t\t # Number of atoms ")
        out.append(f"{m['n_surface_points']}\t\t # Number of surface points")
        out.append(f"{m['surface_type']}\t\t # Surface type")
        out.append(f"{m['epsilon_function_type']}\t\t # Epsilon function type")
        out.append(f"{m['print_level']}\t\t # Print level")
        out.append(f"{m['feps_x_flag']}\t\t # FEps X flag     \n")

        out.append(f"{m['n_lebedev_points']}\t\t # Number of Leb. points  (used if Isodens = 0)   \n")
        out.append(f"{m['isodensity_scheme']}\t\t # Isodensity discr. scheme  \n")

        out.append(f"{self._fmt_meta(m['threshold_h'])}\t # Threshold for H atoms (Ncharges/Angstroem^2) for Isodens = 1")
        out.append(f"{self._fmt_meta(m['threshold_non_h'])}\t # Threshold for non-H atoms (Ncharges/Angstroem^2) for Isodens = 1\n")

        out.append(f"{self._fmt_meta(m['feps_x_parameter'])}\t # FEps X parameter\n")

        out.append(f"{self._fmt_meta(m['cutoff_segment_area'])}\t # cutoff segment area (a.u.)")
        out.append(f"{self._fmt_meta(m['cutoff_switching_function'])}\t # cutoff switching function\n")

        out.append(f"{self._fmt_meta(m['volume'])}\t # Volume")
        out.append(f"{self._fmt_meta(m['area'])}\t # Area\n")

        out.append(f"{self._fmt_meta(m['dielectric_energy'])}\t # CPCM dielectric energy")
        out.append(f"{self._fmt_meta(m['one_electron_energy'])}\t # One-electron operator energy\n")

        return out

    # ------------------------------------------------------------
    # Cartesian coordinates (section-specific)
    # ------------------------------------------------------------
    def _reconstruct_cartesian(self):
        out = []
        out.append("#------------------------------------------------------------")
        out.append("# CARTESIAN COORDINATES (A.U.) + RADII (A.U.) + ATOMIC NUMBER ")
        out.append("#------------------------------------------------------------")

        for x, y, z, r, Z in self.cart:
            # explicit spaces between columns, 9 decimals like your input
            line = (
                f"{self._fmt_cart(x)} {self._fmt_cart(y)} {self._fmt_cart(z)} "
                f"{self._fmt_cart(r)}  {int(Z):2d}"
            )
            out.append(line)

        out.append("")  # blank line
        return out

    # ------------------------------------------------------------
    # Surface points (section-specific)
    # ------------------------------------------------------------
    def _reconstruct_surface_points(self):
        out = []
        out.append("#------------------------------------------------------------")
        out.append("# SURFACE POINTS (A.U.)    (Hint - charge NOT scaled by FEps)")
        out.append("#------------------------------------------------------------")
        out.append(
            "          X                 Y                 Z               area            potential          charge            w_leb             Switch_F          G_width       atom"
        )

        for pt in self.surf:
            x, y, z, area, pot, charge, w_leb, switch_f, g_width, atom = pt

            line = (
                f"{self._fmt_surf_coord(x)}"
                f"{self._fmt_surf_coord(y)}"
                f"{self._fmt_surf_coord(z)}"
                f"{self._fmt_surf(area)}"
                f"{self._fmt_surf(pot)}"
                f"{self._fmt_surf(charge)}"
                f"{self._fmt_wleb(w_leb)} "
                f"{self._fmt_switch(switch_f)} "
                f"{self._fmt_gwidth(g_width)}      {int(atom)}"
            )
            out.append(line)

        return out

    # ------------------------------------------------------------
    # Formatting helpers (section-specific)
    # ------------------------------------------------------------

    def _fmt_surf_coord(self, value: float) -> str:
        # Same visual style as Cartesian coords
        return f"{value:15.9f}"


    def _fmt_surf(self, value: float) -> str:
        # area, potential, charge
        return f"{value:15.9f}"


    def _fmt_wleb(self, value: float) -> str:
        # w_leb: 9 decimals, no extra trailing zeros
        return f"{value:15.9f}"


    def _fmt_switch(self, value: float) -> str:
        # Switch_F: high precision, matches "1.00000000000000" style
        return f"{value:17.14f}"


    def _fmt_gwidth(self, value: float) -> str:
        # G_width: long, high precision
        return f"{value:17.14f}"


    def _fmt_cart(self, value):
        # Cartesian/radii: 9 decimals, width ~14 like your input
        return f"{value:14.9f}"


    def _fmt_meta(self, value):
        # metadata floats: keep 9 decimals
        return f"{value:15.9f}"
