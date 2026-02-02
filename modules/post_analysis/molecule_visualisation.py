#!/usr/bin/env python3
"""
Molecule + COSMO visualisation

Usage:
    python3 -m modules.post_analysis.molecule_visualisation <bundle.json>
"""

import json
import os
import sys
import numpy as np
import pyvista as pv

OUTPUT_DIR = "/home/lunet/cglh4/openCOSMO-RS_Pipeline-LH/post_analysis_results/visualisations"
os.makedirs(OUTPUT_DIR, exist_ok=True)


class MoleculeVisualizer:
    def __init__(self, bundle_path: str):
        if not os.path.exists(bundle_path):
            raise FileNotFoundError(f"Bundle JSON not found: {bundle_path}")

        with open(bundle_path) as f:
            self.bundle = json.load(f)

        paths = self.bundle.get("paths", {})
        self.xyz_path = paths.get("xyz")
        self.cpcm_path = paths.get("cpcm")

        if not self.xyz_path:
            raise ValueError("Bundle JSON missing 'xyz' path.")
        if not self.cpcm_path:
            raise ValueError("Bundle JSON missing 'cpcm' path.")

        self.inchi_key = self.bundle.get("meta", {}).get("inchi_key", "molecule")

        self.atom_coords = None
        self.atom_symbols = None
        self.cosmo_mesh = None

    # ---------------- XYZ loading ----------------
    def load_xyz(self):
        if not os.path.exists(self.xyz_path):
            raise FileNotFoundError(f"XYZ file not found: {self.xyz_path}")

        symbols, coords = [], []
        with open(self.xyz_path) as f:
            lines = f.readlines()

        for line in lines[2:]:
            parts = line.split()
            if len(parts) < 4:
                continue
            sym = parts[0]
            x, y, z = map(float, parts[1:4])
            symbols.append(sym)
            coords.append([x, y, z])

        self.atom_symbols = symbols
        self.atom_coords = np.array(coords)

    # ---------------- COSMO surface loading ----------------
    def load_cpcm_surface(self):
        if not os.path.exists(self.cpcm_path):
            raise FileNotFoundError(f"CPCM JSON not found: {self.cpcm_path}")

        with open(self.cpcm_path) as f:
            cpcm = json.load(f)

        pts, sigma = [], []
        for row in cpcm["surface_points"]:
            pts.append([row[0], row[1], row[2]])
            sigma.append(row[7])

        pts = np.array(pts)
        sigma = np.array(sigma)

        cloud = pv.PolyData(pts)
        cloud["sigma"] = sigma

        # Reconstruct a smooth surface from the shell-like point cloud
        surface = cloud.reconstruct_surface()
        self.cosmo_mesh = surface.interpolate(cloud, sharpness=2.0)

    # ---------------- Bond inference ----------------
    def infer_bonds(self, max_dist=1.8, min_dist=0.4):
        coords = self.atom_coords
        n = len(coords)
        lines = []

        for i in range(n):
            for j in range(i + 1, n):
                d = np.linalg.norm(coords[i] - coords[j])
                if min_dist < d < max_dist:
                    lines.extend([2, i, j])

        return np.array(lines) if lines else None

    # ---------------- Legend spheres ----------------
    def add_legend_spheres(self, plotter, element_colors, used_elements):
        coords = self.atom_coords
        mins = coords.min(axis=0)
        maxs = coords.max(axis=0)
        span = maxs - mins
        span[span == 0.0] = 1.0  # avoid zeros

        # Place legend to the right of the molecule
        base = np.array([maxs[0] + 0.6 * span[0], maxs[1], maxs[2]])
        dy = 0.25 * span[1]
        radius = 0.08 * span[0]

        for i, elem in enumerate(used_elements):
            color = element_colors.get(elem, "white")
            center = base - np.array([0.0, i * dy, 0.0])
            sphere = pv.Sphere(radius=radius, center=center)
            plotter.add_mesh(sphere, color=color)
            plotter.add_point_labels(
                [center + np.array([radius * 1.8, 0.0, 0.0])],
                [elem],
                font_size=12,
                text_color="white",
                point_size=0,
            )

    # ---------------- Main visual ----------------
    def full_visual(self, save_prefix=None):
        if self.atom_coords is None:
            self.load_xyz()
        if self.cosmo_mesh is None:
            self.load_cpcm_surface()

        if save_prefix is None:
            save_prefix = f"{self.inchi_key}_cosmo"

        png_path = os.path.join(OUTPUT_DIR, f"{save_prefix}.png")

        headless = (os.environ.get("DISPLAY", "") == "")
        if headless:
            print("[INFO] No DISPLAY found — enabling off-screen rendering.")
            pv.OFF_SCREEN = True

        plotter = pv.Plotter(off_screen=headless)

        # Element-based colours
        element_colors = {
            "H": "lightgray",
            "C": "black",
            "O": "red",
            "N": "blue",
            "S": "yellow",
            "F": "green",
            "Cl": "green",
            "Br": "brown",
            "I": "purple",
        }

        coords = self.atom_coords
        symbols = self.atom_symbols
        unique_elements = sorted(set(symbols))

        # COSMO surface
        plotter.add_mesh(
            self.cosmo_mesh,
            scalars="sigma",
            cmap="coolwarm",
            smooth_shading=True,
            opacity=0.6,
        )

        # Atoms grouped by element
        for elem in unique_elements:
            idx = [i for i, s in enumerate(symbols) if s == elem]
            if not idx:
                continue
            pts = coords[idx]
            color = element_colors.get(elem, "white")
            plotter.add_points(
                pts,
                color=color,
                render_points_as_spheres=True,
                point_size=30.0,
            )

        # Bonds
        bonds = self.infer_bonds()
        if bonds is not None:
            bond_mesh = pv.PolyData()
            bond_mesh.points = coords
            bond_mesh.lines = bonds
            plotter.add_mesh(bond_mesh, color="gray", line_width=3)

        # Legend with actual coloured spheres
        self.add_legend_spheres(plotter, element_colors, unique_elements)

        plotter.add_axes()
        plotter.show_grid()

        if headless:
            plotter.camera_position = "xy"
            plotter.reset_camera()
            plotter.render()

            plotter.show(auto_close=False)
            plotter.screenshot(png_path)
            print(f"[INFO] Saved COSMO PNG: {png_path}")
            return

        plotter.show()


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 -m modules.post_analysis.molecule_visualisation <bundle.json>")
        sys.exit(1)

    bundle_path = sys.argv[1]
    viz = MoleculeVisualizer(bundle_path)
    viz.full_visual()


if __name__ == "__main__":
    main()
