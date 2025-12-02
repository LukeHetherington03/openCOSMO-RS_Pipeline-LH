"""
conformer_pruning.py

Module for pruning conformers based on different strategies.
Reads a unified _conformer_energies.csv file produced during conformer generation,
applies pruning rules, and writes out a pruned_conformers_lookup_<method>.csv file.

Available methods and parameters
--------------------------------
energy_window
    - energy_window : float (default=5.0)
      Keep conformers within X kcal/mol of the lowest energy.

rmsd to add
    - rmsd_threshold : float (default=0.5)
      Placeholder for RMSD pruning; intended to remove conformers too similar
      based on RMSD clustering.

rot_bond to add
    - max_rot_bonds : int (default=8)
      Keep conformers with <= max_rot_bonds rotatable bonds.

topN
    - top_n : int (default=10)
      Keep the N lowest energy conformers.

percentile
    - lower_pct : float (default=0)
    - upper_pct : float (default=100)
      Keep conformers whose energies fall between the lower and upper percentiles
      of the distribution (values between 0–100).
"""

import os
import pandas as pd


class ConformerPruner:
    def __init__(self, energies_csv, output_dir):
        """
        Parameters
        ----------
        energies_csv : str
            Path to the _conformer_energies.csv file.
        output_dir : str
            Directory to write pruned lookup files.
        """
        self.energies_csv = energies_csv
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.df = pd.read_csv(self.energies_csv)
        
    def prune(self, method="energy_window", output_file=None, **params):
        """
        Dispatch pruning to the chosen method and write a standardised lookup file.

        Parameters
        ----------
        method : str
            Pruning strategy. Supported values:
                - "energy_window" : keep conformers within X kcal/mol of lowest energy
                - "rmsd"          : prune conformers too similar by RMSD
                - "rot_bond"      : keep conformers with <= max_rot_bonds
                - "topN"          : keep N lowest energy conformers
                - "percentile"    : keep conformers within lower/upper percentile bounds
        output_file : str, optional
            Full path to the lookup file to be written. If None, a default name
            will be constructed inside self.output_dir.
        params : dict
            Parameters specific to the pruning method, e.g.:
                - energy_window : float (default=5.0)
                - rmsd_threshold : float (default=0.5)
                - max_rot_bonds : int (default=8)
                - top_n : int (default=10)
                - lower_pct : float (default=0)
                - upper_pct : float (default=100)

        Returns
        -------
        str
            Path to the lookup CSV file containing pruned conformers.
        """

        # --- Dispatch to the correct pruning method ---
        if method == "energy_window":
            df_pruned = self._prune_energy_window(**params)
        elif method == "rmsd":
            df_pruned = self._prune_rmsd(**params)
        elif method == "rot_bond":
            df_pruned = self._prune_rot_bond(**params)
        elif method == "topN":
            df_pruned = self._prune_topN(**params)
        elif method == "percentile":
            df_pruned = self._prune_percentile(**params)
        else:
            raise ValueError(f"Unknown pruning method: {method}")

        # --- Build unique lookup identifiers ---
        # Each conformer is identified by its InChIKey + conf_id
        df_pruned["lookup_id"] = df_pruned["inchi_key"] + "_conf" + df_pruned["conf_id"].astype(str)

        # --- Determine output filename ---
        if output_file is None:
            # Default naming convention: lookup_<method>.csv
            output_file = os.path.join(self.output_dir, f"lookup_{method}.csv")

        # --- Write lookup file ---
        df_pruned[["lookup_id"]].to_csv(output_file, index=False)

        print(f"Pruned conformers written to {output_file}")
        return output_file



    # ------------------ Methods ------------------

    def _prune_topN(self, top_n=10):
        """Keep N lowest energy conformers per molecule."""
        pruned_groups = []
        for inchi, group in self.df.groupby("inchi_key"):
            pruned = group.sort_values("energy").head(top_n)
            pruned_groups.append(pruned)
        return pd.concat(pruned_groups)

    def _prune_energy_window(self, energy_window=5.0):
        """Keep conformers within X kcal/mol of lowest energy per molecule."""
        pruned_groups = []
        for inchi, group in self.df.groupby("inchi_key"):
            min_energy = group["energy"].min()
            pruned = group[group["energy"] <= min_energy + energy_window]
            pruned_groups.append(pruned)
        return pd.concat(pruned_groups)

    def _prune_percentile(self, lower_pct=0, upper_pct=100):
        """Keep conformers within energy percentiles per molecule."""
        pruned_groups = []
        for inchi, group in self.df.groupby("inchi_key"):
            lower_bound = group["energy"].quantile(lower_pct / 100.0)
            upper_bound = group["energy"].quantile(upper_pct / 100.0)
            pruned = group[(group["energy"] >= lower_bound) & (group["energy"] <= upper_bound)]
            pruned_groups.append(pruned)
        return pd.concat(pruned_groups)
