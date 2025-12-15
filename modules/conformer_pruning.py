import os
import pandas as pd

class ConformerPruner:
    def __init__(self, master_csv, output_dir, lookup_csv=None, energy_source=None):
        """
        Initialise ConformerPruner.

        Args:
            master_csv (str): Path to master CSV containing conformer data.
            output_dir (str): Directory to write pruned lookup CSVs.
            lookup_csv (str, optional): Existing lookup CSV to restrict pruning subset.
            energy_source (str, optional): Explicit energy column to use (e.g. 'energy_dft').
                                           If None, falls back through chain.
        """
        self.master_csv = master_csv
        self.output_dir = output_dir
        self.lookup_csv = lookup_csv
        os.makedirs(self.output_dir, exist_ok=True)

        df_master = pd.read_csv(self.master_csv)

        if lookup_csv and os.path.exists(lookup_csv):
            df_lookup = pd.read_csv(lookup_csv)
            df_master = df_master[df_master["lookup_id"].isin(df_lookup["lookup_id"])]

        # --- Energy selection ---
        energy_cols = ["energy_dft", "energy_gxtb", "energy_xtb", "energy_ff", "energy_gen"]

        chosen_col = None
        if energy_source:
            # User override
            if energy_source not in df_master.columns:
                raise ValueError(f"Requested energy_source '{energy_source}' not found in master CSV.")
            candidate = pd.to_numeric(df_master[energy_source], errors="coerce")
            if candidate.notna().any():
                df_master["energy_active"] = candidate
                chosen_col = energy_source
            else:
                raise ValueError(f"Requested energy_source '{energy_source}' has no usable values.")
        else:
            # Fallback chain
            for col in energy_cols:
                if col in df_master.columns:
                    candidate = pd.to_numeric(df_master[col], errors="coerce")
                    if candidate.notna().any():
                        df_master["energy_active"] = candidate
                        chosen_col = col
                        break

        if not chosen_col:
            raise ValueError("No usable energy column found in master CSV.")

        self.energy_source = chosen_col
        df_master = df_master.dropna(subset=["energy_active"])
        df_master["energy_rank"] = df_master.groupby("inchi_key")["energy_active"].rank(method="first")
        self.df = df_master


    def prune(self, method="topN", output_file=None, **params):
        """
        Apply pruning method and write lookup CSV.
        """
        # --- Dispatch ---
        if method == "topN":
            df_pruned = self._prune_topN(**params)
        elif method == "energy_window":
            df_pruned = self._prune_energy_window(**params)
        elif method == "percentile":
            df_pruned = self._prune_percentile(**params)
        elif method == "rmsd":
            df_pruned = self._prune_rmsd(**params)
        elif method == "rot_bond":
            df_pruned = self._prune_rot_bond(**params)
        else:
            raise ValueError(f"Unknown pruning method: {method}")

        # Build lookup IDs
        df_pruned["lookup_id"] = df_pruned["inchi_key"] + "_conf" + df_pruned["conf_id"].astype(str)
        
        # Default output filename
        if output_file is None:
            output_file = os.path.join(self.output_dir, f"lookup_{method}.csv")

        # --- Metadata section ---
        total_in = len(self.df)
        usable = self.df["energy_active"].notna().sum()
        total_out = len(df_pruned)

        metadata_lines = [
            f"# master_csv: {self.master_csv}",
            f"# previous_lookup: {self.lookup_csv if hasattr(self, 'lookup_csv') else 'None'}",
            f"# method: {method}",
            f"# params: {params}",
            f"# energy_source: {self.energy_source}"
            f"# total_input: {total_in}",
            f"# usable_energy: {usable}",
            f"# total_pruned: {total_out}",
            ""
        ]

        # --- Write file with metadata + table ---
        with open(output_file, "w") as f:
            for line in metadata_lines:
                f.write(line + "\n")
            df_pruned[["lookup_id", "energy_active", "energy_rank"]].to_csv(f, index=False)

        print(f"Pruned conformers written to {output_file}")

    
        return output_file, len(df_pruned), self.energy_source



    # ------------------ Methods ------------------

    def _prune_topN(self, top_n=10):
        """Keep N lowest energy conformers per molecule."""
        pruned_groups = []
        for inchi, group in self.df.groupby("inchi_key"):
            pruned = group.sort_values("energy_active").head(top_n)
            pruned_groups.append(pruned)
        return pd.concat(pruned_groups)

    def _prune_energy_window(self, energy_window=5.0):
        """Keep conformers within X kcal/mol of lowest energy per molecule."""
        pruned_groups = []
        for inchi, group in self.df.groupby("inchi_key"):
            min_energy = group["energy_active"].min()
            pruned = group[group["energy_active"] <= min_energy + energy_window]
            pruned_groups.append(pruned)
        return pd.concat(pruned_groups)

    def _prune_percentile(self, lower_pct=0, upper_pct=100):
        """Keep conformers within energy percentiles per molecule."""
        pruned_groups = []
        for inchi, group in self.df.groupby("inchi_key"):
            lower_bound = group["energy_active"].quantile(lower_pct / 100.0)
            upper_bound = group["energy_active"].quantile(upper_pct / 100.0)
            pruned = group[(group["energy_active"] >= lower_bound) & (group["energy_active"] <= upper_bound)]
            pruned_groups.append(pruned)
        return pd.concat(pruned_groups)

    def _prune_rmsd(self, rmsd_threshold=0.5):
        """Placeholder for RMSD pruning."""
        raise NotImplementedError("RMSD pruning not yet implemented")

    def _prune_rot_bond(self, max_rot_bonds=8):
        """Placeholder for rotatable bond pruning."""
        raise NotImplementedError("Rotatable bond pruning not yet implemented")
