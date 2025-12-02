"""
resource_management.py

Module for resource estimation and benchmarking.
Handles RAM/CPU heuristics, timing, and performance logging.
"""

import time
from rdkit import Chem

class ResourceManager:
    @staticmethod
    def estimate_resources(mol):
        """Estimate RAM and cores needed based on molecule size."""
        heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() != "H")
        atoms = Chem.AddHs(mol).GetNumAtoms()

        ram = 2000  # MB
        if heavy_atoms > 15: ram = 4000
        if heavy_atoms > 25: ram = 8000

        cores = 1
        if atoms <= 8: cores = 2
        elif atoms <= 16: cores = 4
        elif atoms >= 32: cores = 8

        return ram, cores


class Benchmark:
    def __init__(self, label="task"):
        self.label = label
        self.start_time = None
        self.end_time = None

    def __enter__(self):
        self.start_time = time.monotonic()
        print(f"[Benchmark] Starting {self.label}...")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = time.monotonic()
        elapsed = self.end_time - self.start_time
        print(f"[Benchmark] Finished {self.label} in {elapsed:.2f} seconds")

    def duration(self):
        if self.start_time and self.end_time:
            return self.end_time - self.start_time
        return None
