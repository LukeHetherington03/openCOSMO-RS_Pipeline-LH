"""
Resource Manager
----------------
Provides stable, predictable resource allocation for the pipeline.

This version:
  - Uses fixed defaults (e.g., 19 cores, 32 GB RAM)
  - Allows optional overrides via Request parameters
  - Does not benchmark or estimate per-molecule resources
  - Keeps the pipeline reproducible and simple
"""

import os
import psutil


class ResourceManager:
    # Default machine profile (your known hardware)
    DEFAULT_CPUS = 19
    DEFAULT_MEMORY_GB = 32

    @staticmethod
    def system_resources():
        """Return total and available system resources."""
        mem = psutil.virtual_memory()
        return {
            "cpus_total": os.cpu_count(),
            "cpus_available": len(psutil.Process().cpu_affinity()),
            "memory_total_gb": mem.total / 1e9,
            "memory_available_gb": mem.available / 1e9,
        }

    @staticmethod
    def get_resources(request_params):
        """
        Merge user-requested resources with defaults.
        If user does not specify resources, use fixed defaults.
        """
        user = request_params.get("resources", {})

        cpus = user.get("cpus", ResourceManager.DEFAULT_CPUS)
        mem = user.get("memory_gb", ResourceManager.DEFAULT_MEMORY_GB)

        return {
            "cpus": cpus,
            "memory_gb": mem
        }
