#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Environment Manager for openCOSMO-RS Pipeline
---------------------------------------------

Validates:
- External executables (crest, gxtb, xtb, ORCA)
- openCOSMO-RS Python + C++ bindings
- Constant files
- Chemistry JSONs
- Required Python packages (pip installs)

Loads everything from config/paths.json.
"""

import os
import sys
import json
import subprocess
from glob import glob
import importlib


# ================================================================
# Required Python packages for the pipeline
# ================================================================
REQUIRED_PIP_PACKAGES = [
    "numpy",
    "scipy",
    "pandas",
    "rdkit",
    "openbabel",
]


class EnvironmentManager:

    def __init__(self, config_path=None):

        if config_path is None:
            project_root = self._find_project_root()
            config_path = os.path.join(project_root, "config", "paths.json")

        self.config_path = config_path
        self.config_source = None

        # Software executables
        self.crest_home = None
        self.crest_executable = None
        self.gxtb_home = None
        self.gxtb_executable = None
        self.xtb_home = None
        self.xtb_executable = None
        self.orca_home = None
        self.orca_executable = None

        # openCOSMO-RS
        self.opencosmo_python_src = None
        self.opencosmo_cpp_bindings = None
        self.opencosmo_python_driver = None

        # Constant files
        self.constant_files_root = None
        self.constant_files_metadata_dir = None
        self.constant_files_solvent_dir = None
        self.chemistry_dir = None

        # Validation results
        self.validation_results = {
            "software": {},
            "resources": {},
            "chemistry": {},
            "pip": {},
            "tests": {},
        }

        self._load_config()

    # ------------------------------------------------------------
    # Config loading
    # ------------------------------------------------------------
    def _find_project_root(self):
        current = os.path.dirname(os.path.abspath(__file__))
        while current != "/":
            if os.path.exists(os.path.join(current, "CONSTANT_FILES")):
                return current
            current = os.path.dirname(current)
        return current

    def _load_config(self):
        if not os.path.exists(self.config_path):
            raise FileNotFoundError(f"paths.json not found: {self.config_path}")

        with open(self.config_path) as f:
            cfg = json.load(f)

        # crest
        crest = cfg.get("crest", {})
        self.crest_home = crest.get("home")
        self.crest_executable = crest.get("executable")

        # g-xTB
        gxtb = cfg.get("gxtb", {})
        self.gxtb_home = gxtb.get("home")
        self.gxtb_executable = gxtb.get("executable")

        # XTB
        xtb = cfg.get("xtb", {})
        self.xtb_home = xtb.get("home")
        self.xtb_executable = xtb.get("executable")

        # ORCA
        orca = cfg.get("orca", {})
        self.orca_home = orca.get("home")
        self.orca_executable = orca.get("executable")

        # openCOSMO-RS
        opencosmo = cfg.get("opencosmo", {})
        self.opencosmo_python_src = opencosmo.get("python_src")
        self.opencosmo_cpp_bindings = opencosmo.get("cpp_bindings")
        self.opencosmo_python_driver = opencosmo.get("python_driver")

        # Constant files
        const = cfg.get("constant_files", {})
        self.constant_files_root = const.get("root")
        self.constant_files_metadata_dir = const.get("metadata_dir")
        self.constant_files_solvent_dir = const.get("solvent_dir")
        self.chemistry_dir = const.get("chemistry_dir")

        self.config_source = self.config_path

    # ------------------------------------------------------------
    # Validation helpers
    # ------------------------------------------------------------
    def _validate_path_exists(self, name, path):
        if not path:
            return {"status": "missing", "path": None, "message": f"{name} not configured"}

        if os.path.exists(path):
            return {"status": "ok", "path": path, "message": f"{name} found"}
        else:
            return {"status": "error", "path": path, "message": f"{name} not found"}

    def _validate_file_exists(self, name, path):
        if path and os.path.isfile(path):
            return {"status": "ok", "path": path, "message": f"{name} found"}
        else:
            return {"status": "error", "path": path, "message": f"{name} not found"}

    def _validate_executable(self, name, path, version_flag="--version"):
        if not path:
            return {"status": "missing", "path": None, "message": f"{name} not configured"}

        if not os.path.isfile(path):
            return {"status": "error", "path": path, "message": f"{name} is not a file"}

        if not os.access(path, os.X_OK):
            return {"status": "error", "path": path, "message": f"{name} is not executable"}

        try:
            if version_flag:
                result = subprocess.run(
                    [path, version_flag],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
            else:
                result = subprocess.run(
                    [path],
                    capture_output=True,
                    text=True,
                    timeout=5
                )

            # ORCA returns exit code 1 normally — treat as OK
            if result.stdout or result.stderr:
                return {
                    "status": "ok",
                    "path": path,
                    "message": f"{name} executable detected"
                }

            return {
                "status": "warning",
                "path": path,
                "message": f"{name} ran but produced no output"
            }

        except Exception as e:
            return {"status": "error", "path": path, "message": f"{name} failed to run: {e}"}

    # ------------------------------------------------------------
    # g-xTB minimal test
    # ------------------------------------------------------------
    def _test_gxtb(self):
        path = self.gxtb_executable

        if not path or not os.path.isfile(path):
            return {
                "status": "error",
                "path": path,
                "message": "g-xTB executable not found"
            }

        import tempfile

        try:
            with tempfile.TemporaryDirectory() as tmp:
                coord_path = os.path.join(tmp, "coord")

                with open(coord_path, "w") as f:
                    f.write("1\n\nH 0.0 0.0 0.0\n")

                result = subprocess.run(
                    [path, "coord", "--gfn", "0"],
                    cwd=tmp,
                    capture_output=True,
                    text=True,
                    timeout=10
                )

                if result.returncode == 0:
                    return {
                        "status": "ok",
                        "path": path,
                        "message": "g-xTB runs successfully"
                    }
                else:
                    return {
                        "status": "warning",
                        "path": path,
                        "message": f"g-xTB returned {result.returncode}"
                    }

        except Exception as e:
            return {
                "status": "error",
                "path": path,
                "message": f"g-xTB test failed: {e}"
            }


    # ------------------------------------------------------------
    # Software validation
    # ------------------------------------------------------------
    def validate_software(self):
        results = {}

        results["crest"] = self._validate_executable("crest", self.crest_executable)
        results["gxtb"] = self._test_gxtb()
        results["xtb"] = self._validate_executable("xtb", self.xtb_executable)
        results["orca"] = self._validate_executable("orca", self.orca_executable, version_flag=None)

        self.validation_results["software"] = results
        return all(r["status"] == "ok" for r in results.values())

    # ------------------------------------------------------------
    # Resource validation
    # ------------------------------------------------------------
    def validate_resources(self):
        results = {}

        results["opencosmo_python_src"] = self._validate_path_exists(
            "opencosmo.python_src", self.opencosmo_python_src
        )
        results["opencosmo_cpp_bindings"] = self._validate_path_exists(
            "opencosmo.cpp_bindings", self.opencosmo_cpp_bindings
        )
        results["opencosmo_python_driver"] = self._validate_file_exists(
            "opencosmo.python_driver", self.opencosmo_python_driver
        )

        results["constant_files_root"] = self._validate_path_exists(
            "constant_files.root", self.constant_files_root
        )
        results["constant_files_metadata_dir"] = self._validate_path_exists(
            "constant_files.metadata_dir", self.constant_files_metadata_dir
        )
        results["constant_files_solvent_dir"] = self._validate_path_exists(
            "constant_files.solvent_dir", self.constant_files_solvent_dir
        )

        self.validation_results["resources"] = results
        return all(r["status"] == "ok" for r in results.values())

    # ------------------------------------------------------------
    # Chemistry validation
    # ------------------------------------------------------------
    def validate_chemistry(self):
        results = {}

        results["chemistry_dir"] = self._validate_path_exists("chemistry_dir", self.chemistry_dir)
        if results["chemistry_dir"]["status"] != "ok":
            self.validation_results["chemistry"] = results
            return False

        json_files = sorted(glob(os.path.join(self.chemistry_dir, "*.json")))
        if not json_files:
            results["json_files"] = {"status": "warning", "message": "No chemistry JSON files found"}
            self.validation_results["chemistry"] = results
            return False

        for jf in json_files:
            name = os.path.basename(jf)
            try:
                with open(jf) as f:
                    json.load(f)
                results[name] = {"status": "ok", "path": jf, "message": "JSON loaded successfully"}
            except Exception as e:
                results[name] = {"status": "error", "path": jf, "message": f"JSON failed: {e}"}

        self.validation_results["chemistry"] = results
        return all(r["status"] == "ok" for r in results.values() if isinstance(r, dict))

    # ------------------------------------------------------------
    # Pip package validation (with version extraction)
    # ------------------------------------------------------------
    def validate_pip_packages(self):
        results = {}

        for pkg in REQUIRED_PIP_PACKAGES:
            try:
                module = importlib.import_module(pkg)
                version = getattr(module, "__version__", "unknown")
                results[pkg] = {"status": "ok", "message": f"version {version}"}
            except ImportError:
                results[pkg] = {"status": "error", "message": "missing package"}

        self.validation_results["pip"] = results
        return all(r["status"] == "ok" for r in results.values())

    # ------------------------------------------------------------
    # Python import tests for openCOSMO
    # ------------------------------------------------------------
    def test_resources(self):
        results = {}

        if self.opencosmo_python_src:
            results["opencosmo_python"] = self._test_python_import(
                self.opencosmo_python_src, "opencosmorspy", "openCOSMO-RS Python"
            )

        if self.opencosmo_cpp_bindings:
            results["opencosmo_cpp"] = self._test_python_import(
                self.opencosmo_cpp_bindings, "openCOSMORS", "openCOSMO-RS C++"
            )

        self.validation_results["tests"]["resources"] = results
        return results

    def _test_python_import(self, path, module_name, description):
        original = sys.path.copy()
        try:
            sys.path.insert(0, path)
            __import__(module_name)
            return {"status": "ok", "message": f"{description} imports successfully"}
        except Exception as e:
            return {"status": "error", "message": f"{description} import failed: {e}"}
        finally:
            sys.path = original

    # ------------------------------------------------------------
    # Full validation
    # ------------------------------------------------------------
    def test_all(self):
        self.validate_software()
        self.validate_resources()
        self.validate_chemistry()
        self.validate_pip_packages()
        self.test_resources()
        return self.validation_results
