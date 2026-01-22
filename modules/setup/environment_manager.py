#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
modules/setup/environment_manager.py

Manages software and resource paths for the pipeline.
Reads from config/paths.json and validates installations.
"""

import os
import sys
import json
import subprocess
from pathlib import Path
from glob import glob


class EnvironmentManager:
    """
    Manages software and resource paths for the pipeline.

    Responsibilities:
    1. Load paths from config/paths.json
    2. Validate that paths exist
    3. Validate executables and modules
    4. Validate chemistry JSON files in chemistry_dir
    5. Add Python paths to sys.path
    6. Report on environment status
    """

    def __init__(self, config_path=None):

        if config_path is None:
            project_root = self._find_project_root()
            config_path = os.path.join(project_root, 'config', 'paths.json')

        self.config_path = config_path
        self.config_source = None

        # Software paths
        self.gxtb_home = None
        self.gxtb_executable = None
        self.xtb_home = None
        self.xtb_executable = None
        self.orca_home = None
        self.orca_executable = None

        # openCOSMO-RS paths
        self.opencosmo_python_home = None
        self.opencosmo_python_src = None
        self.opencosmo_cpp_bindings = None
        self.opencosmo_binary = None
        self.opencosmo_shared_object = None
        self.opencosmo_python_module = None
        self.opencosmo_cpp_module = None

        # Constant files
        self.constant_files_root = None
        self.constant_files_metadata_dir = None
        self.constant_files_solvent_dir = None
        self.chemistry_dir = None

        # Validation results
        self.validation_results = {
            'software': {},
            'resources': {},
            'chemistry': {},
            'tests': {}
        }

        self._load_config()

    # ------------------------------------------------------------
    # Config loading
    # ------------------------------------------------------------
    def _find_project_root(self):
        current = os.path.dirname(os.path.abspath(__file__))
        while current != '/':
            if os.path.exists(os.path.join(current, 'CONSTANT_FILES')):
                return current
            current = os.path.dirname(current)
        return os.path.dirname(os.path.dirname(current))

    def _load_config(self):
        if not os.path.exists(self.config_path):
            raise FileNotFoundError(
                f"Configuration file not found: {self.config_path}"
            )

        with open(self.config_path) as f:
            cfg = json.load(f)

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
        self.opencosmo_python_home = opencosmo.get("python_home")
        self.opencosmo_python_src = opencosmo.get("python_src")
        self.opencosmo_cpp_bindings = opencosmo.get("cpp_bindings")
        self.opencosmo_binary = opencosmo.get("binary")
        self.opencosmo_shared_object = opencosmo.get("shared_object")
        self.opencosmo_python_module = opencosmo.get("python_module", "opencosmorspy")
        self.opencosmo_cpp_module = opencosmo.get("cpp_module", "openCOSMORS")

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
        if path is None:
            return {'status': 'missing', 'path': None, 'message': f'{name} not configured'}

        if os.path.exists(path):
            return {'status': 'ok', 'path': path, 'message': f'{name} found'}
        else:
            return {'status': 'error', 'path': path, 'message': f'{name} not found'}

    def _validate_file_exists(self, name, filepath):
        if filepath and os.path.isfile(filepath):
            return {'status': 'ok', 'path': filepath, 'message': f'{name} found'}
        else:
            return {'status': 'error', 'path': filepath, 'message': f'{name} not found'}

    # ------------------------------------------------------------
    # Chemistry validation
    # ------------------------------------------------------------
    def validate_chemistry(self):
        """
        Validate chemistry folder and all JSON files inside it.
        Expandable: any new JSON file is automatically validated.
        """
        results = {}

        # Validate chemistry directory
        results['chemistry_dir'] = self._validate_path_exists('chemistry_dir', self.chemistry_dir)

        if results['chemistry_dir']['status'] != 'ok':
            self.validation_results['chemistry'] = results
            return False

        # Discover all JSON files
        json_files = sorted(glob(os.path.join(self.chemistry_dir, "*.json")))

        if not json_files:
            results['json_files'] = {
                'status': 'warning',
                'message': 'No chemistry JSON files found'
            }
            self.validation_results['chemistry'] = results
            return False

        # Validate each JSON file loads cleanly
        for jf in json_files:
            name = os.path.basename(jf)
            try:
                with open(jf) as f:
                    json.load(f)
                results[name] = {
                    'status': 'ok',
                    'path': jf,
                    'message': 'JSON loaded successfully'
                }
            except Exception as e:
                results[name] = {
                    'status': 'error',
                    'path': jf,
                    'message': f'JSON failed to load: {e}'
                }

        self.validation_results['chemistry'] = results
        return all(r.get('status') == 'ok' for r in results.values() if isinstance(r, dict))

    # ------------------------------------------------------------
    # Software validation
    # ------------------------------------------------------------
    def validate_software(self):
        results = {}

        # g-xTB
        results['gxtb_home'] = self._validate_path_exists('gxtb.home', self.gxtb_home)
        results['gxtb_executable'] = self._validate_file_exists('gxtb.executable', self.gxtb_executable)

        # XTB
        results['xtb_home'] = self._validate_path_exists('xtb.home', self.xtb_home)
        results['xtb_executable'] = self._validate_file_exists('xtb.executable', self.xtb_executable)

        # ORCA
        results['orca_home'] = self._validate_path_exists('orca.home', self.orca_home)
        results['orca_executable'] = self._validate_file_exists('orca.executable', self.orca_executable)

        self.validation_results['software'] = results
        return all(r.get('status') == 'ok' for r in results.values())

    # ------------------------------------------------------------
    # Resource validation
    # ------------------------------------------------------------
    def validate_resources(self):
        results = {}

        # openCOSMO-RS
        results['opencosmo_python_home'] = self._validate_path_exists('opencosmo.python_home', self.opencosmo_python_home)
        results['opencosmo_python_src'] = self._validate_path_exists('opencosmo.python_src', self.opencosmo_python_src)
        results['opencosmo_cpp_bindings'] = self._validate_path_exists('opencosmo.cpp_bindings', self.opencosmo_cpp_bindings)
        results['opencosmo_binary'] = self._validate_file_exists('opencosmo.binary', self.opencosmo_binary)
        results['opencosmo_shared_object'] = self._validate_file_exists('opencosmo.shared_object', self.opencosmo_shared_object)

        # Constant files
        results['constant_files_root'] = self._validate_path_exists('constant_files.root', self.constant_files_root)
        results['constant_files_metadata_dir'] = self._validate_path_exists('constant_files.metadata_dir', self.constant_files_metadata_dir)
        results['constant_files_solvent_dir'] = self._validate_path_exists('constant_files.solvent_dir', self.constant_files_solvent_dir)

        self.validation_results['resources'] = results
        return all(r.get('status') == 'ok' for r in results.values())

    # ------------------------------------------------------------
    # Functional tests
    # ------------------------------------------------------------
    def test_software(self):
        results = {}

        if self.gxtb_executable:
            results['gxtb'] = self._test_gxtb()

        if self.xtb_executable:
            results['xtb'] = self._test_executable(self.xtb_executable, '--version', 'XTB')

        if self.orca_executable:
            results['orca'] = self._test_executable(self.orca_executable, None, 'ORCA')

        self.validation_results['tests'].setdefault('software', {}).update(results)
        return results

    def _test_executable(self, executable, version_flag, name):
        if not executable or not os.path.isfile(executable):
            return {'status': 'error', 'message': f'{name} executable not found'}

        try:
            if version_flag:
                result = subprocess.run(
                    [executable, version_flag],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
            else:
                if os.access(executable, os.X_OK):
                    return {'status': 'ok', 'message': f'{name} is executable'}
                else:
                    return {'status': 'error', 'message': f'{name} exists but is not executable'}

            if result.returncode == 0:
                return {'status': 'ok', 'message': f'{name} runs successfully'}
            else:
                return {'status': 'warning', 'message': f'{name} returned non-zero exit code'}

        except Exception as e:
            return {'status': 'error', 'message': f'{name} test failed: {e}'}

    def _test_gxtb(self):
        """
        Safely test g-xTB by running it in a temporary directory
        with a minimal coord file. This avoids the usual crash
        when gxtb is run without input.
        """
        import tempfile

        if not self.gxtb_executable or not os.path.isfile(self.gxtb_executable):
            return {
                'status': 'error',
                'message': 'g-xTB executable not found'
            }

        try:
            with tempfile.TemporaryDirectory() as tmp:
                coord_path = os.path.join(tmp, "coord")

                # Minimal valid coord file
                with open(coord_path, "w") as f:
                    f.write("1\n\nH 0.0 0.0 0.0\n")

                # Run a trivial calculation
                result = subprocess.run(
                    [self.gxtb_executable, "coord", "--gfn", "0"],
                    cwd=tmp,
                    capture_output=True,
                    text=True,
                    timeout=10
                )

                if result.returncode == 0:
                    return {
                        'status': 'ok',
                        'message': 'g-xTB runs successfully'
                    }
                else:
                    return {
                        'status': 'warning',
                        'message': f'g-xTB returned non-zero exit code ({result.returncode})'
                    }

        except subprocess.TimeoutExpired:
            return {
                'status': 'warning',
                'message': 'g-xTB test timed out'
            }

        except Exception as e:
            return {
                'status': 'error',
                'message': f'g-xTB test failed: {e}'
            }


    # ------------------------------------------------------------
    # Python import tests
    # ------------------------------------------------------------
    def test_resources(self):
        results = {}

        if self.opencosmo_python_src:
            results['opencosmo_python'] = self._test_python_import(
                self.opencosmo_python_src,
                self.opencosmo_python_module,
                'openCOSMO-RS Python'
            )

        if self.opencosmo_cpp_bindings:
            results['opencosmo_cpp'] = self._test_python_import(
                self.opencosmo_cpp_bindings,
                self.opencosmo_cpp_module,
                'openCOSMO-RS C++'
            )

        self.validation_results['tests'].setdefault('resources', {}).update(results)
        return results

    def _test_python_import(self, path, module_name, description):
        original_path = sys.path.copy()

        try:
            if path not in sys.path:
                sys.path.insert(0, path)

            __import__(module_name)

            return {'status': 'ok', 'message': f'{description} imports successfully'}

        except ImportError as e:
            return {'status': 'error', 'message': f'{description} import failed: {e}'}

        finally:
            sys.path = original_path

    # ------------------------------------------------------------
    # Chemistry loader
    # ------------------------------------------------------------
    def load_chemistry_file(self, filename):
        """
        Load a chemistry JSON file from chemistry_dir.
        Example: load_chemistry_file("cpcm_radii.json")
        """
        if not self.chemistry_dir:
            raise RuntimeError("chemistry_dir is not configured in paths.json")

        full_path = os.path.join(self.chemistry_dir, filename)

        if not os.path.exists(full_path):
            raise FileNotFoundError(f"Chemistry file not found: {full_path}")

        with open(full_path) as f:
            return json.load(f)

    # ------------------------------------------------------------
    # Comprehensive Testing
    # ------------------------------------------------------------
    def test_all(self, verbose=True):
        """
        Run all validations and tests in one go.
        Returns True if everything passes, False otherwise.
        """

        if verbose:
            print("=" * 70)
            print("Running Complete Environment Validation")
            print("=" * 70)

        if verbose:
            print("\n[1/5] Validating Software Paths...")
        software_valid = self.validate_software()

        if verbose:
            print("[2/5] Validating Resource Paths...")
        resources_valid = self.validate_resources()

        if verbose:
            print("[3/5] Validating Chemistry Folder + JSONs...")
        chemistry_valid = self.validate_chemistry()

        if verbose:
            print("[4/5] Testing Software Executables...")
        software_tests = self.test_software()

        if verbose:
            print("[5/5] Testing Python Imports...")
        resource_tests = self.test_resources()

        all_software_tests_ok = all(
            t.get('status') == 'ok'
            for t in software_tests.values()
        )
        all_resource_tests_ok = all(
            t.get('status') == 'ok'
            for t in resource_tests.values()
        )

        overall_pass = (
            software_valid and
            resources_valid and
            chemistry_valid and
            all_software_tests_ok and
            all_resource_tests_ok
        )

        return overall_pass

    # ------------------------------------------------------------
    # Pretty Table Output
    # ------------------------------------------------------------
    def _print_table(self):
        """
        Print all validation results in a clean ASCII table.
        """

        def section_header(title):
            print("\n" + "=" * 70)
            print(title)
            print("=" * 70)

        def print_table(rows):
            if not rows:
                print("(no entries)")
                return

            col1 = max(len(r[0]) for r in rows) + 2
            col2 = max(len(r[1]) for r in rows) + 2
            col3 = max(len(r[2]) for r in rows) + 2

            print(f"{'Item'.ljust(col1)}{'Status'.ljust(col2)}Message")
            print("-" * (col1 + col2 + col3))

            for name, status, msg in rows:
                symbol = "✓" if status == "ok" else ("⚠" if status == "warning" else "✗")
                print(f"{name.ljust(col1)}{symbol + ' ' + status.ljust(col2-2)}{msg}")

        # SOFTWARE
        section_header("SOFTWARE")
        rows = []
        for name, r in self.validation_results['software'].items():
            rows.append((name, r.get('status', 'unknown'), r.get('message', '')))
        print_table(rows)

        # RESOURCES
        section_header("RESOURCES")
        rows = []
        for name, r in self.validation_results['resources'].items():
            rows.append((name, r.get('status', 'unknown'), r.get('message', '')))
        print_table(rows)

        # CHEMISTRY
        section_header("CHEMISTRY")
        rows = []
        for name, r in self.validation_results['chemistry'].items():
            if isinstance(r, dict):
                rows.append((name, r.get('status', 'unknown'), r.get('message', '')))
        print_table(rows)

        # TESTS
        section_header("TESTS")
        rows = []
        for category, tests in self.validation_results['tests'].items():
            for name, r in tests.items():
                rows.append((f"{category}:{name}", r.get('status', 'unknown'), r.get('message', '')))
        print_table(rows)

    # ------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------
    def report(self, verbose=True):
        print("=" * 70)
        print("Environment Manager - Configuration Report")
        print("=" * 70)
        print(f"Config loaded from: {self.config_source}\n")

        print("SOFTWARE")
        print("-" * 70)
        self._report_section('software', verbose)

        print("\nRESOURCES")
        print("-" * 70)
        self._report_section('resources', verbose)

        print("\nCHEMISTRY")
        print("-" * 70)
        self._report_section('chemistry', verbose)

        if self.validation_results.get('tests'):
            print("\nTESTS")
            print("-" * 70)
            self._report_tests(verbose)

        print("\n" + "=" * 70)
        self._report_summary()
        print("=" * 70)

    def _report_section(self, section, verbose):
        results = self.validation_results.get(section, {})
        for name, result in results.items():
            status = result.get('status', 'unknown')
            if status == 'ok':
                symbol = '✓'
            elif status == 'warning':
                symbol = '⚠'
            elif status == 'error':
                symbol = '✗'
            else:
                symbol = '?'
            print(f"{symbol} {name:25s} {result.get('message', '')}")
            if verbose and result.get('path'):
                print(f"  └─ Path: {result['path']}")

    def _report_tests(self, verbose):
        tests = self.validation_results.get('tests', {})
        for category, results in tests.items():
            print(f"\n{category.upper()}:")
            for name, result in results.items():
                status = result.get('status', 'unknown')
                if status == 'ok':
                    symbol = '✓'
                elif status == 'warning':
                    symbol = '⚠'
                elif status == 'error':
                    symbol = '✗'
                else:
                    symbol = '?'
                print(f"  {symbol} {name:20s} {result.get('message', '')}")

    def _report_summary(self):
        software_ok = all(
            r.get('status') == 'ok'
            for r in self.validation_results.get('software', {}).values()
        )
        resources_ok = all(
            r.get('status') == 'ok'
            for r in self.validation_results.get('resources', {}).values()
        )
        chemistry_ok = all(
            r.get('status') == 'ok'
            for r in self.validation_results.get('chemistry', {}).values()
            if isinstance(r, dict)
        )

        print("SUMMARY:")
        print(f"  Software:  {'✓ All OK' if software_ok else '✗ Issues found'}")
        print(f"  Resources: {'✓ All OK' if resources_ok else '✗ Issues found'}")
        print(f"  Chemistry: {'✓ All OK' if chemistry_ok else '✗ Issues found'}")

        if software_ok and resources_ok and chemistry_ok:
            print("\n✓ Environment is ready!")
        else:
            print("\n⚠ Please fix the issues above before running the pipeline")


# ================================================================
# CLI ENTRY POINT (argparse)
# ================================================================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Environment validation tool for openCOSMO-RS Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--software",
        action="store_true",
        help="Validate software paths and executables only"
    )

    parser.add_argument(
        "--resources",
        action="store_true",
        help="Validate resource paths and Python/C++ modules only"
    )

    parser.add_argument(
        "--chemistry",
        action="store_true",
        help="Validate chemistry folder and JSON files only"
    )

    parser.add_argument(
        "--imports",
        action="store_true",
        help="Test Python imports for openCOSMO-RS modules only"
    )

    parser.add_argument(
        "--all",
        action="store_true",
        help="Run full environment validation (default)"
    )

    parser.add_argument(
        "--table",
        action="store_true",
        help="Output results in table format"
    )

    args = parser.parse_args()

    # Default behaviour: run all tests
    if not (args.software or args.resources or args.chemistry or args.imports or args.all):
        args.all = True

    env = EnvironmentManager()

    # Run selected tests
    if args.all or args.software:
        env.validate_software()
        env.test_software()

    if args.all or args.resources:
        env.validate_resources()
        env.test_resources()

    if args.all or args.chemistry:
        env.validate_chemistry()

    if args.all or args.imports:
        env.test_resources()

    # Output
    if args.table:
        env._print_table()
    else:
        env.report(verbose=True)

    # Exit code
    ok = (
        all(r.get('status') == 'ok' for r in env.validation_results['software'].values()) and
        all(r.get('status') == 'ok' for r in env.validation_results['resources'].values()) and
        all(r.get('status') == 'ok' for r in env.validation_results['chemistry'].values()
            if isinstance(r, dict))
    )

    sys.exit(0 if ok else 1)
