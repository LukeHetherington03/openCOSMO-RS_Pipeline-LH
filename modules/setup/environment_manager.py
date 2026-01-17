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
    3. Validate that executables and modules work (unit-test style)
    4. Add Python paths to sys.path
    5. Report on environment status
    """
    
    def __init__(self, config_path=None):
        """
        Initialize the environment manager.
        
        Args:
            config_path: Path to paths.json (default: auto-detect from project root)
        """
        # Find project root and config
        if config_path is None:
            project_root = self._find_project_root()
            config_path = os.path.join(project_root, 'config', 'paths.json')
        
        self.config_path = config_path
        self.config_source = None  # Will track where config was loaded from
        
        # Software paths (compiled binaries)
        self.gxtb_home = None
        self.gxtb_executable = None
        self.xtb_home = None
        self.xtb_executable = None
        self.orca_home = None
        self.orca_executable = None
        
        # Resource paths (Python modules, C++ bindings)
        self.opencosmo_home = None
        self.opencosmo_python_src = None
        self.opencosmo_cpp_bindings = None
        self.opencosmo_binary = None

        # Optional: explicit Python module name (else default)
        self.opencosmo_python_module = None
        # Optional: explicit C++ module name (else default)
        self.opencosmo_cpp_module = None
        
        # Other paths
        self.constant_files_dir = None
        
        # Validation results
        self.validation_results = {
            'software': {},
            'resources': {},
            'tests': {}
        }
        
        # Load configuration
        self._load_config()
    
    def _find_project_root(self):
        """Find project root by looking for CONSTANT_FILES directory."""
        current = os.path.dirname(os.path.abspath(__file__))
        while current != '/':
            if os.path.exists(os.path.join(current, 'CONSTANT_FILES')):
                return current
            current = os.path.dirname(current)
        
        # Fallback to assuming we're in modules/setup/
        return os.path.dirname(os.path.dirname(current))
    
    def _load_config(self):
        """Load paths from config/paths.json."""
        if not os.path.exists(self.config_path):
            raise FileNotFoundError(
                f"Configuration file not found: {self.config_path}\n"
                f"Please create it from the template:\n"
                f"  cp config/paths.json.template config/paths.json"
            )
        
        with open(self.config_path) as f:
            config = json.load(f)
        
        # Software paths
        self.gxtb_home = config.get('gxtb_home')
        self.gxtb_executable = config.get('gxtb_executable')
        self.xtb_home = config.get('xtb_home')
        self.xtb_executable = config.get('xtb_executable')
        self.orca_home = config.get('orca_home')
        self.orca_executable = config.get('orca_executable')
        
        # Resource paths
        self.opencosmo_home = config.get('opencosmo_home')
        self.opencosmo_python_src = config.get('opencosmo_python_src')
        self.opencosmo_cpp_bindings = config.get('opencosmo_cpp_bindings')
        self.opencosmo_binary = config.get('opencosmo_binary')

        # Optional module names
        self.opencosmo_python_module = config.get('opencosmo_python_module', 'opencosmorspy')
        self.opencosmo_cpp_module = config.get('opencosmo_cpp_module', 'openCOSMORS')
        
        # Other paths
        self.constant_files_dir = config.get('constant_files_dir')
        
        self.config_source = self.config_path
    
    # ================================================================
    # Path Validation (checks if paths exist)
    # ================================================================
    
    def validate_software(self):
        """
        Validate all software installations.
        Checks if paths exist and executables are present.
        """
        results = {}
        
        # g-xTB validation
        results['gxtb'] = self._validate_path_exists('gxtb_home', self.gxtb_home)
        results['gxtb_executable'] = self._validate_file_exists('gxtb executable', self.gxtb_executable)
        
        # XTB validation
        results['xtb'] = self._validate_path_exists('xtb_home', self.xtb_home)
        results['xtb_executable'] = self._validate_file_exists('xtb executable', self.xtb_executable)
        
        # ORCA validation
        results['orca'] = self._validate_path_exists('orca_home', self.orca_home)
        results['orca_executable'] = self._validate_file_exists('orca executable', self.orca_executable)
        
        self.validation_results['software'] = results
        return all(r.get('status') == 'ok' for r in results.values())
    
    def validate_resources(self):
        """
        Validate all resource paths.
        Checks if Python modules and C++ bindings are accessible.
        """
        results = {}
        
        # openCOSMO-RS home
        results['opencosmo_home'] = self._validate_path_exists('opencosmo_home', self.opencosmo_home)
        
        # Python source
        results['opencosmo_python_src'] = self._validate_path_exists('opencosmo_python_src', self.opencosmo_python_src)
        
        # C++ bindings directory
        results['opencosmo_cpp_bindings'] = self._validate_path_exists('opencosmo_cpp_bindings', self.opencosmo_cpp_bindings)
        
        # C++ binary
        results['opencosmo_binary'] = self._validate_file_exists('openCOSMORS binary', self.opencosmo_binary)

        # C++ .so (auto-detect in bindings dir)
        cpp_so_result = self._detect_cpp_so()
        results['opencosmo_cpp_so'] = cpp_so_result
        
        # Constant files
        results['constant_files'] = self._validate_path_exists('constant_files_dir', self.constant_files_dir)
        
        self.validation_results['resources'] = results
        return all(r.get('status') == 'ok' for r in results.values())
    
    def _detect_cpp_so(self):
        """Detect the openCOSMORS .so file in the bindings directory."""
        if not self.opencosmo_cpp_bindings or not os.path.isdir(self.opencosmo_cpp_bindings):
            return {
                'status': 'error',
                'path': self.opencosmo_cpp_bindings,
                'message': 'C++ bindings directory not found for .so detection'
            }
        
        pattern = os.path.join(self.opencosmo_cpp_bindings, 'openCOSMORS*.so')
        matches = glob(pattern)
        
        if not matches:
            return {
                'status': 'error',
                'path': self.opencosmo_cpp_bindings,
                'message': 'No openCOSMORS .so file found in bindings directory'
            }
        
        # Take the first match
        so_path = matches[0]
        if os.path.isfile(so_path):
            return {
                'status': 'ok',
                'path': so_path,
                'message': 'openCOSMORS .so file found'
            }
        else:
            return {
                'status': 'error',
                'path': so_path,
                'message': 'Detected .so path is not a file'
            }
    
    def _validate_path_exists(self, name, path):
        """Check if a path exists."""
        if path is None:
            return {'status': 'missing', 'path': None, 'message': f'{name} not configured'}
        
        if os.path.exists(path):
            return {'status': 'ok', 'path': path, 'message': f'{name} found'}
        else:
            return {'status': 'error', 'path': path, 'message': f'{name} not found at path'}
    
    def _validate_file_exists(self, name, filepath):
        """Check if a specific file exists."""
        if filepath and os.path.isfile(filepath):
            return {'status': 'ok', 'path': filepath, 'message': f'{name} found'}
        else:
            return {'status': 'error', 'path': filepath, 'message': f'{name} not found'}
    
    # ================================================================
    # Functional Testing (unit-test style validation)
    # ================================================================
    
    def test_software(self):
        """
        Run functional tests on software (unit-test style).
        Actually executes commands to verify they work.
        """
        results = {}
        
        # Test g-xTB
        if self.gxtb_executable:
            results['gxtb'] = self._test_executable(self.gxtb_executable, '--version', 'g-xTB')
        
        # Test XTB
        if self.xtb_executable:
            results['xtb'] = self._test_executable(self.xtb_executable, '--version', 'XTB')
        
        # Test ORCA
        if self.orca_executable:
            results['orca'] = self._test_executable(self.orca_executable, None, 'ORCA')
        
        self.validation_results['tests']['software'] = results
        return results
    
    def test_resources(self):
        """
        Run functional tests on resources (unit-test style).
        Tests if Python modules and C++ bindings can be imported.
        """
        results = {}
        
        # Test openCOSMO-RS Python import
        if self.opencosmo_python_src:
            results['opencosmo_python'] = self._test_python_import(
                self.opencosmo_python_src,
                self.opencosmo_python_module,
                'openCOSMO-RS Python'
            )
        
        # Test openCOSMO-RS C++ import (Python extension)
        if self.opencosmo_cpp_bindings:
            results['opencosmo_cpp'] = self._test_python_import(
                self.opencosmo_cpp_bindings,
                self.opencosmo_cpp_module,
                'openCOSMO-RS C++'
            )
        
        self.validation_results['tests']['resources'] = results
        return results
    
    def _test_executable(self, executable, version_flag, name):
        """Test if an executable runs."""
        if not executable or not os.path.isfile(executable):
            return {
                'status': 'error',
                'message': f'{name} executable not found',
                'executable': executable
            }
        
        try:
            # Try to run with version flag
            if version_flag:
                result = subprocess.run(
                    [executable, version_flag],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
            else:
                # For ORCA, just check if it exists and is executable
                if os.access(executable, os.X_OK):
                    return {
                        'status': 'ok',
                        'message': f'{name} is executable',
                        'executable': executable
                    }
                else:
                    return {
                        'status': 'error',
                        'message': f'{name} exists but is not executable',
                        'executable': executable
                    }
            
            # Check if command succeeded
            if result.returncode == 0 or version_flag:
                return {
                    'status': 'ok',
                    'message': f'{name} runs successfully',
                    'executable': executable,
                    'output': result.stdout[:200] if result.stdout else result.stderr[:200]
                }
            else:
                return {
                    'status': 'warning',
                    'message': f'{name} ran but returned non-zero exit code',
                    'executable': executable,
                    'returncode': result.returncode
                }
        
        except subprocess.TimeoutExpired:
            return {
                'status': 'warning',
                'message': f'{name} command timed out',
                'executable': executable
            }
        except Exception as e:
            return {
                'status': 'error',
                'message': f'{name} test failed: {str(e)}',
                'executable': executable
            }
    
    def _test_python_import(self, path, module_name, description):
        """Test if a Python module can be imported."""
        original_path = sys.path.copy()
        
        try:
            if path not in sys.path:
                sys.path.insert(0, path)
            
            __import__(module_name)
            
            return {
                'status': 'ok',
                'message': f'{description} imports successfully',
                'module': module_name,
                'path': path
            }
        
        except ImportError as e:
            return {
                'status': 'error',
                'message': f'{description} import failed: {str(e)}',
                'module': module_name,
                'path': path
            }
        
        finally:
            sys.path = original_path
    
    # ================================================================
    # Comprehensive Testing
    # ================================================================
    
    def test_all(self, verbose=True):
        """
        Run all validations and tests in one go.
        Returns True if everything passes, False otherwise.
        """
        if verbose:
            print("=" * 70)
            print("Running Complete Environment Validation")
            print("=" * 70)
        
        # Step 1: Validate paths
        if verbose:
            print("\n[1/4] Validating Software Paths...")
        software_valid = self.validate_software()
        
        if verbose:
            print("[2/4] Validating Resource Paths...")
        resources_valid = self.validate_resources()
        
        # Step 2: Run functional tests
        if verbose:
            print("[3/4] Testing Software Executables...")
        software_tests = self.test_software()
        
        if verbose:
            print("[4/4] Testing Python Imports...")
        resource_tests = self.test_resources()
        
        # Step 3: Add to Python path
        paths_added = self.add_to_python_path()
        
        # Determine overall status
        all_software_tests_ok = all(
            t.get('status') == 'ok' 
            for t in software_tests.values()
        )
        all_resource_tests_ok = all(
            t.get('status') == 'ok' 
            for t in resource_tests.values()
        )
        
        overall_pass = (software_valid and resources_valid and 
                        all_software_tests_ok and all_resource_tests_ok)
        
        if verbose:
            print("\n" + "=" * 70)
            if overall_pass:
                print("✓ ALL TESTS PASSED - Environment is ready!")
            else:
                print("✗ SOME TESTS FAILED - See details above")
            print("=" * 70)
        
        return overall_pass
    
    # ================================================================
    # Environment Setup
    # ================================================================
    
    def add_to_python_path(self):
        """Add necessary paths to sys.path for Python imports."""
        paths_added = []
        
        # Add openCOSMO-RS Python source
        if self.opencosmo_python_src and self.opencosmo_python_src not in sys.path:
            sys.path.insert(0, self.opencosmo_python_src)
            paths_added.append(('opencosmo_python_src', self.opencosmo_python_src))
        
        # Add openCOSMO-RS C++ bindings
        if self.opencosmo_cpp_bindings and self.opencosmo_cpp_bindings not in sys.path:
            sys.path.insert(0, self.opencosmo_cpp_bindings)
            paths_added.append(('opencosmo_cpp_bindings', self.opencosmo_cpp_bindings))
        
        return paths_added
    
    def get_executable_path(self, software):
        """
        Get the full path to an executable.
        
        Args:
            software: 'gxtb', 'xtb', or 'orca'
        """
        if software == 'gxtb':
            return self.gxtb_executable
        elif software == 'xtb':
            return self.xtb_executable
        elif software == 'orca':
            return self.orca_executable
        else:
            return None
    
    # ================================================================
    # Reporting
    # ================================================================
    
    def report(self, verbose=True):
        """
        Print a detailed report of the environment status.
        """
        print("=" * 70)
        print("Environment Manager - Configuration Report")
        print("=" * 70)
        print(f"Config loaded from: {self.config_source}")
        print()
        
        # Software section
        print("SOFTWARE (Compiled Binaries)")
        print("-" * 70)
        self._report_section('software', verbose)
        
        # Resources section
        print("\nRESOURCES (Python Modules & Libraries)")
        print("-" * 70)
        self._report_section('resources', verbose)
        
        # Test results if available
        if self.validation_results.get('tests'):
            print("\nFUNCTIONAL TESTS")
            print("-" * 70)
            self._report_tests(verbose)
        
        # Summary
        print("\n" + "=" * 70)
        self._report_summary()
        print("=" * 70)
    
    def _report_section(self, section, verbose):
        """Report on a validation section."""
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
            
            print(f"{symbol} {name:25s} {result.get('message', 'No message')}")
            
            if verbose and result.get('path'):
                print(f"  └─ Path: {result['path']}")
    
    def _report_tests(self, verbose):
        """Report on functional tests."""
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
                
                print(f"  {symbol} {name:20s} {result.get('message', 'No message')}")
                
                if verbose and result.get('output'):
                    print(f"    Output: {result['output'][:100]}...")
    
    def _report_summary(self):
        """Print summary of validation status."""
        software_ok = all(
            r.get('status') == 'ok' 
            for r in self.validation_results.get('software', {}).values()
        )
        resources_paths_ok = all(
            r.get('status') == 'ok' 
            for r in self.validation_results.get('resources', {}).values()
        )
        resources_tests_ok = all(
            r.get('status') == 'ok'
            for r in self.validation_results.get('tests', {}).get('resources', {}).values()
        )
        
        resources_ok = resources_paths_ok and resources_tests_ok
        
        print("SUMMARY:")
        print(f"  Software:  {'✓ All OK' if software_ok else '✗ Issues found'}")
        print(f"  Resources: {'✓ All OK' if resources_ok else '✗ Issues found'}")
        
        if software_ok and resources_ok:
            print("\n✓ Environment is ready!")
        else:
            print("\n⚠ Please fix the issues above before running the pipeline")
