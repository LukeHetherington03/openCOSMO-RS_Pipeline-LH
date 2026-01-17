#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
modules/utils/external_paths.py

Manages paths to external dependencies (openCOSMO-RS, etc.)
Reads from environment variables or config file.
"""

import os
import sys
import json


class ExternalPaths:
    """
    Singleton class to manage external dependency paths.
    
    Priority order:
    1. Environment variables
    2. config/paths.json
    3. Default fallback paths
    """
    
    _instance = None
    _initialized = False
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(ExternalPaths, cls).__new__(cls)
        return cls._instance
    
    def __init__(self):
        if not self._initialized:
            self._load_paths()
            self._initialized = True
    
    def _load_paths(self):
        """Load paths from environment or config file."""
        
        # Try to find project root (where CONSTANT_FILES is)
        current_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = current_dir
        while project_root != '/':
            if os.path.exists(os.path.join(project_root, 'CONSTANT_FILES')):
                break
            project_root = os.path.dirname(project_root)
        
        # Try to load from config file
        config_file = os.path.join(project_root, 'config', 'paths.json')
        config_data = {}
        
        if os.path.exists(config_file):
            with open(config_file) as f:
                config_data = json.load(f)
        
        # Load paths with priority: env vars > config file > defaults
        self.opencosmo_python_path = (
            os.environ.get('OPENCOSMO_PYTHON_PATH') or
            config_data.get('opencosmo_python_path') or
            '/home/lunet/cglh4/resources/raw_modules/openCOSMO-RS_py-master/src'
        )
        
        self.opencosmo_cpp_path = (
            os.environ.get('OPENCOSMO_CPP_PATH') or
            config_data.get('opencosmo_cpp_path') or
            '/home/lunet/cglh4/resources/openCOSMO-RS_cpp/bindings'
        )
        
        self.orca_command = (
            os.environ.get('ORCA_COMMAND') or
            config_data.get('orca_command') or
            'orca'
        )
        
        self.constant_files_dir = (
            os.environ.get('CONSTANT_FILES_DIR') or
            config_data.get('constant_files_dir') or
            os.path.join(project_root, 'CONSTANT_FILES')
        )
    
    def add_to_path(self):
        """Add external paths to sys.path for imports."""
        
        if self.opencosmo_python_path not in sys.path:
            sys.path.insert(0, self.opencosmo_python_path)
        
        if self.opencosmo_cpp_path not in sys.path:
            sys.path.insert(0, self.opencosmo_cpp_path)
    
    def verify(self):
        """Verify that paths exist and are accessible."""
        
        issues = []
        
        # Check Python wrapper
        if not os.path.exists(self.opencosmo_python_path):
            issues.append(f"openCOSMO Python wrapper not found: {self.opencosmo_python_path}")
        
        # Check C++ binary
        binary_path = os.path.join(self.opencosmo_cpp_path, 'openCOSMORS.so')
        if not os.path.exists(binary_path):
            issues.append(f"openCOSMO C++ binary not found: {binary_path}")
        
        # Check CONSTANT_FILES
        if not os.path.exists(self.constant_files_dir):
            issues.append(f"CONSTANT_FILES directory not found: {self.constant_files_dir}")
        
        return issues
    
    def print_status(self):
        """Print current path configuration."""
        
        print("=" * 60)
        print("External Dependency Paths")
        print("=" * 60)
        print(f"openCOSMO Python: {self.opencosmo_python_path}")
        print(f"openCOSMO C++:    {self.opencosmo_cpp_path}")
        print(f"ORCA command:     {self.orca_command}")
        print(f"Constant files:   {self.constant_files_dir}")
        print("=" * 60)
        
        issues = self.verify()
        if issues:
            print("\n⚠ Issues found:")
            for issue in issues:
                print(f"  - {issue}")
        else:
            print("\n✓ All paths verified")


# Singleton instance
paths = ExternalPaths()


# Convenience function for stages to use
def get_paths():
    """Get the ExternalPaths singleton."""
    return paths