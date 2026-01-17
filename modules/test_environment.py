#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_environment.py

Simplified test script using the test_all() method.
Run from project root: python3 test_environment.py
"""

import sys
import os

# Add modules to path
project_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(project_root, 'modules'))

# Import the environment manager
from modules.setup.environment_manager import EnvironmentManager

def main():
    print("\n" + "=" * 70)
    print("Testing openCOSMO-RS Pipeline Environment")
    print("=" * 70 + "\n")
    
    # Create environment manager
    try:
        env = EnvironmentManager()
        print(f"✓ Configuration loaded from: {env.config_source}\n")
    except FileNotFoundError as e:
        print(f"✗ Error: {e}\n")
        print("Please create config/paths.json from the template:")
        print("  cp config/paths.json.template config/paths.json")
        return 1
    
    # Run all tests
    print("Running comprehensive validation...\n")
    all_pass = env.test_all(verbose=True)
    
    # Show detailed report
    print("\n")
    env.report(verbose=True)
    
    # Return exit code
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())