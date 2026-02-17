#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CLI commands for environment validation.
Integrates with EnvironmentManager.
"""

from modules.setup.environment_manager import EnvironmentManager


class EnvCommands:

    @staticmethod
    def dispatch(args):
        if not args:
            return EnvCommands._help()

        sub = args[0]

        if sub == "check":
            return EnvCommands.check()

        if sub == "software":
            return EnvCommands.software()

        if sub == "resources":
            return EnvCommands.resources()

        if sub == "chemistry":
            return EnvCommands.chemistry()

        if sub == "pip":
            return EnvCommands.pip_packages()

        if sub == "table":
            return EnvCommands.table()

        print(f"Unknown env command: {sub}")
        return EnvCommands._help()

    # ------------------------------------------------------------
    # HELP
    # ------------------------------------------------------------
    @staticmethod
    def _help():
        print("""
Environment Commands:
    pl env check         Run full environment validation
    pl env software      Validate external executables
    pl env resources     Validate openCOSMO paths + constant files
    pl env chemistry     Validate chemistry JSON files
    pl env pip           Validate required Python packages
    pl env table         Pretty table output
""")

    # ------------------------------------------------------------
    # LOAD ENVIRONMENT MANAGER
    # ------------------------------------------------------------
    @staticmethod
    def _load_env():
        return EnvironmentManager()

    # ------------------------------------------------------------
    # COMMANDS
    # ------------------------------------------------------------
    @staticmethod
    def check():
        env = EnvCommands._load_env()
        results = env.test_all()
        EnvCommands._print_summary(results)
        return results

    @staticmethod
    def software():
        env = EnvCommands._load_env()
        env.validate_software()
        EnvCommands._print_section(env.validation_results["software"], "Software Validation")

    @staticmethod
    def resources():
        env = EnvCommands._load_env()
        env.validate_resources()
        EnvCommands._print_section(env.validation_results["resources"], "Resource Validation")

    @staticmethod
    def chemistry():
        env = EnvCommands._load_env()
        env.validate_chemistry()
        EnvCommands._print_section(env.validation_results["chemistry"], "Chemistry Validation")

    @staticmethod
    def pip_packages():
        env = EnvCommands._load_env()
        env.validate_pip_packages()
        EnvCommands._print_section(env.validation_results["pip"], "Pip Package Validation")

    @staticmethod
    def table():
        env = EnvCommands._load_env()
        env.test_all()
        EnvCommands._print_table(env.validation_results)

    # ------------------------------------------------------------
    # OUTPUT HELPERS
    # ------------------------------------------------------------
    @staticmethod
    def _print_section(section, title):
        print("=" * 70)
        print(title)
        print("=" * 70)

        for name, result in section.items():
            status = result.get("status", "unknown")
            msg = result.get("message", "")
            path = result.get("path", None)
            symbol = EnvCommands._symbol(status)

            print(f"{symbol} {name:25s} {msg}")
            if path:
                print(f"    Path: {path}")

        print()

    @staticmethod
    def _print_summary(results):
        print("=" * 70)
        print("Environment Validation Summary")
        print("=" * 70)

        for category, section in results.items():
            print(f"\n[{category.upper()}]")
            for name, result in section.items():
                status = result.get("status", "unknown")
                msg = result.get("message", "")
                path = result.get("path", None)
                symbol = EnvCommands._symbol(status)

                print(f"  {symbol} {name:25s} {msg}")
                if path:
                    print(f"      Path: {path}")

        print("\nDone.\n")

    @staticmethod
    def _print_table(results):
        print("=" * 70)
        print("Environment Validation (Table View)")
        print("=" * 70)

        for category, section in results.items():
            print(f"\n{category.upper()}")
            print("-" * 70)

            for name, result in section.items():
                status = result.get("status", "unknown")
                msg = result.get("message", "")
                path = result.get("path", None)
                symbol = EnvCommands._symbol(status)

                print(f"{name:25s} {symbol} {msg}")
                if path:
                    print(f"    Path: {path}")

        print()

    @staticmethod
    def _symbol(status):
        if status == "ok":
            return "✓"
        if status == "warning":
            return "⚠"
        if status == "error":
            return "✗"
        return "?"
