import subprocess

def get_git_version():
    """
    Returns a short, human-readable Git version string.
    Example: 'main@abc1234' or 'unknown' if Git is unavailable.
    """
    try:
        sha = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.STDOUT
        ).decode().strip()

        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            stderr=subprocess.STDOUT
        ).decode().strip()

        return f"{branch}@{sha}"
    except Exception:
        return "unknown"
