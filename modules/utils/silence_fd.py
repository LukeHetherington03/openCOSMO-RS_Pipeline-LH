# modules/utils/silence_fd.py

import os
import sys

class SilenceFD:
    """
    Temporarily redirect process-level stdout/stderr to /dev/null.
    This silences noisy C/C++ libraries (e.g. OpenBabel GA output).
    """

    def __enter__(self):
        self._stdout_fd = sys.stdout.fileno()
        self._stderr_fd = sys.stderr.fileno()

        # Save originals
        self._saved_stdout = os.dup(self._stdout_fd)
        self._saved_stderr = os.dup(self._stderr_fd)

        # Open /dev/null
        self._null = os.open(os.devnull, os.O_WRONLY)

        # Redirect both to /dev/null
        os.dup2(self._null, self._stdout_fd)
        os.dup2(self._null, self._stderr_fd)

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Restore originals
        os.dup2(self._saved_stdout, self._stdout_fd)
        os.dup2(self._saved_stderr, self._stderr_fd)

        os.close(self._null)
        os.close(self._saved_stdout)
        os.close(self._saved_stderr)
