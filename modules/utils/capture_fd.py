# modules/utils/capture_fd.py

import os
import sys
import tempfile

class CaptureFD:
    """
    Capture all output written to the process-level stdout/stderr
    (including C/C++ libraries like OpenBabel).
    """

    def __enter__(self):
        # File descriptors for Python stdout/stderr
        self._stdout_fd = sys.stdout.fileno()
        self._stderr_fd = sys.stderr.fileno()

        # Save original FDs
        self._saved_stdout = os.dup(self._stdout_fd)
        self._saved_stderr = os.dup(self._stderr_fd)

        # Temporary file to capture output
        self._tmp = tempfile.TemporaryFile(mode="w+b")

        # Redirect stdout/stderr to temp file
        os.dup2(self._tmp.fileno(), self._stdout_fd)
        os.dup2(self._tmp.fileno(), self._stderr_fd)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Restore original stdout/stderr
        os.dup2(self._saved_stdout, self._stdout_fd)
        os.dup2(self._saved_stderr, self._stderr_fd)

        os.close(self._saved_stdout)
        os.close(self._saved_stderr)

        # Read captured output
        self._tmp.seek(0)
        self.captured = self._tmp.read().decode(errors="replace")
        self._tmp.close()
