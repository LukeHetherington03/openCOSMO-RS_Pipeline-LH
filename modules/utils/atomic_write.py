"""
AtomicWriter: a class-based atomic file writing utility.

Guarantees:
- Writes go to a temporary file in the same directory.
- Temp file is atomically replaced with the final file on success.
- If an exception occurs, the temp file is removed and the original file is untouched.
"""

import os
import tempfile


class AtomicWriter:
    """
    Class-based atomic file writer.

    Usage:
        with AtomicWriter(path, mode="w") as f:
            f.write("data")
    """

    def __init__(self, path, mode="w", encoding=None):
        self.final_path = os.fspath(path)
        self.mode = mode
        self.encoding = encoding

        # Temp file path will be created in __enter__
        self._tmp_fd = None
        self._tmp_path = None

    def __enter__(self):
        directory = os.path.dirname(self.final_path)

        # Create a temporary file in the same directory
        self._tmp_fd, self._tmp_path = tempfile.mkstemp(dir=directory)

        # Open the file descriptor as a Python file object
        self._file = os.fdopen(self._tmp_fd, self.mode, encoding=self.encoding)
        return self._file

    def __exit__(self, exc_type, exc, tb):
        self._file.close()

        if exc_type is None:
            # Success → atomically replace final file
            os.replace(self._tmp_path, self.final_path)
        else:
            # Failure → remove temp file
            try:
                os.remove(self._tmp_path)
            except OSError:
                pass

        # Returning False propagates exceptions
        return False


# Convenience alias to match your previous API
def atomic_write(path, mode="w", encoding=None):
    return AtomicWriter(path, mode=mode, encoding=encoding)
