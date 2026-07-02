from __future__ import annotations

import gzip
from contextlib import contextmanager


@contextmanager
def open_text(path: str, encoding: str = "utf-8", errors: str = "strict"):
    if path.endswith(".gz"):
        fh = gzip.open(path, "rt", encoding=encoding, errors=errors)
    else:
        fh = open(path, "r", encoding=encoding, errors=errors)
    try:
        yield fh
    finally:
        fh.close()
