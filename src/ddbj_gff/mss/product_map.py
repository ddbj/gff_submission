from __future__ import annotations


def load_product_map(path: str) -> dict[str, str]:
    """Read a 2-column TSV (id<TAB>product) into a dict. Blank lines skipped."""
    result: dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            key, product = parts[0].strip(), parts[1].strip()
            if key and product:
                result[key] = product
    return result
