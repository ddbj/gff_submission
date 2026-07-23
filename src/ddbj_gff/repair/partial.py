from __future__ import annotations


def is_partial(feature) -> bool:
    return feature.attributes.get("partial") == ["true"]


def partial_attrs(five_prime: bool, three_prime: bool, strand: str,
                  start: int, end: int) -> dict[str, list[str]]:
    """INSDC partial attributes for the given partial ends.

    5' maps to col4 (start) on +/./? strand and to col5 (end) on - strand;
    3' maps to the other. start_range applies to col4, end_range to col5.
    Value form: '.,<col4>' for start_range, '<col5>,.' for end_range.
    """
    if strand == "-":
        start_partial, end_partial = three_prime, five_prime
    else:
        start_partial, end_partial = five_prime, three_prime
    attrs: dict[str, list[str]] = {}
    if start_partial or end_partial:
        attrs["partial"] = ["true"]
    if start_partial:
        attrs["start_range"] = [f".,{start}"]
    if end_partial:
        attrs["end_range"] = [f"{end},."]
    return attrs
