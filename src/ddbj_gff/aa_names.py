"""Amino-acid abbreviation <-> INSDC full-name mappings for special-case canonicalization."""
from __future__ import annotations

ABBREV_TO_NAME: dict[str, str] = {
    "Ala": "alanine", "Arg": "arginine", "Asn": "asparagine", "Asp": "aspartic acid",
    "Cys": "cysteine", "Gln": "glutamine", "Glu": "glutamic acid", "Gly": "glycine",
    "His": "histidine", "Ile": "isoleucine", "Leu": "leucine", "Lys": "lysine",
    "Met": "methionine", "Phe": "phenylalanine", "Pro": "proline", "Ser": "serine",
    "Thr": "threonine", "Trp": "tryptophan", "Tyr": "tyrosine", "Val": "valine",
    "Sec": "selenocysteine", "Pyl": "pyrrolysine",
}
NAME_TO_ABBREV: dict[str, str] = {v: k for k, v in ABBREV_TO_NAME.items()}
_STOP = {"Term", "TERM", "*"}


def full_name(abbrev: str) -> str:
    return ABBREV_TO_NAME.get(abbrev, abbrev)


def to_abbrev(name: str) -> str:
    return NAME_TO_ABBREV.get(name, name)


def is_stop(aa: str) -> bool:
    return aa in _STOP
