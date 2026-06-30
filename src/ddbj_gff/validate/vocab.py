from __future__ import annotations

import functools
from dataclasses import dataclass, field
from pathlib import Path

_DATA = Path(__file__).parent / "data"


@dataclass(frozen=True)
class Vocab:
    feature_types: frozenset[str]
    insdc_map: dict[str, str]
    dbxref_dbtags: frozenset[str]
    feature_qualifiers: dict[str, tuple[str, ...]] = field(default_factory=dict)


def _is_concrete(quals: tuple[str, ...]) -> bool:
    return all("<" not in q and ">" not in q and "*" not in q for q in quals)


def _read_feature_mapping() -> tuple[frozenset[str], dict[str, str], dict[str, tuple[str, ...]]]:
    terms: set[str] = set()
    mapping: dict[str, str] = {}
    qualifiers: dict[str, tuple[str, ...]] = {}
    with open(_DATA / "feature-mapping.tsv", encoding="utf-8") as fh:
        next(fh, None)  # header
        for line in fh:
            if not line.strip():
                continue
            cols = [c.strip() for c in line.rstrip("\n").split("\t")]
            if len(cols) < 2 or not cols[1]:
                continue
            so_term = cols[1]
            terms.add(so_term)
            insdc = cols[3] if len(cols) > 3 and cols[3] else ""
            quals = tuple(c for c in cols[4:6] if c)  # 列5,6 (Qualifier 1,2)
            if not insdc:
                continue
            if so_term not in mapping:
                mapping[so_term] = insdc
                qualifiers[so_term] = quals
            elif not _is_concrete(qualifiers.get(so_term, ())) and _is_concrete(quals):
                qualifiers[so_term] = quals  # 重複: 具体値の行を優先
    return frozenset(terms), mapping, qualifiers


def _read_dbxref() -> frozenset[str]:
    tags: set[str] = set()
    with open(_DATA / "dbxref.tsv", encoding="utf-8") as fh:
        for line in fh:
            t = line.strip()
            if t and not t.startswith("#"):
                tags.add(t)
    return frozenset(tags)


@functools.lru_cache(maxsize=1)
def load_vocab() -> Vocab:
    terms, mapping, quals = _read_feature_mapping()
    return Vocab(
        feature_types=terms,
        insdc_map=mapping,
        dbxref_dbtags=_read_dbxref(),
        feature_qualifiers=quals,
    )
