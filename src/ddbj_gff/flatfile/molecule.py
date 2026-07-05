from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass
class MoleculeInfo:
    taxid: int | None
    organism: str | None
    division: str | None
    topology: str          # "linear" | "circular"
    compartment: str       # "nuclear" | "organelle"
    hierarchy: str         # "three_level" | "two_level"
    transl_table: int


def _source_feature(rec):
    for f in rec.features:
        if f.type == "source":
            return f
    return None


def detect_molecule(rec) -> MoleculeInfo:
    src = _source_feature(rec)
    quals = src.qualifiers if src is not None else {}
    taxid = None
    for xref in quals.get("db_xref", []):
        m = re.match(r"taxon:(\d+)", xref)
        if m:
            taxid = int(m.group(1))
    organism = (quals.get("organism") or [rec.annotations.get("organism")])[0] \
        if (quals.get("organism") or rec.annotations.get("organism")) else None
    division = rec.annotations.get("data_file_division")
    topology = rec.annotations.get("topology") or "linear"
    compartment = "organelle" if quals.get("organelle") else "nuclear"
    hierarchy = "two_level" if compartment == "organelle" else "three_level"
    # transl_table: primary from a CDS qualifier, default 1
    tt = 1
    for f in rec.features:
        if f.type == "CDS" and f.qualifiers.get("transl_table"):
            tt = int(f.qualifiers["transl_table"][0])
            break
    return MoleculeInfo(taxid=taxid, organism=organism, division=division,
                        topology=topology, compartment=compartment,
                        hierarchy=hierarchy, transl_table=tt)
