from __future__ import annotations

import tomllib
from dataclasses import dataclass


@dataclass
class NormalizeConfig:
    taxid: int | None = None
    transl_table: int = 1
    insdc_gff_version: str = "1.0.0"
    coerce_transcript_to_mrna: bool = True


def load_normalize_config(path: str) -> NormalizeConfig:
    with open(path, "rb") as fh:
        data = tomllib.load(fh)
    n = data.get("normalize", {})
    return NormalizeConfig(
        taxid=n.get("taxid"),
        transl_table=n.get("transl_table", 1),
        insdc_gff_version=n.get("insdc_gff_version", "1.0.0"),
        coerce_transcript_to_mrna=n.get("coerce_transcript_to_mrna", True),
    )
