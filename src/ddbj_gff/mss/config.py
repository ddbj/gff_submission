from __future__ import annotations

import tomllib
from dataclasses import dataclass, field

from ..errors import Diagnostic, GffParseError, Severity

_KNOWN_SECTIONS = {"source", "locus_tag", "cds", "assembly_gap", "product"}


@dataclass
class MssConfig:
    source: dict[str, str]
    chromosome_pattern: str | None = None
    ff_def_chromosome: str | None = None
    ff_def_default: str | None = None
    locus_tag_prefix: str = ""
    locus_tag_width: int = 6
    locus_tag_start: int = 10
    locus_tag_step: int = 10
    transl_table: int = 1
    gap_min_length: int = 10
    gap_type: str = "within scaffold"
    gap_linkage_evidence: str = "align genus"
    gap_estimated_length: str = "known"
    product_default: str = "hypothetical protein"


def load_config(path: str) -> tuple[MssConfig, list[Diagnostic]]:
    with open(path, "rb") as fh:
        data = tomllib.load(fh)
    diags: list[Diagnostic] = []
    for key in data:
        if key not in _KNOWN_SECTIONS:
            diags.append(Diagnostic(Severity.WARNING, None, "unknown-config-key",
                                    f"unknown config section [{key}] ignored"))

    if "source" not in data or not isinstance(data["source"], dict):
        raise GffParseError(Diagnostic(Severity.ERROR, None, "config-missing",
                                       "config requires a [source] table"))
    source_tbl = data["source"]
    chromosome = source_tbl.get("chromosome", {})
    ff = source_tbl.get("ff_definition", {})
    # plain source qualifiers = scalar entries of [source] (exclude the nested tables)
    source_qual = {k: str(v) for k, v in source_tbl.items()
                   if k not in ("chromosome", "ff_definition")}

    lt = data.get("locus_tag", {})
    if "prefix" not in lt:
        raise GffParseError(Diagnostic(Severity.ERROR, None, "config-missing",
                                       "config requires [locus_tag].prefix"))
    cds = data.get("cds", {})
    gap = data.get("assembly_gap", {})
    product = data.get("product", {})

    cfg = MssConfig(
        source=source_qual,
        chromosome_pattern=chromosome.get("pattern"),
        ff_def_chromosome=ff.get("chromosome"),
        ff_def_default=ff.get("default"),
        locus_tag_prefix=lt["prefix"],
        locus_tag_width=lt.get("width", 6),
        locus_tag_start=lt.get("start", 10),
        locus_tag_step=lt.get("step", 10),
        transl_table=cds.get("transl_table", 1),
        gap_min_length=gap.get("min_length", 10),
        gap_type=gap.get("gap_type", "within scaffold"),
        gap_linkage_evidence=gap.get("linkage_evidence", "align genus"),
        gap_estimated_length=gap.get("estimated_length", "known"),
        product_default=product.get("default", "hypothetical protein"),
    )
    return cfg, diags


def load_common(path: str) -> list[str]:
    with open(path, encoding="ascii") as fh:
        rows = [line.rstrip("\r\n") for line in fh]
    while rows and rows[-1] == "":
        rows.pop()
    if not rows or not rows[0].startswith("COMMON"):
        raise GffParseError(Diagnostic(Severity.ERROR, None, "common-invalid",
                                       "common metadata must be non-empty and start with COMMON"))
    return rows
