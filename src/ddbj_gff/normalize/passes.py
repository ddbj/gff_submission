from __future__ import annotations

from dataclasses import dataclass

from ..model import Directive
from .report import Change

_SPECIES_URL = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}"


@dataclass
class NormalizeContext:
    vocab: object             # validate.vocab.Vocab
    seq_lengths: dict | None
    config: object            # NormalizeConfig


def pass_directives(doc, ctx) -> list:
    cfg = ctx.config
    changes: list = []

    if doc.gff_version is None:
        doc.directives.insert(0, Directive("##gff-version 3", "gff-version", "3"))
        changes.append(Change("add-directive", "gff-version", "added ##gff-version 3"))

    if doc.insdc_gff_version is None:
        v = cfg.insdc_gff_version
        doc.directives.append(Directive(f"#!insdc-gff-version {v}", "insdc-gff-version", v))
        changes.append(Change("add-directive", "insdc-gff-version", f"added #!insdc-gff-version {v}"))

    if not isinstance(doc.species, int):
        if cfg.taxid is not None:
            url = _SPECIES_URL.format(taxid=cfg.taxid)
            doc.directives.append(Directive(f"##species {url}", "species", cfg.taxid))
            changes.append(Change("add-directive", "species", f"added ##species (taxid {cfg.taxid})"))
        else:
            changes.append(Change("no-taxid", "species",
                                  "##species not added: no taxid (set [normalize].taxid or --taxid)"))

    have = set(doc.sequence_regions)
    seqids: list = []
    for f in doc.features:
        for s in f.spans:
            if s.seqid not in have and s.seqid not in seqids:
                seqids.append(s.seqid)
    for seqid in seqids:
        length = ctx.seq_lengths.get(seqid) if ctx.seq_lengths else None
        approx = length is None
        if approx:
            length = max((s.end for f in doc.features for s in f.spans if s.seqid == seqid), default=1)
        doc.directives.append(
            Directive(f"##sequence-region {seqid} 1 {length}", "sequence-region", (seqid, 1, length)))
        if approx:
            changes.append(Change("approx-region", seqid,
                                  f"added ##sequence-region {seqid} 1 {length} "
                                  f"(length approximated from max feature end; provide --fasta for true length)"))
        else:
            changes.append(Change("add-directive", seqid, f"added ##sequence-region {seqid} 1 {length}"))

    if doc.transl_table_map is None:
        if any(f.type == "CDS" for f in doc.features):
            vals = {f.transl_table for f in doc.features if f.type == "CDS" and f.transl_table is not None}
            n = vals.pop() if len(vals) == 1 else cfg.transl_table
            doc.directives.append(Directive(f"#!transl_table primary:{n}", "transl_table", {"primary": n}))
            changes.append(Change("add-directive", "transl_table", f"added #!transl_table primary:{n}"))

    return changes
