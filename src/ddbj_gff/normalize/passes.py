from __future__ import annotations

import re
from dataclasses import dataclass

from ..model import Directive, Feature, Span
from .. import aa_names
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


def _is_placeholder(qual: str) -> bool:
    return "<" in qual or ">" in qual or "*" in qual


def _qualifier_to_attr(qual: str) -> tuple[str, str | None]:
    body = qual.lstrip("/")
    if "=" in body:
        key, val = body.split("=", 1)
        return key.strip(), val.strip().strip('"')
    return body.strip(), None  # valueless flag (e.g. /pseudo)


_COLLAPSE_TARGETS = {"gene", "mRNA", "CDS", "exon", "intron"}


def pass_so_terms(doc, ctx) -> list:
    vocab = ctx.vocab
    changes: list = []
    for f in doc.features:
        target = vocab.insdc_map.get(f.type)
        if target is None:
            changes.append(Change("unmapped-type", f.id or "?",
                                  f"feature type {f.type!r} is not a known SO term; left unchanged"))
            continue
        if target == f.type or target not in _COLLAPSE_TARGETS:
            # already a core type, or maps to a non-core INSDC feature that Phase 2
            # handles during generation (ncRNA/precursor_RNA/5'UTR/...) -> leave as-is
            continue
        old = f.type
        f.type = target
        changes.append(Change("rename-type", f.id or "?", f"{old} -> {target}"))
        for qual in vocab.feature_qualifiers.get(old, ()):
            if _is_placeholder(qual):
                changes.append(Change("needs-manual", f.id or "?",
                                      f"qualifier {qual} for {old} needs a manual value (not added)"))
                continue
            key, val = _qualifier_to_attr(qual)
            if key in f.attributes:
                continue
            f.attributes[key] = [val if val is not None else "true"]
            changes.append(Change("add-qualifier", f.id or "?",
                                  f"added {key}={f.attributes[key][0]}"))
    return changes


def _parse_pos_spec(spec: str) -> dict | None:
    """Parse '(pos:139..141,aa:Sec[,seq:ttc])' (already URL-decoded). Single range only."""
    if "join" in spec.lower():
        return None
    m = re.search(r"pos:(complement\()?(\d+)(?:\.\.(\d+))?\)?", spec)
    if not m:
        return None
    start = int(m.group(2))
    end = int(m.group(3)) if m.group(3) else start
    strand = "-" if m.group(1) else "+"
    aa_m = re.search(r"aa:([A-Za-z*]+)", spec)
    seq_m = re.search(r"seq:([A-Za-z]+)", spec)
    return {"start": start, "end": end, "strand": strand,
            "aa": aa_m.group(1) if aa_m else None,
            "seq": seq_m.group(1) if seq_m else None}


def _attach_children(doc, pending) -> None:
    for child, parent in pending:
        doc.features.append(child)
        if child.id:
            doc.feature_index[child.id] = child
        parent.children.append(child)


def pass_transl_except(doc, ctx) -> list:
    changes: list = []
    pending: list = []
    for f in list(doc.features):
        if f.type != "CDS":
            continue
        specs = f.attributes.get("transl_except")
        if not specs or not f.spans:
            continue
        lo = min(s.start for s in f.spans)
        hi = max(s.end for s in f.spans)
        seqid = f.spans[0].seqid
        kept: list = []
        made = 0
        for spec in specs:
            p = _parse_pos_spec(spec)
            if p is None or p["aa"] is None or not (lo <= p["start"] and p["end"] <= hi):
                changes.append(Change("needs-manual", f.id or "?",
                                      f"CDS {f.id!r} transl_except {spec!r} not a single in-bounds pos; kept as attribute"))
                kept.append(spec)
                continue
            made += 1
            child_id = f"{f.id}_recoded_{made}"
            if aa_names.is_stop(p["aa"]):
                ctype = "stop_codon"
                attrs = {"ID": [child_id], "Parent": [f.id], "Note": ["stop codon completed"]}
            else:
                ctype = "recoded_codon"
                attrs = {"ID": [child_id], "Parent": [f.id],
                         "codon_redefined": [aa_names.full_name(p["aa"])]}
            child = Feature(child_id, f.source, ctype,
                            [Span(seqid, p["start"], p["end"], p["strand"], 0)], attrs, [f.id])
            pending.append((child, f))
            changes.append(Change("add-child-feature", child_id,
                                  f"CDS {f.id!r}: transl_except -> {ctype} ({p['start']}..{p['end']})"))
        if made:
            if kept:
                f.attributes["transl_except"] = kept
            else:
                del f.attributes["transl_except"]
    _attach_children(doc, pending)
    return changes


def pass_anticodon(doc, ctx) -> list:
    changes: list = []
    pending: list = []
    for f in list(doc.features):
        if f.type != "tRNA":
            continue
        specs = f.attributes.get("anticodon")
        if not specs or not f.spans:
            continue
        lo = min(s.start for s in f.spans)
        hi = max(s.end for s in f.spans)
        seqid = f.spans[0].seqid
        kept: list = []
        made = 0
        for spec in specs:
            p = _parse_pos_spec(spec)
            if p is None or not (lo <= p["start"] and p["end"] <= hi):
                changes.append(Change("needs-manual", f.id or "?",
                                      f"tRNA {f.id!r} anticodon {spec!r} not a single in-bounds pos; kept as attribute"))
                kept.append(spec)
                continue
            made += 1
            child_id = f"{f.id}_anticodon_{made}"
            attrs = {"ID": [child_id], "Parent": [f.id]}
            if p["aa"]:
                attrs["amino_acid"] = [aa_names.full_name(p["aa"])]
            if p["seq"]:
                attrs["sequence"] = [p["seq"]]
            child = Feature(child_id, f.source, "anticodon",
                            [Span(seqid, p["start"], p["end"], p["strand"], None)], attrs, [f.id])
            pending.append((child, f))
            changes.append(Change("add-child-feature", child_id,
                                  f"tRNA {f.id!r}: anticodon -> anticodon child ({p['start']}..{p['end']})"))
        if made:
            if kept:
                f.attributes["anticodon"] = kept
            else:
                del f.attributes["anticodon"]
    _attach_children(doc, pending)
    return changes
