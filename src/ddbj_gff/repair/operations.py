"""Concrete repair operations, registered on import.

Each operation is a `detect(doc, ctx) -> list[Candidate]` (non-destructive scan)
paired with an `apply(doc, ctx, selection) -> list[Change]` (mutating; `selection`
is either a caller-provided subset of candidates or `None` to re-detect).
"""
from __future__ import annotations

from ..normalize.report import Change
from .cds import coding_sequence
from .context import RepairContext
from .registry import Operation, register
from .report import Candidate
from .partial import is_partial, partial_attrs


def _child_spans(mrna, ftype):
    return [s for c in mrna.children if c.type == ftype for s in c.spans]


# --- utr-absent-to-partial-mrna -------------------------------------------

def _detect_utr(doc, ctx: RepairContext) -> list[Candidate]:
    out: list[Candidate] = []
    for f in doc.features:
        if f.type != "mRNA" or is_partial(f):
            continue
        exon = _child_spans(f, "exon")
        cds = _child_spans(f, "CDS")
        if not exon or not cds:
            continue
        exon_lo, exon_hi = min(s.start for s in exon), max(s.end for s in exon)
        cds_lo, cds_hi = min(s.start for s in cds), max(s.end for s in cds)
        strand = exon[0].strand
        left_partial = exon_lo == cds_lo
        right_partial = exon_hi == cds_hi
        five, three = ((right_partial, left_partial) if strand == "-"
                       else (left_partial, right_partial))
        if not (five or three):
            continue
        m_lo = min(s.start for s in f.spans) if f.spans else exon_lo
        m_hi = max(s.end for s in f.spans) if f.spans else exon_hi
        ends = ", ".join(e for e, on in (("5'", five), ("3'", three)) if on)
        out.append(Candidate(
            "utr-absent-to-partial-mrna", f.id, exon[0].seqid,
            detail=f"mRNA {f.id!r} missing UTR on {ends} end -> mark partial",
            payload={"five": five, "three": three, "strand": strand,
                     "start": m_lo, "end": m_hi}))
    return out


def _apply_utr(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_utr(doc, ctx)
    changes: list[Change] = []
    for c in cands:
        f = doc.feature_index.get(c.feature_id)
        if f is None or f.type != "mRNA" or is_partial(f):
            continue
        p = c.payload
        attrs = partial_attrs(p["five"], p["three"], p["strand"], p["start"], p["end"])
        f.attributes.update(attrs)
        changes.append(Change("mark-partial", f.id or "?",
                              f"mRNA marked partial ({', '.join(sorted(attrs))})"))
    return changes


register(Operation("utr-absent-to-partial-mrna",
                   "Mark an mRNA partial on ends where a UTR is absent (structural).",
                   requires_sequence=False, detect=_detect_utr, apply=_apply_utr))


# --- missing-start-stop-to-partial-cds ------------------------------------

def _detect_startstop(doc, ctx: RepairContext) -> list[Candidate]:
    out: list[Candidate] = []
    if ctx.sequences is None:
        return out
    for f in doc.features:
        if f.type != "CDS" or is_partial(f) or not f.spans:
            continue
        if f.is_trans_spliced:
            continue
        seqid = f.spans[0].seqid
        if seqid not in ctx.sequences:
            continue
        coding, table = coding_sequence(f, ctx)
        codon_start = f.codon_start or 1
        strand = f.spans[0].strand
        first = coding[:3]
        last = coding[-3:]
        five = codon_start > 1 or len(first) < 3 or first not in table.start_codons
        three = len(last) < 3 or last not in table.stop_codons
        if not (five or three):
            continue
        c_lo = min(s.start for s in f.spans)
        c_hi = max(s.end for s in f.spans)
        ends = ", ".join(e for e, on in (("5' (no start codon)", five),
                                         ("3' (no stop codon)", three)) if on)
        out.append(Candidate(
            "missing-start-stop-to-partial-cds", f.id, seqid,
            detail=f"CDS {f.id!r} partial: {ends}",
            payload={"five": five, "three": three, "strand": strand,
                     "start": c_lo, "end": c_hi}))
    return out


def _apply_startstop(doc, ctx: RepairContext, selection):
    cands = selection if selection is not None else _detect_startstop(doc, ctx)
    changes: list[Change] = []
    for c in cands:
        f = doc.feature_index.get(c.feature_id)
        if f is None or f.type != "CDS" or is_partial(f):
            continue
        p = c.payload
        attrs = partial_attrs(p["five"], p["three"], p["strand"], p["start"], p["end"])
        f.attributes.update(attrs)
        changes.append(Change("mark-partial", f.id or "?",
                              f"CDS marked partial ({', '.join(sorted(attrs))})"))
    return changes


register(Operation("missing-start-stop-to-partial-cds",
                   "Mark a CDS partial when its sequence lacks a start or stop codon.",
                   requires_sequence=True, detect=_detect_startstop, apply=_apply_startstop))
