from __future__ import annotations

import re
from dataclasses import dataclass

from ..model import Directive, Feature, Span
from .. import aa_names
from .report import Change
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.SeqIO.InsdcIO import _insdc_location_string

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


def pass_coerce_transcript_to_mrna(doc, ctx) -> list:
    changes: list = []
    if not getattr(ctx.config, "coerce_transcript_to_mrna", True):
        return changes
    for f in doc.features:
        if f.type != "transcript":
            continue
        if any(c.type == "CDS" for c in f.children):
            f.type = "mRNA"
            changes.append(Change("rename-type", f.id or "?", "transcript -> mRNA (has CDS)"))
    return changes


def pass_wrap_cds_in_mrna(doc, ctx) -> list:
    """Insert an mRNA level for genes whose CDS/exon hang directly off the gene.

    Common in organelle LiftOff output (gene -> CDS). tRNA/ncRNA genes
    (gene -> tRNA/ncRNA -> exon) already have a transcript level and are left alone.
    """
    changes: list = []
    if not getattr(ctx.config, "wrap_cds_in_mrna", True):
        return changes
    struct = {"CDS", "exon"}
    new_mrnas: list = []
    for gene in list(doc.features):
        if gene.type != "gene":
            continue
        if any(c.type in ("mRNA", "transcript") for c in gene.children):
            continue
        wrapped = [c for c in gene.children if c.type in struct]
        if not any(c.type == "CDS" for c in wrapped):
            continue
        seqid = wrapped[0].spans[0].seqid
        strand = wrapped[0].spans[0].strand
        lo = min(s.start for c in wrapped for s in c.spans)
        hi = max(s.end for c in wrapped for s in c.spans)
        mid = f"{gene.id}.mrna"
        mrna = Feature(mid, gene.source, "mRNA", [Span(seqid, lo, hi, strand)],
                       {"ID": [mid], "Parent": [gene.id]}, [gene.id])
        for c in wrapped:
            c.parent_ids = [mid]
            c.parents = [mrna]
            if "Parent" in c.attributes:
                c.attributes["Parent"] = [mid]
            mrna.children.append(c)
        mrna.parents = [gene]
        gene.children = [c for c in gene.children if c.type not in struct] + [mrna]
        new_mrnas.append(mrna)
        changes.append(Change("add-child-feature", mid,
                              f"wrapped {len(wrapped)} CDS/exon of {gene.id!r} in an mRNA"))
    for m in new_mrnas:
        doc.features.append(m)
        if m.id:
            doc.feature_index[m.id] = m
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


def pass_circular_origin(doc, ctx) -> list:
    """Propagate Is_circular=true onto origin-spanning features (a span with
    end>seqlen) on a circular landmark. Coordinates are left as-is (canonical keeps
    the INSDC end>seqlen convention); the flag lets validate/downstream treat the
    feature as circular."""
    changes: list = []
    circular = doc.circular_seqids
    if not circular:
        return changes
    regions = doc.sequence_regions
    for f in doc.features:
        for s in f.spans:
            if s.seqid not in circular:
                continue
            seqlen = regions.get(s.seqid, (None, None))[1]
            if seqlen is None and ctx.seq_lengths:
                seqlen = ctx.seq_lengths.get(s.seqid)
            if seqlen is not None and s.end > seqlen:
                if f.attributes.get("Is_circular") != ["true"]:
                    f.attributes["Is_circular"] = ["true"]
                    changes.append(Change("add-qualifier", f.id or "?",
                                          f"propagated Is_circular=true to origin-spanning "
                                          f"feature (span {s.start}..{s.end} > seqlen {seqlen})"))
                break
    return changes


def _find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def pass_merge_overlapping_loci(doc, ctx) -> list:
    """Merge same-strand gene loci whose mRNAs overlap into one gene (opt-in).

    Edge between two mRNAs when overlap_bp/min(len) >= config.merge_overlap_min_fraction
    (default 0.0 = any overlap). Connected components form a locus; the gene of the
    lowest-coordinate mRNA is the representative, the others' mRNAs are reparented to it,
    and its span becomes the union. Trans-spliced transcripts (mRNA or a CDS child with
    exception=trans-splicing) are excluded entirely."""
    changes: list = []
    if not getattr(ctx.config, "merge_overlapping_loci", False):
        return changes
    frac = getattr(ctx.config, "merge_overlap_min_fraction", 0.0)

    def _is_trans(m):
        return m.is_trans_spliced or any(
            c.type == "CDS" and c.is_trans_spliced for c in m.children)

    def _gene_of(m):
        pid = m.parent_ids[0] if m.parent_ids else None
        return doc.feature_index.get(pid) if pid else None

    groups: dict = {}
    for f in doc.features:
        if f.type != "mRNA" or not f.spans or _is_trans(f):
            continue
        lo = min(s.start for s in f.spans)
        hi = max(s.end for s in f.spans)
        groups.setdefault((f.spans[0].seqid, f.spans[0].strand), []).append((lo, hi, f))

    touched: set = set()
    for (seqid, strand), items in groups.items():
        items.sort(key=lambda t: (t[0], t[1], t[2].id or ""))
        n = len(items)
        parent = list(range(n))
        for i in range(n):
            lo_i, hi_i, _ = items[i]
            for j in range(i + 1, n):
                lo_j, hi_j, _ = items[j]
                if lo_j > hi_i:
                    break
                ov = min(hi_i, hi_j) - max(lo_i, lo_j) + 1
                if ov > 0 and ov / min(hi_i - lo_i + 1, hi_j - lo_j + 1) >= frac:
                    ri, rj = _find(parent, i), _find(parent, j)
                    if ri != rj:
                        parent[rj] = ri
        comps: dict = {}
        for i in range(n):
            comps.setdefault(_find(parent, i), []).append(items[i][2])
        for members in comps.values():
            genes = {}
            for m in members:
                g = _gene_of(m)
                if g is not None:
                    genes[g.id] = g
            if len(genes) < 2:
                continue
            members.sort(key=lambda m: (min(s.start for s in m.spans),
                                        max(s.end for s in m.spans), m.id or ""))
            rep = next((_gene_of(m) for m in members if _gene_of(m) is not None), None)
            if rep is None:
                continue
            u_lo = min(min(s.start for s in m.spans) for m in members)
            u_hi = max(max(s.end for s in m.spans) for m in members)
            for m in members:
                g = _gene_of(m)
                if g is rep:
                    continue
                if g is not None:
                    g.children = [c for c in g.children if c is not m]
                    touched.add(g.id)
                m.parent_ids = [rep.id]
                m.attributes["Parent"] = [rep.id]
                m.parents = [rep]
                rep.children.append(m)
            rep.spans = [Span(seqid, u_lo, u_hi, strand)]
            changes.append(Change("merge-loci", rep.id or "?",
                                  f"merged {len(genes)} loci into {rep.id!r} "
                                  f"({len(members)} mRNAs, {seqid}:{u_lo}..{u_hi})"))
    dead: set = set()
    for gid in touched:
        g = doc.feature_index.get(gid)
        if g is None:
            continue
        if not g.children:
            dead.add(gid)
        else:                                  # survived with remaining children -> recompute span
            cs = [s for c in g.children for s in c.spans]
            if cs:
                strand0 = g.spans[0].strand if g.spans else cs[0].strand
                g.spans = [Span(cs[0].seqid, min(s.start for s in cs),
                                max(s.end for s in cs), strand0)]
    if dead:
        doc.features = [f for f in doc.features if not (f.type == "gene" and f.id in dead)]
        doc.roots = [f for f in doc.roots if not (f.type == "gene" and f.id in dead)]
        for gid in dead:
            doc.feature_index.pop(gid, None)
    return changes


def pass_trans_splicing_location(doc, ctx) -> list:
    """Build the canonical INSDC location=join(...) attribute for trans-spliced
    multi-part features. Per-part strand: '-' -> complement, '+'/'?'/'.' -> forward;
    parts kept in feature.ordered_spans() order (honors is_ordered). An existing
    location= is authoritative (may be remote) and is left untouched."""
    changes: list = []
    regions = doc.sequence_regions
    for f in doc.features:
        if not f.is_trans_spliced or len(f.spans) < 2 or f.attributes.get("location"):
            continue
        parts = f.ordered_spans()
        seqlen = regions.get(parts[0].seqid, (None, None))[1]
        if seqlen is None and ctx.seq_lengths:
            seqlen = ctx.seq_lengths.get(parts[0].seqid)
        seqlen = seqlen or max(s.end for s in parts)
        locs = [FeatureLocation(s.start - 1, s.end, strand=(-1 if s.strand == "-" else 1))
                for s in parts]
        loc_str = _insdc_location_string(CompoundLocation(locs), seqlen)
        f.attributes["location"] = [loc_str]
        changes.append(Change("add-qualifier", f.id or "?",
                              f"built location={loc_str} for trans-spliced feature"))
    return changes
