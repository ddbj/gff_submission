from __future__ import annotations

from ..model import Span

_STRAND = {1: "+", -1: "-", 0: "?", None: "."}

# INSDC qualifier -> canonical GFF attribute key
_QUAL_MAP = {
    "gene": "gene", "locus_tag": "locus_tag", "product": "product",
    "note": "Note", "protein_id": "protein_id", "gene_synonym": "gene_synonym",
    "transl_table": "transl_table", "db_xref": "Dbxref", "pseudo": "pseudo",
    "ncRNA_class": "ncRNA_class",
}
_DROP_QUALS = {"translation", "codon_start"}


def bio_location_to_spans(location, seqid, *, is_cds, codon_start=1):
    """BioPython location -> 1-based Spans in biological (5'->3') order.
    parts come from BioPython in transcription order. For CDS, per-segment phase
    is derived from codon_start; else phase is None."""
    parts = list(location.parts)
    spans = []
    phase = (codon_start - 1) if is_cds else None
    for p in parts:
        start = int(p.start) + 1
        end = int(p.end)
        strand = _STRAND.get(p.strand, ".")
        spans.append(Span(seqid, start, end, strand, phase=phase))
        if is_cds:
            seg_len = end - start + 1
            phase = (3 - ((seg_len - phase) % 3)) % 3
    return spans


def qualifiers_to_attrs(feature) -> dict:
    attrs: dict[str, list[str]] = {}
    for k, vals in feature.qualifiers.items():
        if k in _DROP_QUALS:
            continue
        gk = _QUAL_MAP.get(k, k)
        attrs[gk] = list(vals)
    return attrs


from collections import OrderedDict
from ..model import Feature

_RNA_TYPES = {"mRNA", "tRNA", "rRNA", "ncRNA", "misc_RNA"}
_BIOTYPE = {"CDS": "protein_coding", "mRNA": "protein_coding",
            "tRNA": "tRNA", "rRNA": "rRNA", "ncRNA": "ncRNA"}


def _locus_tag(f):
    return f.qualifiers.get("locus_tag", [None])[0]


def _cds_within(cds_spans, mrna_spans) -> bool:
    for cs in cds_spans:
        if not any(ms.strand == cs.strand and ms.start <= cs.start and cs.end <= ms.end
                   for ms in mrna_spans):
            return False
    return True


def _shared_boundaries(cds_spans, mrna_spans) -> int:
    edges = {(s.start) for s in mrna_spans} | {(s.end) for s in mrna_spans}
    return sum((s.start in edges) + (s.end in edges) for s in cds_spans)


def synthesize_features(rec, seqid) -> list:
    """Group flatfile mRNA/CDS/RNA by locus_tag into canonical gene->mRNA->exon/CDS."""
    bio = [f for f in rec.features if f.type in _RNA_TYPES or f.type == "CDS"]
    groups = OrderedDict()
    for f in bio:
        groups.setdefault(_locus_tag(f), []).append(f)

    out: list = []
    for lt, members in groups.items():
        gene_id = f"gene-{lt}" if lt else f"gene-{len(out)}"
        mrnas = [f for f in members if f.type == "mRNA"]
        cdss = [f for f in members if f.type == "CDS"]
        rnas = [f for f in members if f.type in _RNA_TYPES and f.type != "mRNA"]

        # pair each CDS to a containing mRNA; synthesize an mRNA if none
        transcripts = []  # list of (mrna_feature_or_None_synth, [cds])
        used = {id(m): [] for m in mrnas}
        for c in cdss:
            cspans = bio_location_to_spans(c.location, seqid, is_cds=True)
            cand = [m for m in mrnas
                    if _cds_within(cspans, bio_location_to_spans(m.location, seqid, is_cds=False))]
            if cand:
                best = max(cand, key=lambda m: _shared_boundaries(
                    cspans, bio_location_to_spans(m.location, seqid, is_cds=False)))
                used[id(best)].append(c)
            else:
                transcripts.append((None, [c]))   # synth mRNA = CDS
        for m in mrnas:
            transcripts.append((m, used[id(m)]))
        for r in rnas:                              # tRNA/rRNA: transcript is the RNA itself
            transcripts.append((r, []))

        # gene span = union of all member spans
        all_spans = []
        for m in members:
            all_spans += bio_location_to_spans(m.location, seqid, is_cds=(m.type == "CDS"))
        g_lo, g_hi = min(s.start for s in all_spans), max(s.end for s in all_spans)
        g_strand = all_spans[0].strand
        biotype = _BIOTYPE.get(members[0].type, "other")
        gene_attrs = {"ID": [gene_id]}
        if _locus_tag(members[0]):
            gene_attrs["locus_tag"] = [_locus_tag(members[0])]
        if members[0].qualifiers.get("gene"):
            gene_attrs["gene"] = list(members[0].qualifiers["gene"])
            gene_attrs["Name"] = list(members[0].qualifiers["gene"])
        if members[0].qualifiers.get("gene_synonym"):
            gene_attrs["gene_synonym"] = list(members[0].qualifiers["gene_synonym"])
        gene_attrs["gene_biotype"] = [biotype]
        out.append(Feature(gene_id, "DDBJ", "gene", [Span(seqid, g_lo, g_hi, g_strand)],
                           gene_attrs, []))

        for i, (mfeat, member_cds) in enumerate(transcripts, 1):
            is_rna = mfeat is not None and mfeat.type in _RNA_TYPES and mfeat.type != "mRNA"
            tx_type = mfeat.type if is_rna else "mRNA"
            tx_id = f"{('rna' if is_rna else 'mrna')}-{lt}-{i}"
            if mfeat is not None:
                tx_spans = bio_location_to_spans(mfeat.location, seqid, is_cds=False)
                tx_attrs = qualifiers_to_attrs(mfeat)
            else:                                   # synth mRNA from the CDS
                cspans = bio_location_to_spans(member_cds[0].location, seqid, is_cds=True)
                tx_spans = [Span(seqid, s.start, s.end, s.strand) for s in cspans]
                tx_attrs = {k: v for k, v in qualifiers_to_attrs(member_cds[0]).items()
                            if k in ("locus_tag", "gene", "product", "Note")}
            tx_attrs["ID"] = [tx_id]
            tx_attrs["Parent"] = [gene_id]
            out.append(Feature(tx_id, "DDBJ", tx_type, tx_spans, tx_attrs, [gene_id]))

            if not is_rna:                          # exons for mRNA transcripts
                for j, sp in enumerate(tx_spans, 1):
                    ex_id = f"exon-{lt}-{i}-{j}"
                    out.append(Feature(ex_id, "DDBJ", "exon",
                                       [Span(seqid, sp.start, sp.end, sp.strand)],
                                       {"ID": [ex_id], "Parent": [tx_id]}, [tx_id]))
            for c in member_cds:
                cspans = bio_location_to_spans(c.location, seqid, is_cds=True,
                    codon_start=int(c.qualifiers.get("codon_start", ["1"])[0]))
                c_attrs = qualifiers_to_attrs(c)
                c_id = f"cds-{c.qualifiers.get('protein_id', [tx_id])[0]}"
                c_attrs["ID"] = [c_id]
                c_attrs["Parent"] = [tx_id]
                out.append(Feature(c_id, "DDBJ", "CDS", cspans, c_attrs, [tx_id]))
    return out
