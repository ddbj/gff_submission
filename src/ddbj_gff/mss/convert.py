from __future__ import annotations

import re

from Bio.SeqFeature import (AfterPosition, BeforePosition, CompoundLocation,
                            FeatureLocation)
from Bio.SeqIO.InsdcIO import _insdc_location_string

from .config import MssConfig
from .model import MssFeature, MssQualifier

_STRAND = {"+": 1, "-": -1}


def collect_spans(parent, feature_type: str) -> list:
    spans = []
    for child in parent.children:
        if child.type == feature_type:
            spans.extend(child.spans)
    return spans


def _ordered(spans: list):
    strand = spans[0].strand
    return sorted(spans, key=lambda s: s.start, reverse=(strand == "-"))


def build_insdc_location(spans: list, seqlen: int,
                         five_prime_partial: bool = False,
                         three_prime_partial: bool = False) -> str:
    bio = _STRAND.get(spans[0].strand, 0)
    ordered = _ordered(spans)
    locs = [FeatureLocation(s.start - 1, s.end, strand=bio) for s in ordered]
    if five_prime_partial:
        f = locs[0]
        if bio == -1:
            locs[0] = FeatureLocation(f.start, AfterPosition(int(f.end)), strand=bio)
        else:
            locs[0] = FeatureLocation(BeforePosition(int(f.start)), f.end, strand=bio)
    if three_prime_partial:
        l = locs[-1]
        if bio == -1:
            locs[-1] = FeatureLocation(BeforePosition(int(l.start)), l.end, strand=bio)
        else:
            locs[-1] = FeatureLocation(l.start, AfterPosition(int(l.end)), strand=bio)
    compound = locs[0] if len(locs) == 1 else CompoundLocation(locs)
    return _insdc_location_string(compound, seqlen)


def extract_seq(spans: list, genome_seq):
    bio = _STRAND.get(spans[0].strand, 0)
    ordered = _ordered(spans)
    locs = [FeatureLocation(s.start - 1, s.end, strand=bio) for s in ordered]
    compound = locs[0] if len(locs) == 1 else CompoundLocation(locs)
    return compound.extract(genome_seq)


def build_source_feature(seqid: str, seqlen: int, cfg: MssConfig) -> MssFeature:
    quals = [MssQualifier(k, v) for k, v in cfg.source.items()]
    chromosome = None
    if cfg.chromosome_pattern:
        m = re.match(cfg.chromosome_pattern, seqid)
        if m:
            chromosome = m.group(1)
    if chromosome is not None:
        quals.append(MssQualifier("chromosome", chromosome))
        if cfg.ff_def_chromosome:
            quals.append(MssQualifier("ff_definition", cfg.ff_def_chromosome))
    else:
        quals.append(MssQualifier("submitter_seqid", seqid))
        if cfg.ff_def_default:
            quals.append(MssQualifier("ff_definition", cfg.ff_def_default))
    return MssFeature("source", f"1..{seqlen}", quals)


def _submitter_note(gene, mrna) -> MssQualifier:
    return MssQualifier(
        "note", f"submitter_gene_id: {gene.id}, submitter_transcript_id: {mrna.id}")


def mrna_partial_flags(mrna) -> tuple[bool, bool]:
    exon = collect_spans(mrna, "exon")
    cds = collect_spans(mrna, "CDS")
    if not exon or not cds:
        return (False, False)
    exon_lo, exon_hi = min(s.start for s in exon), max(s.end for s in exon)
    cds_lo, cds_hi = min(s.start for s in cds), max(s.end for s in cds)
    strand = exon[0].strand
    # no UTR on a genomic side -> that side is partial; map to 5'/3' by strand
    left_partial = exon_lo == cds_lo
    right_partial = exon_hi == cds_hi
    if strand == "-":
        return (right_partial, left_partial)   # 5' = genomic right on minus strand
    return (left_partial, right_partial)


def build_mrna_feature(mrna, gene, locus_tag: str, seqlen: int) -> MssFeature:
    spans = collect_spans(mrna, "exon") or collect_spans(mrna, "CDS")
    fp, tp = mrna_partial_flags(mrna)
    location = build_insdc_location(spans, seqlen, fp, tp)
    quals = [MssQualifier("locus_tag", locus_tag)]
    if gene.gene:
        quals.append(MssQualifier("gene", gene.gene))
    quals.append(_submitter_note(gene, mrna))
    return MssFeature("mRNA", location, quals)
