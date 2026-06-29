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
