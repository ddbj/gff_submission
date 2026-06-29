from __future__ import annotations

import re

from .config import MssConfig
from .model import MssFeature, MssQualifier


def assembly_gap_features(seq: str, cfg: MssConfig) -> list[MssFeature]:
    pattern = re.compile("n{%d,}" % cfg.gap_min_length)
    feats: list[MssFeature] = []
    for m in pattern.finditer(seq.lower()):
        start = m.start() + 1  # 1-based inclusive
        end = m.end()
        feats.append(MssFeature("assembly_gap", f"{start}..{end}", [
            MssQualifier("estimated_length", cfg.gap_estimated_length),
            MssQualifier("gap_type", cfg.gap_type),
            MssQualifier("linkage_evidence", cfg.gap_linkage_evidence),
        ]))
    return feats
