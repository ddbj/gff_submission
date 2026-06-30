from __future__ import annotations

from .severities import make_diagnostic


def _all_strings(feature) -> list[str]:
    out = [feature.type, feature.source]
    out += [s.seqid for s in feature.spans]
    for key, values in feature.attributes.items():
        out.append(key)
        out.extend(values)
    return out


def rule_directives(doc, vocab) -> list:
    diags = []
    if doc.gff_version is None:
        diags.append(make_diagnostic("missing-gff-version", "##gff-version directive is missing"))
    if doc.insdc_gff_version is None:
        diags.append(make_diagnostic("missing-insdc-gff-version",
                                     "#!insdc-gff-version directive is missing"))
    if not isinstance(doc.species, int):
        diags.append(make_diagnostic("missing-species-taxid",
                                     "##species directive with NCBI taxid is missing"))
    seqids = [d.value[0] for d in doc.directives
              if d.kind == "sequence-region" and d.value]
    if not seqids:
        diags.append(make_diagnostic("missing-sequence-region",
                                     "##sequence-region directive is missing"))
    seen = set()
    for sid in seqids:
        if sid in seen:
            diags.append(make_diagnostic("duplicate-sequence-region",
                                         f"duplicate ##sequence-region for seqid {sid!r}"))
        seen.add(sid)
    return diags


def rule_ascii(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        if any(not s.isascii() for s in _all_strings(f)):
            diags.append(make_diagnostic("non-ascii",
                                         f"feature {f.id!r} contains non-ASCII characters"))
    return diags


def rule_seqid_bounds(doc, vocab) -> list:
    diags = []
    regions = doc.sequence_regions
    for f in doc.features:
        circular = f.is_circular
        for s in f.spans:
            if s.seqid not in regions:
                diags.append(make_diagnostic("undefined-seqid",
                                             f"feature {f.id!r} references seqid {s.seqid!r} "
                                             f"with no ##sequence-region"))
                continue
            lo, hi = regions[s.seqid]
            if (s.start < lo or s.end > hi) and not circular:
                diags.append(make_diagnostic("feature-outside-region",
                                             f"feature {f.id!r} span {s.start}..{s.end} is outside "
                                             f"sequence-region {s.seqid}:{lo}..{hi}"))
    return diags


def rule_start_gt_end(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        for s in f.spans:
            if s.start > s.end:
                diags.append(make_diagnostic("start-gt-end",
                                             f"feature {f.id!r} has start>end ({s.start}>{s.end})"))
    return diags
