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


def rule_feature_type(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        if f.type not in vocab.feature_types:
            diags.append(make_diagnostic("feature-type-not-insdc",
                                         f"feature {f.id!r} type {f.type!r} is not an INSDC-supported "
                                         f"SO term"))
    return diags


def rule_parents(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        if len(f.parent_ids) > 1:
            diags.append(make_diagnostic("multiple-parents",
                                         f"feature {f.id!r} has {len(f.parent_ids)} parents "
                                         f"(INSDC allows a single parent per row)"))
        for pid in f.parent_ids:
            if pid not in doc.feature_index:
                diags.append(make_diagnostic("dangling-parent",
                                             f"feature {f.id!r} references missing Parent {pid!r}"))
    return diags


def rule_cds(doc, vocab) -> list:
    diags = []
    has_file_table = doc.transl_table_map is not None
    for f in doc.features:
        if f.type != "CDS":
            continue
        if f.transl_table is None and not has_file_table:
            diags.append(make_diagnostic("cds-missing-transl-table",
                                         f"CDS {f.id!r} lacks transl_table and no file-level "
                                         f"#!transl_table is present"))
        for s in f.spans:
            if s.phase not in (0, 1, 2):
                diags.append(make_diagnostic("cds-invalid-phase",
                                             f"CDS {f.id!r} has invalid phase {s.phase!r}"))
    return diags


def rule_gene_locus_tag(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        if f.type == "gene" and not f.locus_tag:
            diags.append(make_diagnostic("gene-missing-locus-tag",
                                         f"gene {f.id!r} has no locus_tag"))
    return diags


def rule_dbxref(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        for xref in f.dbxref:
            dbtag = xref.split(":", 1)[0]
            if dbtag and dbtag not in vocab.dbxref_dbtags:
                diags.append(make_diagnostic("dbxref-unknown-dbtag",
                                             f"feature {f.id!r} Dbxref DBTAG {dbtag!r} is not in the "
                                             f"INSDC vocabulary"))
    return diags


def rule_special_case(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        if f.is_trans_spliced and "location" not in f.attributes:
            diags.append(make_diagnostic("noncanonical-special-case",
                                         f"feature {f.id!r} uses non-canonical trans-splicing "
                                         f"representation (no location= attribute)"))
        if "transl_except" in f.attributes:
            diags.append(make_diagnostic("noncanonical-special-case",
                                         f"feature {f.id!r} uses transl_except attribute "
                                         f"(canonical form is a recoded_codon child feature)"))
        if "anticodon" in f.attributes:
            diags.append(make_diagnostic("noncanonical-special-case",
                                         f"feature {f.id!r} uses anticodon attribute "
                                         f"(canonical form is an anticodon child feature)"))
    return diags


ALL_RULES = [
    rule_directives,
    rule_ascii,
    rule_seqid_bounds,
    rule_start_gt_end,
    rule_feature_type,
    rule_parents,
    rule_cds,
    rule_gene_locus_tag,
    rule_dbxref,
    rule_special_case,
]
