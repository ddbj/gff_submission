from ddbj_gff.model import Feature, Span
from ddbj_gff.mss.convert import build_mrna_feature, mrna_partial_flags


def build_gene_mrna(exon_spans, cds_spans, strand="+", gene_attr=None):
    gene = Feature("Mp1g1", "S", "gene", [Span("chr1", 1, 999, strand)], gene_attr or {}, [])
    mrna = Feature("Mp1g1.1", "S", "mRNA", [Span("chr1", 1, 999, strand)], {}, [])
    children = []
    for i, (s, e) in enumerate(exon_spans):
        children.append(Feature(f"ex{i}", "S", "exon", [Span("chr1", s, e, strand)], {}, []))
    for i, (s, e) in enumerate(cds_spans):
        children.append(Feature(f"cds{i}", "S", "CDS", [Span("chr1", s, e, strand)], {}, []))
    mrna.children = children
    return gene, mrna


def test_mrna_complete_when_utr_present_plus():
    # exon extends beyond CDS at both ends -> complete (no partial)
    gene, mrna = build_gene_mrna([(1, 30), (40, 70)], [(11, 30), (40, 60)], strand="+")
    assert mrna_partial_flags(mrna) == (False, False)
    f = build_mrna_feature(mrna, gene, "PFX_000010", 10000)
    assert f.key == "mRNA"
    assert f.location == "join(1..30,40..70)"
    assert ("locus_tag", "PFX_000010") in [(q.key, q.value) for q in f.qualifiers]
    assert any(q.key == "note" and "submitter_gene_id: Mp1g1" in q.value
               and "submitter_transcript_id: Mp1g1.1" in q.value for q in f.qualifiers)


def test_mrna_partial_when_no_utr_plus():
    # exon start == CDS start (no 5'UTR) and exon end == CDS end (no 3'UTR)
    gene, mrna = build_gene_mrna([(1, 30), (40, 70)], [(1, 30), (40, 70)], strand="+")
    assert mrna_partial_flags(mrna) == (True, True)
    f = build_mrna_feature(mrna, gene, "PFX_000010", 10000)
    assert f.location == "join(<1..30,40..>70)"


def test_mrna_gene_qualifier_when_present():
    gene, mrna = build_gene_mrna([(1, 30)], [(11, 20)], strand="+", gene_attr={"gene": ["MpX"]})
    f = build_mrna_feature(mrna, gene, "PFX_000010", 10000)
    assert ("gene", "MpX") in [(q.key, q.value) for q in f.qualifiers]


def test_mrna_partial_minus_strand_swap():
    # minus strand: genomic-left (low coord) lacks UTR -> that is the biological 3' end.
    # exon [(1,30),(40,70)] vs CDS [(1,30),(40,60)]: left_partial=True (exon_lo==cds_lo),
    # right_partial=False (exon_hi 70 > cds_hi 60). On '-' strand -> (right, left) = (False, True).
    gene, mrna = build_gene_mrna([(1, 30), (40, 70)], [(1, 30), (40, 60)], strand="-")
    assert mrna_partial_flags(mrna) == (False, True)
    f = build_mrna_feature(mrna, gene, "PFX_000010", 10000)
    assert f.location == "complement(join(<1..30,40..70))"
