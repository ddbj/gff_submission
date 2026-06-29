from ddbj_gff.model import Feature, Span
from ddbj_gff.mss.locus_tag import LocusTagAssigner


def gene(attrs):
    return Feature("g", "S", "gene", [Span("chr1", 1, 9, "+")], attrs, [])


def test_prefers_gff_locus_tag_attribute():
    a = LocusTagAssigner("PFX", width=6, start=10, step=10)
    assert a.assign(gene({"locus_tag": ["ABC_123"]})) == "ABC_123"


def test_sequential_fallback_increments():
    a = LocusTagAssigner("PFX", width=6, start=10, step=10)
    assert a.assign(gene({})) == "PFX_000010"
    assert a.assign(gene({})) == "PFX_000020"
    # an attribute-bearing gene does not consume a number
    assert a.assign(gene({"locus_tag": ["X_1"]})) == "X_1"
    assert a.assign(gene({})) == "PFX_000030"
