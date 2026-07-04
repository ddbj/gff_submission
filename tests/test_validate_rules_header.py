from ddbj_gff import parse
from ddbj_gff.model import Feature, Span, Directive, GffDocument
from ddbj_gff.validate.vocab import Vocab
from ddbj_gff.validate import rules

V = Vocab(frozenset({"gene", "mRNA", "CDS", "exon"}), {}, frozenset({"GenBank"}))


def codes(diags):
    return {d.code for d in diags}


def test_directives_all_missing():
    doc = GffDocument()  # no directives
    c = codes(rules.rule_directives(doc, V))
    assert {"missing-gff-version", "missing-insdc-gff-version",
            "missing-species-taxid", "missing-sequence-region"} <= c


def test_directives_present_ok_and_duplicate_region():
    doc = GffDocument(directives=[
        Directive("##gff-version 3", "gff-version", "3"),
        Directive("#!insdc-gff-version 1.0.0", "insdc-gff-version", "1.0.0"),
        Directive("##species ...", "species", 4530),
        Directive("##sequence-region c 1 100", "sequence-region", ("c", 1, 100)),
        Directive("##sequence-region c 1 200", "sequence-region", ("c", 1, 200)),  # dup seqid
    ])
    c = codes(rules.rule_directives(doc, V))
    assert "missing-gff-version" not in c and "missing-insdc-gff-version" not in c
    assert "missing-species-taxid" not in c and "missing-sequence-region" not in c
    assert "duplicate-sequence-region" in c


def test_ascii_flags_non_ascii_attribute():
    f = Feature("g", "S", "gene", [Span("c", 1, 9, "+")], {"product": ["protéin"]}, [])
    doc = GffDocument(features=[f])
    assert "non-ascii" in codes(rules.rule_ascii(doc, V))


def test_seqid_bounds_undefined_and_outside():
    seq_dir = Directive("x", "sequence-region", ("c", 1, 100))
    f_out = Feature("a", "S", "gene", [Span("c", 50, 150, "+")], {}, [])      # 150 > 100
    f_undef = Feature("b", "S", "gene", [Span("z", 1, 9, "+")], {}, [])       # seqid z undefined
    doc = GffDocument(directives=[seq_dir], features=[f_out, f_undef])
    c = codes(rules.rule_seqid_bounds(doc, V))
    assert "feature-outside-region" in c
    assert "undefined-seqid" in c


def test_seqid_bounds_circular_origin_spanning_allowed():
    seq_dir = Directive("x", "sequence-region", ("c", 1, 100))
    f = Feature("a", "S", "CDS", [Span("c", 90, 130, "+")], {"Is_circular": ["true"]}, [])  # end>len ok if circular
    doc = GffDocument(directives=[seq_dir], features=[f])
    assert "feature-outside-region" not in codes(rules.rule_seqid_bounds(doc, V))


def test_start_gt_end():
    f = Feature("a", "S", "gene", [Span("c", 50, 10, "+")], {}, [])
    doc = GffDocument(features=[f])
    assert "start-gt-end" in codes(rules.rule_start_gt_end(doc, V))


def test_seqid_bounds_circular_landmark_allows_origin_spanning():
    region = Feature("r", "S", "region", [Span("c", 1, 100, "+")], {"Is_circular": ["true"]}, [])
    seq_dir = Directive("x", "sequence-region", ("c", 1, 100))
    cds = Feature("a", "S", "CDS", [Span("c", 90, 130, "+")], {}, [])  # flag on region, not cds
    doc = GffDocument(directives=[seq_dir], features=[region, cds])
    assert "feature-outside-region" not in codes(rules.rule_seqid_bounds(doc, V))


def test_seqid_bounds_noncircular_end_beyond_region_flagged():
    seq_dir = Directive("x", "sequence-region", ("c", 1, 100))
    cds = Feature("a", "S", "CDS", [Span("c", 90, 130, "+")], {}, [])
    doc = GffDocument(directives=[seq_dir], features=[cds])
    assert "feature-outside-region" in codes(rules.rule_seqid_bounds(doc, V))
