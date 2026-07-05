from ddbj_gff import parse
from ddbj_gff.model import Feature, Span, Directive, GffDocument
from ddbj_gff.normalize.passes import pass_trans_splicing_location, NormalizeContext
from ddbj_gff.normalize.normalize import normalize
from ddbj_gff.validate import validate


def _ctx():
    return NormalizeContext(vocab=None, seq_lengths=None, config=None)


def test_builds_location_for_mixed_strand_trans_cds():
    # 3 parts, part order 1,2,3, strands -,+,+
    spans = [Span("c", 1641, 1754, "-", 0, part=1),
             Span("c", 93, 324, "+", 0, part=2),
             Span("c", 829, 854, "+", 2, part=3)]
    cds = Feature("cds", "S", "CDS", spans, {"exception": ["trans-splicing"]}, [])
    seq_dir = Directive("x", "sequence-region", ("c", 1, 1754))
    doc = GffDocument(directives=[seq_dir], features=[cds])
    changes = pass_trans_splicing_location(doc, _ctx())
    assert cds.attributes["location"] == ["join(complement(1641..1754),93..324,829..854)"]
    assert len(changes) == 1


def test_two_part_trans_intron_and_qmark_strand():
    spans = [Span("c", 855, 1640, "-", None, part=1),
             Span("c", 1, 92, "?", None, part=2)]
    intron = Feature("i", "S", "intron", spans, {"exception": ["trans-splicing"]}, [])
    seq_dir = Directive("x", "sequence-region", ("c", 1, 1754))
    doc = GffDocument(directives=[seq_dir], features=[intron])
    pass_trans_splicing_location(doc, _ctx())
    assert intron.attributes["location"] == ["join(complement(855..1640),1..92)"]


def test_skips_non_trans_and_single_span_and_existing_location():
    non_trans = Feature("g", "S", "gene",
                        [Span("c", 855, 1754, "-", None, part=1), Span("c", 1, 854, "?", None, part=2)],
                        {"is_ordered": ["true"]}, [])           # multi-part but not trans-spliced
    cis = Feature("c1", "S", "intron", [Span("c", 325, 828, "+")], {"exception": ["trans-splicing"]}, [])  # single span
    preset = Feature("p", "S", "CDS",
                     [Span("c", 1641, 1754, "-", 0, part=1), Span("c", 93, 324, "+", 0, part=2)],
                     {"exception": ["trans-splicing"], "location": ["join(remote:1..9,1..2)"]}, [])
    doc = GffDocument(directives=[Directive("x", "sequence-region", ("c", 1, 1754))],
                      features=[non_trans, cis, preset])
    pass_trans_splicing_location(doc, _ctx())
    assert "location" not in non_trans.attributes
    assert "location" not in cis.attributes
    assert preset.attributes["location"] == ["join(remote:1..9,1..2)"]   # preserved


def test_full_normalize_and_validate_clean_on_fixture():
    with open("tests/normalize_fixtures/trans_splicing_rps12.gff3") as fh:
        doc = parse(fh.read())
    work, _report = normalize(doc)
    cds = next(f for f in work.features if f.type == "CDS")
    intron1 = next(f for f in work.features if f.type == "intron" and len(f.spans) == 2)
    assert cds.attributes["location"] == ["join(complement(1641..1754),93..324,829..854)"]
    assert intron1.attributes["location"] == ["join(complement(855..1640),1..92)"]
    # canonical form set -> no noncanonical-special-case for trans-splicing
    codes = {d.code for d in validate(work)}
    assert "noncanonical-special-case" not in codes
