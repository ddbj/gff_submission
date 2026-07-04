from ddbj_gff import parse
from ddbj_gff.model import Feature, Span, Directive, GffDocument
from ddbj_gff.normalize.passes import pass_circular_origin, NormalizeContext
from ddbj_gff.normalize.normalize import normalize


def _ctx():
    return NormalizeContext(vocab=None, seq_lengths=None, config=None)


def test_pass_flags_origin_spanning_feature_on_circular_seqid():
    region = Feature("r", "S", "region", [Span("CP", 1, 100, "+")], {"Is_circular": ["true"]}, [])
    wrap = Feature("w", "S", "gene", [Span("CP", 90, 130, "+")], {}, [])   # end 130 > 100
    inside = Feature("i", "S", "gene", [Span("CP", 10, 40, "+")], {}, [])  # within bounds
    seq_dir = Directive("x", "sequence-region", ("CP", 1, 100))
    doc = GffDocument(directives=[seq_dir], features=[region, wrap, inside])
    changes = pass_circular_origin(doc, _ctx())
    assert wrap.attributes.get("Is_circular") == ["true"]
    assert "Is_circular" not in inside.attributes
    assert len(changes) == 1


def test_pass_noop_when_seqid_not_circular():
    wrap = Feature("w", "S", "gene", [Span("MT", 90, 130, "+")], {}, [])
    seq_dir = Directive("x", "sequence-region", ("MT", 1, 100))
    doc = GffDocument(directives=[seq_dir], features=[wrap])
    assert pass_circular_origin(doc, _ctx()) == []
    assert "Is_circular" not in wrap.attributes


def test_full_normalize_flags_moda_on_cp187952():
    with open("tests/normalize_fixtures/cp187952_origin.gff3") as fh:
        doc = parse(fh.read())
    work, _report = normalize(doc)
    moda = [f for f in work.features
            if f.type in ("gene", "CDS") and f._first("locus_tag") == "ACPZ3T_00005"]
    assert moda, "modA gene/CDS not found after normalize"
    assert all(f.is_circular for f in moda)
