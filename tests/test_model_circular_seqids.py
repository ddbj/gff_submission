from ddbj_gff.model import Feature, Span, GffDocument


def test_circular_seqids_from_region_landmark():
    region = Feature("r", "S", "region", [Span("CP", 1, 100, "+")], {"Is_circular": ["true"]}, [])
    gene = Feature("g", "S", "gene", [Span("CP", 90, 130, "+")], {}, [])
    linear = Feature("l", "S", "region", [Span("MT", 1, 50, "+")], {}, [])
    doc = GffDocument(features=[region, gene, linear])
    assert doc.circular_seqids == {"CP"}


def test_circular_seqids_empty_when_no_landmark_flag():
    gene = Feature("g", "S", "gene", [Span("CP", 1, 9, "+")], {"Is_circular": ["true"]}, [])
    # flag on a non-landmark feature type does not make the seqid circular
    doc = GffDocument(features=[gene])
    assert doc.circular_seqids == set()
