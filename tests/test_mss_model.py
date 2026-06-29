from ddbj_gff.mss.model import MssQualifier, MssFeature, MssEntry, MssDocument


def test_model_construction():
    q = MssQualifier("locus_tag", "MPTK1_000010")
    f = MssFeature("CDS", "join(1..10,20..30)", [q])
    e = MssEntry("chr1", [f])
    d = MssDocument(["COMMON\tDBLINK"], [e])
    assert f.qualifiers[0].key == "locus_tag"
    assert e.features[0].location == "join(1..10,20..30)"
    assert d.entries[0].name == "chr1"
    assert d.common_rows == ["COMMON\tDBLINK"]


def test_defaults_independent():
    a = MssEntry("a")
    b = MssEntry("b")
    a.features.append(MssFeature("source", "1..9", []))
    assert b.features == []
