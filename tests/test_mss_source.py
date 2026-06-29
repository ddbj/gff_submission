from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import build_source_feature


def cfg():
    return MssConfig(
        source={"organism": "Foo bar", "mol_type": "genomic DNA"},
        chromosome_pattern="^chr(.+)$",
        ff_def_chromosome="@@[organism]@@ chromosome: @@[chromosome]@@",
        ff_def_default="@@[organism]@@ @@[entry]@@",
    )


def test_source_chromosome_match():
    f = build_source_feature("chr1", 1000, cfg())
    assert f.key == "source" and f.location == "1..1000"
    quals = [(q.key, q.value) for q in f.qualifiers]
    assert ("organism", "Foo bar") in quals
    assert ("mol_type", "genomic DNA") in quals
    assert ("chromosome", "1") in quals
    assert ("ff_definition", "@@[organism]@@ chromosome: @@[chromosome]@@") in quals
    assert quals[0] == ("organism", "Foo bar")  # config order preserved, organism first


def test_source_non_chromosome_uses_submitter_seqid():
    f = build_source_feature("scaffold_7", 500, cfg())
    quals = [(q.key, q.value) for q in f.qualifiers]
    assert ("submitter_seqid", "scaffold_7") in quals
    assert ("ff_definition", "@@[organism]@@ @@[entry]@@") in quals
    assert not any(k == "chromosome" for k, _ in quals)
