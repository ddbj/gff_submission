from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.gaps import assembly_gap_features


def cfg(min_length=10):
    return MssConfig(source={}, gap_min_length=min_length, gap_type="within scaffold",
                     gap_linkage_evidence="align genus", gap_estimated_length="known")


def test_detects_runs_at_or_above_min_length():
    # 10 N's at positions 6..15 (1-based)
    seq = "ACGTA" + "N" * 10 + "ACGTA"
    feats = assembly_gap_features(seq, cfg(10))
    assert len(feats) == 1
    f = feats[0]
    assert f.key == "assembly_gap"
    assert f.location == "6..15"
    assert [(q.key, q.value) for q in f.qualifiers] == [
        ("estimated_length", "known"),
        ("gap_type", "within scaffold"),
        ("linkage_evidence", "align genus"),
    ]


def test_lowercase_and_below_min_ignored():
    seq = "AAA" + "n" * 12 + "AAA" + "N" * 3 + "AAA"  # 12-run counts (lowercase), 3-run ignored
    feats = assembly_gap_features(seq, cfg(10))
    assert len(feats) == 1
    assert feats[0].location == "4..15"


def test_no_gaps():
    assert assembly_gap_features("ACGTACGT", cfg(10)) == []
