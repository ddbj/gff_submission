from ddbj_gff.normalize.config import NormalizeConfig


def test_merge_config_defaults_off():
    c = NormalizeConfig()
    assert c.merge_overlapping_loci is False
    assert c.merge_overlap_min_fraction == 0.0
