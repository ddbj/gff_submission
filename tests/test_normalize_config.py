from ddbj_gff.normalize.config import NormalizeConfig, load_normalize_config


def test_defaults():
    c = NormalizeConfig()
    assert c.taxid is None
    assert c.transl_table == 1
    assert c.insdc_gff_version == "1.0.0"


def test_load_from_toml(tmp_path):
    p = tmp_path / "n.toml"
    p.write_text('[normalize]\ntaxid = 3702\ntransl_table = 11\ninsdc_gff_version = "1.0.0"\n')
    c = load_normalize_config(str(p))
    assert c.taxid == 3702
    assert c.transl_table == 11


def test_load_missing_section_uses_defaults(tmp_path):
    p = tmp_path / "n.toml"
    p.write_text('[other]\nx = 1\n')
    c = load_normalize_config(str(p))
    assert c.taxid is None and c.transl_table == 1
