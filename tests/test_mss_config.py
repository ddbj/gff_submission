import pytest
from ddbj_gff.errors import Severity, GffParseError
from ddbj_gff.mss.config import load_config, load_common

CONFIG = """\
[source]
organism = "Foo bar"
mol_type = "genomic DNA"
[source.chromosome]
pattern = "^chr(.+)$"
[source.ff_definition]
chromosome = "@@[organism]@@ chromosome: @@[chromosome]@@"
default = "@@[organism]@@ @@[entry]@@"
[locus_tag]
prefix = "PFX"
[unknown_section]
x = 1
"""


def test_load_config_defaults_and_unknown_key(tmp_path):
    p = tmp_path / "c.toml"
    p.write_text(CONFIG)
    cfg, diags = load_config(str(p))
    assert cfg.source == {"organism": "Foo bar", "mol_type": "genomic DNA"}
    assert cfg.chromosome_pattern == "^chr(.+)$"
    assert cfg.locus_tag_prefix == "PFX"
    assert cfg.locus_tag_width == 6 and cfg.locus_tag_start == 10 and cfg.locus_tag_step == 10
    assert cfg.transl_table == 1
    assert cfg.gap_min_length == 10 and cfg.gap_estimated_length == "known"
    assert cfg.product_default == "hypothetical protein"
    assert any(d.code == "unknown-config-key" and d.severity == Severity.WARNING for d in diags)


def test_missing_required_raises(tmp_path):
    p = tmp_path / "c.toml"
    p.write_text("[locus_tag]\nprefix='X'\n")  # no [source]
    with pytest.raises(GffParseError):
        load_config(str(p))


def test_load_common_verbatim(tmp_path):
    p = tmp_path / "common.tsv"
    p.write_text("COMMON\tDBLINK\t\tproject\tPRJDB1\n\t\t\tbiosample\tSAMD1\n")
    rows = load_common(str(p))
    assert rows == ["COMMON\tDBLINK\t\tproject\tPRJDB1", "\t\t\tbiosample\tSAMD1"]


def test_load_common_rejects_non_common(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("not common\n")
    with pytest.raises(GffParseError):
        load_common(str(p))


def test_missing_recommended_source_qualifiers_warn(tmp_path):
    p = tmp_path / "c.toml"
    p.write_text('[source]\nstrain = "X"\n[locus_tag]\nprefix = "P"\n')  # no organism/mol_type
    cfg, diags = load_config(str(p))
    codes = [d.code for d in diags]
    assert codes.count("source-missing-qualifier") == 2
