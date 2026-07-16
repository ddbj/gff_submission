from ddbj_gff.normalize.config import NormalizeConfig, load_normalize_config


def test_reparent_config_default_on():
    assert NormalizeConfig().reparent_gene_children is True


def test_reparent_config_loader_reads_flag(tmp_path):
    p = tmp_path / "n.toml"
    p.write_text("[normalize]\nreparent_gene_children = false\n")
    assert load_normalize_config(str(p)).reparent_gene_children is False


def test_reparent_config_loader_defaults_on(tmp_path):
    p = tmp_path / "n.toml"
    p.write_text("[normalize]\ntaxid = 1\n")
    assert load_normalize_config(str(p)).reparent_gene_children is True
