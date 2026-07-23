import re
from pathlib import Path

import yaml  # PyYAML, present in the mss_tools env

SKILL = Path(__file__).resolve().parents[1] / "SKILL.md"


def _frontmatter_and_body(text):
    m = re.match(r"^---\n(.*?)\n---\n(.*)$", text, re.DOTALL)
    assert m, "SKILL.md must start with a YAML frontmatter block delimited by ---"
    return yaml.safe_load(m.group(1)), m.group(2)


def test_frontmatter_required_and_hermes_config():
    fm, _ = _frontmatter_and_body(SKILL.read_text())
    assert fm["name"] == "gff-to-ddbj"
    assert isinstance(fm["description"], str) and fm["description"].strip()
    assert fm["platforms"] == ["linux"]
    cfg = {c["key"]: c for c in fm["metadata"]["hermes"]["config"]}
    assert set(cfg) == {"gff_to_ddbj.env_bin", "gff_to_ddbj.validator_dir"}
    assert cfg["gff_to_ddbj.env_bin"]["default"] == \
        "/lustre9/open/home/yt/micromamba/envs/mss_tools/bin"
    assert cfg["gff_to_ddbj.validator_dir"]["default"] == \
        "/home/w3const/ddbj-validator-production"
    for c in cfg.values():
        assert {"key", "description", "default", "prompt"} <= set(c)
    assert "terminal" in fm["metadata"]["hermes"]["requires_toolsets"]


def test_description_has_no_workflow_summary():
    fm, _ = _frontmatter_and_body(SKILL.read_text())
    # description is triggering-conditions only; must start with "Use when"
    assert fm["description"].lstrip().startswith("Use when")


def test_body_has_five_sections_in_order():
    _, body = _frontmatter_and_body(SKILL.read_text())
    heads = re.findall(r"^## (.+)$", body, re.MULTILINE)
    assert heads[:5] == ["When to use", "Quick reference", "Procedure",
                         "Pitfalls", "Verification"]
