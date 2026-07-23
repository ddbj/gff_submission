# Hermes skills (in-repo source)

Source for Hermes agent skills that expose this repo's tools. Skills are authored per
the Hermes format (agentskills.io + `metadata.hermes`).

## gff-to-ddbj

End-to-end GFF3 → DDBJ MSS workflow (normalize → validate → repair → gff2mss →
ddbj-validator).

### Install

Copy or symlink the skill into the Hermes skills tree (category `bioinformatics`):

```bash
mkdir -p ~/.hermes/skills/bioinformatics
ln -s "$(pwd)/hermes/skills/gff-to-ddbj" ~/.hermes/skills/bioinformatics/gff-to-ddbj
```

### Config keys (Hermes `config.yaml` → `skills.config`)

| key | default (NIG/DDBJ cluster) | purpose |
|---|---|---|
| `gff_to_ddbj.env_bin` | `/lustre9/open/home/yt/micromamba/envs/mss_tools/bin` | env with `python` (ddbj_gff) + `gff2mss` |
| `gff_to_ddbj.validator_dir` | `/home/w3const/ddbj-validator-production` | `ddbj-validator` install |

Override per install if the tools live elsewhere.
