import json
import os
import shutil
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "gff_to_ddbj.sh"
ENV_BIN = "/lustre9/open/home/yt/micromamba/envs/mss_tools/bin"
FIX = Path("/lustre9/open/home/yt/ddbj/ddbj_mss_tools/tests/mss_fixtures")

COMMON_JSON = {
    "DBLINK": {"project": "PRJDB00000", "sample": "SAMD00000000"},
    "SOURCE": {"organism": "Testus fixtureus", "mol_type": "genomic DNA"},
    "SOURCE_IDENTIFIER": "strain",
}


def test_driver_produces_ann_and_fasta(tmp_path):
    # stage a known-good minimal gff2mss fixture set
    for f in ("mini.gff3", "mini.fa", "config.toml"):
        shutil.copy(FIX / f, tmp_path / f)
    common_path = tmp_path / "common.json"
    common_path.write_text(json.dumps(COMMON_JSON))
    out_prefix = tmp_path / "submission" / "mini"
    r = subprocess.run(
        ["bash", str(SCRIPT),
         "--gff", str(tmp_path / "mini.gff3"),
         "--fasta", str(tmp_path / "mini.fa"),
         "--mss-config", str(tmp_path / "config.toml"),
         "--common", str(common_path),
         "--out-prefix", str(out_prefix),
         "--workdir", str(tmp_path / "work")],
        env={**os.environ, "GFF_TO_DDBJ_ENV_BIN": ENV_BIN},
        capture_output=True, text=True)
    assert r.returncode == 0, f"driver failed:\nSTDOUT{r.stdout}\nSTDERR{r.stderr}"
    ann = out_prefix.with_suffix(".ann")
    fasta = out_prefix.with_suffix(".fasta")
    assert ann.exists() and ann.stat().st_size > 0
    assert fasta.exists() and fasta.stat().st_size > 0
