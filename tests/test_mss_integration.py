from pathlib import Path
import pytest
from Bio import SeqIO
from ddbj_gff import parse
from ddbj_gff.errors import Severity
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import convert

ROOT = Path(__file__).resolve().parents[1]
GFF = ROOT / "examples" / "marchantia" / "MpTak1_v7.1.marpolbase.gff"
FASTA = ROOT / "examples" / "marchantia" / "Tak-1_hifi_v2.nextpolish.fa"

pytestmark = pytest.mark.slow


def _require(p: Path):
    if not p.exists():
        pytest.skip(f"example file not present: {p}")


def test_marchantia_converts_without_errors():
    _require(GFF)
    _require(FASTA)
    doc = parse(GFF.read_text(errors="replace"))
    seqs = {rec.id: rec.seq for rec in SeqIO.parse(str(FASTA), "fasta")}
    cfg = MssConfig(source={"organism": "Marchantia polymorpha subsp. ruderalis",
                            "mol_type": "genomic DNA"},
                    chromosome_pattern="^chr(.+)$",
                    ff_def_chromosome="@@[organism]@@ DNA, chromosome: @@[chromosome]@@",
                    ff_def_default="@@[organism]@@ DNA, @@[entry]@@",
                    locus_tag_prefix="MPTK1")
    mss, diags = convert(doc, seqs, cfg, ["COMMON\tDBLINK"])
    errors = [d for d in diags if d.severity == Severity.ERROR]
    assert errors == [], f"unexpected ERROR diagnostics: {errors[:5]}"
    # at least one entry with mRNA+CDS produced
    total_cds = sum(1 for e in mss.entries for f in e.features if f.key == "CDS")
    assert total_cds > 1000
    # every CDS feature carries locus_tag, codon_start, product, transl_table
    for e in mss.entries:
        for f in e.features:
            if f.key == "CDS":
                keys = {q.key for q in f.qualifiers}
                assert {"locus_tag", "codon_start", "product", "transl_table"} <= keys
