from Bio import SeqIO  # noqa: F401  (ensures biopython importable)
from ddbj_gff.mss.cli import main

GFF = (
    "##gff-version 3\n"
    "chr1\tS\tgene\t1\t9\t.\t+\t.\tID=g1\n"
    "chr1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g1.1;Parent=g1\n"
    "chr1\tS\texon\t1\t9\t.\t+\t.\tID=e1;Parent=g1.1\n"
    "chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g1.1\n"
)
FASTA = ">chr1\nATGAAATAA\n"
CONFIG = '[source]\norganism = "Foo bar"\n[locus_tag]\nprefix = "PFX"\n'
COMMON = "COMMON\tDBLINK\t\tproject\tPRJDB1\n"


def test_cli_writes_ann_and_fasta(tmp_path):
    (tmp_path / "g.gff").write_text(GFF)
    (tmp_path / "g.fa").write_text(FASTA)
    (tmp_path / "c.toml").write_text(CONFIG)
    (tmp_path / "common.tsv").write_text(COMMON)
    out = tmp_path / "result"
    rc = main(["--gff", str(tmp_path / "g.gff"), "--fasta", str(tmp_path / "g.fa"),
               "--config", str(tmp_path / "c.toml"), "--common", str(tmp_path / "common.tsv"),
               "--out", str(out)])
    assert rc == 0
    ann = (tmp_path / "result.ann").read_text()
    assert ann.startswith("COMMON\tDBLINK")
    assert "\tsource\t1..9\t" in ann
    assert "\tCDS\t1..9\tlocus_tag\tPFX_000010" in ann
    fasta = (tmp_path / "result.fasta").read_text()
    assert fasta.startswith(">chr1") and fasta.rstrip().endswith("//")


def test_cli_missing_sequence_sets_error_exit_and_stderr(tmp_path, capsys):
    # FASTA lacks chr1 -> convert emits a missing-sequence ERROR -> rc 1 + stderr summary
    (tmp_path / "g.gff").write_text(GFF)
    (tmp_path / "g.fa").write_text(">other\nACGTACGT\n")
    (tmp_path / "c.toml").write_text(CONFIG)
    (tmp_path / "common.tsv").write_text(COMMON)
    out = tmp_path / "result"
    rc = main(["--gff", str(tmp_path / "g.gff"), "--fasta", str(tmp_path / "g.fa"),
               "--config", str(tmp_path / "c.toml"), "--common", str(tmp_path / "common.tsv"),
               "--out", str(out)])
    assert rc == 1
    assert "ERROR" in capsys.readouterr().err
