import gzip
from pathlib import Path

import pytest

from ddbj_gff import parse, write

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"

pytestmark = pytest.mark.slow


def _read_text(path: Path) -> str:
    if path.suffix == ".gz":
        return gzip.decompress(path.read_bytes()).decode("ascii", errors="replace")
    return path.read_text(errors="replace")


def _require(path: Path) -> str:
    if not path.exists():
        pytest.skip(f"example file not present: {path}")
    return _read_text(path)


def test_rice_cp_roundtrip_and_trans_splicing():
    text = _require(EXAMPLES / "rice_cp" / "rice_cp.gff3")
    doc = parse(text)
    # rps12 trans-spliced CDS has 3 segments sharing one ID
    cds = doc.get("cds-YP_009305283.1")
    assert cds is not None and len(cds.spans) == 3
    assert parse(write(doc)).semantically_equals(doc)


def test_chloroplast_ddbj_roundtrip():
    text = _require(EXAMPLES / "marchantia" / "chloroplast.gff3")
    doc = parse(text)
    assert parse(write(doc)).semantically_equals(doc)


def test_ecoli_transl_except_present_and_roundtrip():
    text = _require(EXAMPLES / "ecoli" / "GCF_000005845.2_ASM584v2_genomic.gff.gz")
    doc = parse(text)
    has_transl_except = any("transl_except" in f.attributes for f in doc.features)
    assert has_transl_except
    assert parse(write(doc)).semantically_equals(doc)


def test_arabidopsis_parses_within_budget():
    import time

    text = _require(EXAMPLES / "arabidopsis" / "AT_chr1.gff3")
    t0 = time.perf_counter()
    doc = parse(text)
    elapsed = time.perf_counter() - t0
    assert len(doc.features) > 1000
    assert elapsed < 120  # loose ceiling: 121MB must parse in-memory in reasonable time
