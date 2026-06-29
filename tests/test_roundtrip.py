from pathlib import Path

import pytest

from ddbj_gff import parse, write

FIXTURES = Path(__file__).parent / "fixtures"
FILES = [
    "canonical_gene.gff3",
    "discontinuous_cds.gff3",
    "trans_splicing.gff3",
    "transl_except.gff3",
    "circular.gff3",
    "attributes_escaping.gff3",
    "peptide_fasta.gff3",
]


@pytest.mark.parametrize("name", FILES)
def test_semantic_roundtrip(name):
    text = (FIXTURES / name).read_text()
    once = parse(text)
    twice = parse(write(once))
    assert twice.semantically_equals(once)


def test_trans_splicing_structure_preserved():
    text = (FIXTURES / "trans_splicing.gff3").read_text()
    doc = parse(text)
    gene = doc.get("gene-Mp_Cg00010")
    assert len(gene.spans) == 2
    assert gene.is_trans_spliced and gene.is_ordered
    intron = doc.get("id-Mp_Cg00010")
    assert intron.number == 1  # number is the qualifier, not the part order
    assert [s.part for s in intron.ordered_spans()] == [1, 2]


def test_discontinuous_cds_has_three_spans():
    text = (FIXTURES / "discontinuous_cds.gff3").read_text()
    assert len(parse(text).get("cds-ycf3").spans) == 3
