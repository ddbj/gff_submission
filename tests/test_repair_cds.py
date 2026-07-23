from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from ddbj_gff.repair.context import RepairContext
from ddbj_gff.repair.cds import (collect_transl_excepts, protein_of,
                                  has_internal_stop, coding_sequence)


def _cds(seqid, start, end, strand="+", phase=0, **attrs):
    a = {"ID": ["c1"]}
    a.update({k: (v if isinstance(v, list) else [v]) for k, v in attrs.items()})
    return Feature("c1", "src", "CDS", [Span(seqid, start, end, strand, phase)], a)


def test_has_internal_stop():
    assert has_internal_stop("MKV*") is False        # trailing stop only
    assert has_internal_stop("MK*V") is True          # internal stop
    assert has_internal_stop("MKV") is False


def test_protein_of_clean_cds():
    # ATG AAA GTT TAA -> M K V (stop stripped by translate cds path)
    seq = Seq("ATGAAAGTTTAA")
    ctx = RepairContext(sequences={"s": seq}, transl_table=1)
    cds = _cds("s", 1, 12, "+", 0, transl_table="1")
    prot = protein_of(cds, ctx)
    assert prot.startswith("MKV")
    assert has_internal_stop(prot) is False


def test_protein_of_internal_stop():
    # ATG TAA AAA TAA -> M * K (internal stop at aa index 1)
    seq = Seq("ATGTAAAAATAA")
    ctx = RepairContext(sequences={"s": seq}, transl_table=1)
    cds = _cds("s", 1, 12, "+", 0, transl_table="1")
    prot = protein_of(cds, ctx)
    assert has_internal_stop(prot) is True


def test_coding_sequence_respects_codon_start():
    seq = Seq("GATGAAAGTTTAA")   # leading G, codon_start=2
    ctx = RepairContext(sequences={"s": seq})
    cds = _cds("s", 1, 13, "+", 1)   # phase 1 -> codon_start 2
    coding, table = coding_sequence(cds, ctx)
    assert coding.startswith("ATGAAAGTT")
