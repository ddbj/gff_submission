from ddbj_gff import parse
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.normalize.passes import NormalizeContext, pass_directives

GFF = (
    "##gff-version 3\n"
    "chr1\tS\tgene\t100\t900\t.\t+\t.\tID=g;locus_tag=ABC_1\n"
    "chr1\tS\tCDS\t130\t870\t.\t+\t0\tID=c;Parent=g\n"
)  # missing insdc-gff-version, species, sequence-region, transl_table


def _ctx(seq_lengths=None, **cfg):
    return NormalizeContext(vocab=None, seq_lengths=seq_lengths, config=NormalizeConfig(**cfg))


def test_adds_insdc_version_and_transl_table():
    doc = parse(GFF)
    pass_directives(doc, _ctx(taxid=3702))
    assert doc.insdc_gff_version == "1.0.0"
    assert doc.transl_table_map == {"primary": 1}


def test_adds_species_from_taxid():
    doc = parse(GFF)
    pass_directives(doc, _ctx(taxid=3702))
    assert doc.species == 3702


def test_no_taxid_reports_unresolved_and_skips_species():
    doc = parse(GFF)
    changes = pass_directives(doc, _ctx())  # no taxid
    assert doc.species is None
    assert any(c.action == "no-taxid" for c in changes)


def test_sequence_region_from_seq_lengths():
    doc = parse(GFF)
    pass_directives(doc, _ctx(seq_lengths={"chr1": 10000}, taxid=3702))
    assert doc.sequence_regions["chr1"] == (1, 10000)


def test_sequence_region_approx_when_no_fasta():
    doc = parse(GFF)
    changes = pass_directives(doc, _ctx(taxid=3702))
    assert doc.sequence_regions["chr1"] == (1, 900)  # max feature end
    assert any(c.action == "approx-region" for c in changes)


def test_transl_table_promotes_consistent_cds_value():
    g = (
        "##gff-version 3\n"
        "chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=c;transl_table=11\n"
    )
    doc = parse(g)
    pass_directives(doc, _ctx(taxid=3702))
    assert doc.transl_table_map == {"primary": 11}  # promoted, not default 1


def test_idempotent():
    doc = parse(GFF)
    pass_directives(doc, _ctx(seq_lengths={"chr1": 10000}, taxid=3702))
    n1 = len(doc.directives)
    pass_directives(doc, _ctx(seq_lengths={"chr1": 10000}, taxid=3702))
    assert len(doc.directives) == n1  # no duplicate directives
