from ddbj_gff.parser import parse_directive


def test_gff_version():
    d = parse_directive("##gff-version 3")
    assert d.kind == "gff-version" and d.value == "3"


def test_sequence_region():
    d = parse_directive("##sequence-region NC_031333.1 1 134502")
    assert d.kind == "sequence-region"
    assert d.value == ("NC_031333.1", 1, 134502)


def test_species_extracts_taxid():
    d = parse_directive("##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4530")
    assert d.kind == "species" and d.value == 4530


def test_insdc_version_and_spec_version():
    assert parse_directive("#!insdc-gff-version 1.2.3").value == "1.2.3"
    assert parse_directive("#!gff-spec-version 1.21").kind == "gff-spec-version"


def test_transl_table_map():
    d = parse_directive("#!transl_table primary:1,chloroplast:11")
    assert d.kind == "transl_table"
    assert d.value == {"primary": 1, "chloroplast": 11}


def test_fasta_and_boundary_and_unknown():
    assert parse_directive("##FASTA").kind == "FASTA"
    assert parse_directive("###").kind == "resolution-boundary"
    assert parse_directive("##something-else foo").kind == "unknown"
