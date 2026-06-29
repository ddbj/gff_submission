from ddbj_gff.attributes import (
    decode_value,
    encode_value,
    parse_attributes,
    serialize_attributes,
)


def test_encode_decode_reserved_roundtrip():
    raw = "a;b=c,d&e%f\tg"
    enc = encode_value(raw)
    assert ";" not in enc and "=" not in enc and "," not in enc
    assert "%3B" in enc and "%3D" in enc and "%2C" in enc and "%25" in enc
    assert "%26" in enc  # & encoded
    assert "%09" in enc  # tab encoded
    assert decode_value(enc) == raw


def test_parse_simple_and_order_preserved():
    attrs = parse_attributes("ID=gene1;Name=psbA;locus_tag=AKK66")
    assert attrs == {"ID": ["gene1"], "Name": ["psbA"], "locus_tag": ["AKK66"]}
    assert list(attrs.keys()) == ["ID", "Name", "locus_tag"]


def test_parse_multivalue_dbxref():
    attrs = parse_attributes("Dbxref=GenBank:YP_1,GeneID:29")
    assert attrs["Dbxref"] == ["GenBank:YP_1", "GeneID:29"]


def test_parse_decodes_percent_and_keeps_literal_tilde():
    attrs = parse_attributes("Note=LSC%3B~large single-copy region")
    assert attrs["Note"] == ["LSC;~large single-copy region"]


def test_parse_empty_or_dot():
    assert parse_attributes("") == {}
    assert parse_attributes(".") == {}


def test_serialize_roundtrip():
    col9 = "ID=cds1;Note=has%3Bsemicolon;Dbxref=A:1,B:2"
    assert serialize_attributes(parse_attributes(col9)) == col9
