# src/ddbj_gff/repair/cds.py
from __future__ import annotations

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature

from .. import aa_names
from . import translate as _translate


def collect_transl_excepts(cds) -> list[str]:
    """transl_except specs from the CDS attribute and recoded_codon/stop_codon children."""
    specs = list(cds.attributes.get("transl_except", []))
    for child in cds.children:
        if child.type not in ("recoded_codon", "stop_codon"):
            continue
        sp = child.spans[0]
        loc = f"{sp.start}..{sp.end}" if sp.start != sp.end else f"{sp.start}"
        if sp.strand == "-":
            loc = f"complement({loc})"
        if child.type == "stop_codon":
            aa = "Term"
        else:
            aa = aa_names.to_abbrev((child.attributes.get("codon_redefined") or [""])[0])
        specs.append(f"(pos:{loc},aa:{aa})")
    return specs


def _table_id(cds, ctx) -> int:
    return cds.transl_table or ctx.transl_table


def coding_sequence(cds, ctx):
    """Return (coding_sequence_after_codon_start_upper, codon_table)."""
    seqid = cds.spans[0].seqid
    parent = ctx.sequences[seqid]
    coding = str(cds.to_biopython_location().extract(parent)).upper()
    codon_start = cds.codon_start or 1
    table = CodonTable.ambiguous_generic_by_id[int(_table_id(cds, ctx))]
    return coding[codon_start - 1:], table


def protein_of(cds, ctx) -> str:
    """Translate a CDS (transl_except-aware) to protein, trailing stop stripped."""
    table_id = _table_id(cds, ctx)
    excepts = collect_transl_excepts(cds)
    if excepts:
        sf = SeqFeature(cds.to_biopython_location(), type="CDS",
                        qualifiers={"transl_table": [str(table_id)],
                                    "codon_start": [str(cds.codon_start or 1)],
                                    "transl_except": excepts})
        return str(_translate.translate_cds_with_transl_except(sf, ctx.sequences[cds.spans[0].seqid]))
    coding, _ = coding_sequence(cds, ctx)
    protein = str(Seq(coding).translate(table=int(table_id)))
    return protein[:-1] if protein.endswith("*") else protein


def has_internal_stop(protein: str) -> bool:
    body = protein[:-1] if protein.endswith("*") else protein
    return "*" in body
