from __future__ import annotations

import re

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import (AfterPosition, BeforePosition, CompoundLocation,
                            FeatureLocation)
from Bio.SeqIO.InsdcIO import _insdc_location_string

from Bio.SeqFeature import SeqFeature

from ..errors import Diagnostic, GffParseError, Severity
from .. import aa_names
from .config import MssConfig
from .gaps import assembly_gap_features
from .locus_tag import LocusTagAssigner
from .model import MssDocument, MssEntry, MssFeature, MssQualifier
from .translate import translate_cds_with_transl_except

_STRAND = {"+": 1, "-": -1}


def _collect_transl_excepts(cds_feat) -> list:
    """Gather transl_except specs from the CDS attribute (raw) and recoded_codon/stop_codon children."""
    specs = list(cds_feat.attributes.get("transl_except", []))
    for child in cds_feat.children:
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


def collect_spans(parent, feature_type: str) -> list:
    spans = []
    for child in parent.children:
        if child.type == feature_type:
            spans.extend(child.spans)
    return spans


def _ordered(spans: list):
    strand = spans[0].strand
    return sorted(spans, key=lambda s: s.start, reverse=(strand == "-"))


def build_insdc_location(spans: list, seqlen: int,
                         five_prime_partial: bool = False,
                         three_prime_partial: bool = False) -> str:
    bio = _STRAND.get(spans[0].strand, 0)
    ordered = _ordered(spans)
    locs = [FeatureLocation(s.start - 1, s.end, strand=bio) for s in ordered]
    if five_prime_partial:
        f = locs[0]
        if bio == -1:
            locs[0] = FeatureLocation(f.start, AfterPosition(int(f.end)), strand=bio)
        else:
            locs[0] = FeatureLocation(BeforePosition(int(f.start)), f.end, strand=bio)
    if three_prime_partial:
        l = locs[-1]
        if bio == -1:
            locs[-1] = FeatureLocation(BeforePosition(int(l.start)), l.end, strand=bio)
        else:
            locs[-1] = FeatureLocation(l.start, AfterPosition(int(l.end)), strand=bio)
    compound = locs[0] if len(locs) == 1 else CompoundLocation(locs)
    return _insdc_location_string(compound, seqlen)


def extract_seq(spans: list, genome_seq):
    bio = _STRAND.get(spans[0].strand, 0)
    ordered = _ordered(spans)
    locs = [FeatureLocation(s.start - 1, s.end, strand=bio) for s in ordered]
    compound = locs[0] if len(locs) == 1 else CompoundLocation(locs)
    return compound.extract(genome_seq)


def build_source_feature(seqid: str, seqlen: int, cfg: MssConfig) -> MssFeature:
    quals = [MssQualifier(k, v) for k, v in cfg.source.items()]
    chromosome = None
    if cfg.chromosome_pattern:
        m = re.match(cfg.chromosome_pattern, seqid)
        if m:
            chromosome = m.group(1)
    if chromosome is not None:
        quals.append(MssQualifier("chromosome", chromosome))
        if cfg.ff_def_chromosome:
            quals.append(MssQualifier("ff_definition", cfg.ff_def_chromosome))
    else:
        quals.append(MssQualifier("submitter_seqid", seqid))
        if cfg.ff_def_default:
            quals.append(MssQualifier("ff_definition", cfg.ff_def_default))
    return MssFeature("source", f"1..{seqlen}", quals)


def _submitter_note(gene, mrna) -> MssQualifier:
    return MssQualifier(
        "note", f"submitter_gene_id: {gene.id}, submitter_transcript_id: {mrna.id}")


def mrna_partial_flags(mrna) -> tuple[bool, bool]:
    exon = collect_spans(mrna, "exon")
    cds = collect_spans(mrna, "CDS")
    if not exon or not cds:
        return (False, False)
    exon_lo, exon_hi = min(s.start for s in exon), max(s.end for s in exon)
    cds_lo, cds_hi = min(s.start for s in cds), max(s.end for s in cds)
    strand = exon[0].strand
    # no UTR on a genomic side -> that side is partial; map to 5'/3' by strand
    left_partial = exon_lo == cds_lo
    right_partial = exon_hi == cds_hi
    if strand == "-":
        return (right_partial, left_partial)   # 5' = genomic right on minus strand
    return (left_partial, right_partial)


def build_mrna_feature(mrna, gene, locus_tag: str, seqlen: int) -> MssFeature:
    spans = collect_spans(mrna, "exon") or collect_spans(mrna, "CDS")
    fp, tp = mrna_partial_flags(mrna)
    location = build_insdc_location(spans, seqlen, fp, tp)
    quals = [MssQualifier("locus_tag", locus_tag)]
    if gene.gene:
        quals.append(MssQualifier("gene", gene.gene))
    quals.append(_submitter_note(gene, mrna))
    return MssFeature("mRNA", location, quals)


def _product(mrna, gene, cfg: MssConfig) -> str:
    pmap = cfg.product_map or {}
    hit = pmap.get(mrna.id) or (pmap.get(gene.id) if gene and gene.id else None)
    if hit:
        return hit
    vals = mrna.attributes.get("product")
    if vals and vals[0]:
        return vals[0]
    # col-9 product may live on the CDS child (e.g. organelle LiftOff annotations)
    for c in mrna.children:
        if c.type == "CDS":
            cvals = c.attributes.get("product")
            if cvals and cvals[0]:
                return cvals[0]
    return cfg.product_default


def build_cds_feature(mrna, gene, locus_tag: str, genome_seq, cfg: MssConfig,
                      diagnostics: list) -> "MssFeature | None":
    spans = collect_spans(mrna, "CDS")
    if not spans:
        diagnostics.append(Diagnostic(Severity.WARNING, None, "no-cds",
                                      f"mRNA {mrna.id!r} has no CDS; skipped"))
        return None
    ordered = _ordered(spans)
    phase = ordered[0].phase
    codon_start = 1 if phase is None else phase + 1

    # transl_table: CDS attribute wins, else config default
    table_id = cfg.transl_table
    for child in mrna.children:
        if child.type == "CDS" and child.transl_table is not None:
            table_id = child.transl_table
            break
    table = CodonTable.unambiguous_dna_by_id[table_id]

    cds_seq = extract_seq(spans, genome_seq)
    coding = str(cds_seq[codon_start - 1:]).upper()
    coding_full = coding[: len(coding) - len(coding) % 3]
    first_codon = coding[:3]
    # last codon must be evaluated in-frame (coding_full), not the trailing 3 bases:
    # a CDS whose length is not a multiple of 3 otherwise misjudges 3' completeness.
    last_codon = coding_full[-3:] if len(coding_full) >= 3 else ""
    five_prime_partial = codon_start != 1 or first_codon not in table.start_codons
    three_prime_partial = last_codon not in table.stop_codons

    # translation validation (report only)
    if len(coding) % 3 != 0:
        diagnostics.append(Diagnostic(Severity.WARNING, None, "translation-not-multiple-of-3",
                                      f"CDS {mrna.id!r} coding length not a multiple of 3"))
    cds_feat = next((c for c in mrna.children if c.type == "CDS"), None)
    excepts = _collect_transl_excepts(cds_feat) if cds_feat is not None else []
    if excepts:
        sf = SeqFeature(cds_feat.to_biopython_location(), type="CDS",
                        qualifiers={"transl_table": [str(table_id)],
                                    "codon_start": [str(codon_start)],
                                    "transl_except": excepts})
        protein = str(translate_cds_with_transl_except(sf, genome_seq))
    else:
        protein = str(Seq(coding_full).translate(table=table_id))
    body = protein[:-1] if protein.endswith("*") else protein
    if "*" in body:
        diagnostics.append(Diagnostic(Severity.WARNING, None, "translation-internal-stop",
                                      f"CDS {mrna.id!r} has an internal stop codon"))
        loc = build_insdc_location(spans, len(genome_seq))
        note = f"internal stop codon(s) detected in CDS {mrna.id}; not translated"
        quals = [MssQualifier("locus_tag", locus_tag), MssQualifier("note", note)]
        return MssFeature("misc_feature", loc, quals)
    if not five_prime_partial and not protein.startswith("M"):
        diagnostics.append(Diagnostic(Severity.WARNING, None, "translation-no-start",
                                      f"CDS {mrna.id!r} does not start with M"))

    seqlen = len(genome_seq)
    location = build_insdc_location(spans, seqlen, five_prime_partial, three_prime_partial)
    quals = [
        MssQualifier("locus_tag", locus_tag),
        MssQualifier("transl_table", str(table_id)),
        MssQualifier("codon_start", str(codon_start)),
        MssQualifier("product", _product(mrna, gene, cfg)),
    ]
    if gene.gene or mrna.gene:
        quals.append(MssQualifier("gene", gene.gene or mrna.gene))
    inference = mrna.attributes.get("inference")
    if inference:
        quals.append(MssQualifier("inference", inference[0]))
    quals.append(_submitter_note(gene, mrna))
    return MssFeature("CDS", location, quals)


def _representative_mrna(gene, diagnostics):
    mrnas = [c for c in gene.children if c.type == "mRNA"]
    if not mrnas:
        diagnostics.append(Diagnostic(Severity.WARNING, None, "no-transcript",
                                      f"gene {gene.id!r} has no mRNA; skipped"))
        return None
    if len(mrnas) > 1:
        diagnostics.append(Diagnostic(Severity.WARNING, None, "multi-transcript",
                                      f"gene {gene.id!r} has {len(mrnas)} transcripts; keeping one representative transcript"))
    for m in mrnas:
        if m.id and m.id.endswith(".1"):
            return m
    return mrnas[0]


def _span_start(feature) -> int:
    return min(s.start for s in feature.spans) if feature.spans else 0


_RNA_MAP = {
    "pre_miRNA": "precursor_RNA",
    "miRNA": "ncRNA",
    "ncRNA": "ncRNA",
    "snRNA": "ncRNA",
    "snoRNA": "ncRNA",
    "tRNA": "tRNA",
    "rRNA": "rRNA",
    "tmRNA": "tmRNA",
}

_STRUCTURAL = {"exon", "CDS", "intron", "five_prime_UTR", "three_prime_UTR",
               "start_codon", "stop_codon"}

_PARENTLESS_RNA_TYPES = set(_RNA_MAP) | {"ncRNA", "tRNA", "rRNA"}
_NCRNA_KNOWN = {"snRNA", "snoRNA", "miRNA", "siRNA", "scRNA", "antisense_RNA",
                "ribozyme", "RNase_P_RNA", "telomerase_RNA", "lncRNA", "SRP_RNA",
                "guide_RNA", "vault_RNA", "Y_RNA", "autocatalytically_spliced_intron"}


def _submitter_note_ids(gene_id, tx_id) -> MssQualifier:
    return MssQualifier("note", f"submitter_gene_id: {gene_id}, submitter_transcript_id: {tx_id}")


def build_rna_feature(rna, locus_tag: str, seqlen: int, gene_id: str, tx_id: str,
                      gene_name: str | None = None) -> MssFeature:
    feat_key = _RNA_MAP.get(rna.type, "misc_RNA")
    spans = collect_spans(rna, "exon") or rna.spans
    location = build_insdc_location(spans, seqlen)
    quals = [MssQualifier("locus_tag", locus_tag)]
    if feat_key == "ncRNA":
        klass = rna.type if rna.type in _NCRNA_KNOWN else "other"
        quals.append(MssQualifier("ncRNA_class", klass))
    product = rna.product
    if not product and rna.type == "tRNA":
        iso = rna._first("isotype")
        if iso:
            product = f"tRNA-{iso}"
    if product:
        quals.append(MssQualifier("product", product))
    if gene_name:
        quals.append(MssQualifier("gene", gene_name))
    for x in rna.dbxref:
        quals.append(MssQualifier("db_xref", x))
    for note_val in rna.note:
        quals.append(MssQualifier("note", note_val))
    for note_val in rna.attributes.get("note", []):
        quals.append(MssQualifier("note", note_val))
    quals.append(_submitter_note_ids(gene_id, tx_id))
    return MssFeature(feat_key, location, quals)


def build_noncoding_features(gene, locus_tag: str, seqlen: int, cfg) -> list:
    features = []
    for rna in gene.children:
        if rna.type in _STRUCTURAL:
            continue
        features.append(build_rna_feature(rna, locus_tag, seqlen, gene.id, rna.id,
                                           gene.gene or rna.gene))
    return features


def _set_submitter_transcripts(cds, gene, transcript_ids) -> None:
    value = (f"submitter_gene_id: {gene.id}, "
             f"submitter_transcript_id: {', '.join(transcript_ids)}")
    for q in cds.qualifiers:
        if q.key == "note" and q.value.startswith("submitter_gene_id:"):
            q.value = value
            return
    cds.qualifiers.append(MssQualifier("note", value))


def build_gene_features(gene, mode, assigner, genome_seq, cfg, diagnostics) -> list:
    transcripts = [c for c in gene.children if c.type == "mRNA"]
    if not transcripts:
        rna_children = [c for c in gene.children if c.type not in _STRUCTURAL]
        if not rna_children:
            diagnostics.append(Diagnostic(Severity.WARNING, None, "no-rna",
                                          f"gene {gene.id!r} has no mRNA or RNA child; skipped"))
            return []
        return build_noncoding_features(gene, assigner.assign(gene), len(genome_seq), cfg)

    transcripts = sorted(transcripts, key=lambda m: m.id or "")
    if mode == "minimal":
        rep = _representative_mrna(gene, diagnostics)
        transcripts = [rep] if rep is not None else []

    features = []
    misc_feats: list = []
    locus_tag = None
    cds_index = {}
    cds_order = []
    for mrna in transcripts:
        if not collect_spans(mrna, "exon") and not collect_spans(mrna, "CDS"):
            diagnostics.append(Diagnostic(Severity.WARNING, None, "no-exon",
                                          f"mRNA {mrna.id!r} has no exon or CDS; skipped"))
            continue
        if locus_tag is None:
            locus_tag = assigner.assign(gene)
        cds = build_cds_feature(mrna, gene, locus_tag, genome_seq, cfg, diagnostics)
        if cds is not None and cds.key == "misc_feature":
            misc_feats.append(cds)
            continue
        features.append(build_mrna_feature(mrna, gene, locus_tag, len(genome_seq)))
        if cds is None:
            continue
        if mode == "nonredundant":
            if cds.location in cds_index:
                cds_index[cds.location][1].append(mrna.id)
            else:
                cds_index[cds.location] = [cds, [mrna.id]]
                cds_order.append(cds.location)
        else:
            features.append(cds)
    if mode == "nonredundant":
        for loc in cds_order:
            cds, tids = cds_index[loc]
            if len(tids) > 1:
                _set_submitter_transcripts(cds, gene, tids)
            features.append(cds)
    if not features and misc_feats:
        features.append(misc_feats[0])
    return features


def _seqids_in_order(doc) -> list:
    seen: list = []
    for feat in doc.features:
        for s in feat.spans:
            if s.seqid not in seen:
                seen.append(s.seqid)
    return seen


def build_entry_features(doc, seqs, cfg, diagnostics: list) -> dict:
    """Per-seqid feature blocks (gene/RNA/misc). No source, no assembly_gap."""
    assigner = LocusTagAssigner.from_config(cfg)
    result: dict = {}
    for seqid in _seqids_in_order(doc):
        if seqid not in seqs:
            diagnostics.append(Diagnostic(Severity.ERROR, None, "missing-sequence",
                                          f"seqid {seqid!r} not found in FASTA; entry skipped"))
            continue
        genome_seq = seqs[seqid]
        genes = [f for f in doc.roots if f.type == "gene"
                 and any(s.seqid == seqid for s in f.spans)]
        parentless = [f for f in doc.roots
                      if f.type in _PARENTLESS_RNA_TYPES
                      and any(s.seqid == seqid for s in f.spans)]

        def _is_pseudogene(f):
            return f.type == "pseudogene" or f._first("gene_biotype") == "pseudogene"

        pseudo_roots = [f for f in doc.roots
                        if _is_pseudogene(f) and any(s.seqid == seqid for s in f.spans)]
        skipped = len(pseudo_roots)
        genes = [g for g in genes if not _is_pseudogene(g)]
        parentless = [r for r in parentless if not _is_pseudogene(r)]
        if skipped:
            diagnostics.append(Diagnostic(Severity.WARNING, None, "pseudogene-skipped",
                                          f"{seqid}: skipped {skipped} pseudogene feature(s)"))
        items = [(_span_start(g), g) for g in genes] + [(_span_start(r), r) for r in parentless]
        items.sort(key=lambda t: t[0])
        feats: list = []
        for _, feat in items:
            if feat.type == "gene":
                feats.extend(build_gene_features(feat, cfg.transcript_mode, assigner,
                                                 genome_seq, cfg, diagnostics))
            else:
                feats.append(build_rna_feature(feat, assigner.assign(feat),
                                                len(genome_seq), feat.id, feat.id, feat.gene))
        result[seqid] = feats
    return result


def convert(doc, seqs, cfg, common_rows, *, strict: bool = False):
    diagnostics: list = []
    per_entry = build_entry_features(doc, seqs, cfg, diagnostics)
    entries: list = []
    for seqid, feats in per_entry.items():
        genome_seq = seqs[seqid]
        entry_feats = [build_source_feature(seqid, len(genome_seq), cfg)]
        entry_feats.extend(assembly_gap_features(str(genome_seq), cfg))
        entry_feats.extend(feats)
        entries.append(MssEntry(seqid, entry_feats))
    if strict:
        for d in diagnostics:
            if d.severity == Severity.ERROR:
                raise GffParseError(d)
    return MssDocument(common_rows, entries), diagnostics
