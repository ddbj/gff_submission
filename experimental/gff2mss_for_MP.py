#!/usr/bin/env python

import re
import sys
from Bio import SeqIO
from BCBio import GFF
from Bio.SeqIO.InsdcIO import _insdc_location_string
from Bio.SeqFeature import FeatureLocation, SeqFeature, BeforePosition, AfterPosition, ExactPosition, CompoundLocation
import Bio.Data.CodonTable

if len(sys.argv) > 1 and sys.argv[1] == "tak2":
    common_metadata_file = "common.metadata_tak2.tsv"
    input_gff_file = "MpTak2_v7.1.ddbj.gff"
    input_fasta_file = "MpTak2_v7.1.ddbj.fa"

    output_ann_file = "SAMD00647144_marchantia.ann"
    output_seq_file = "SAMD00647144_marchantia.fasta"

    # Qualifiers for source features.
    # The values described here is commonly applied to each entry.
    souce_qualifiers = {
        "organism": "Marchantia polymorpha subsp. ruderalis",
        "sub_species": "ruderalis",
        "strain": "Tak-2",
        "mol_type": "genomic DNA",
        "sex": "female",
        "collection_date": "missing: lab stock",
        "country": "missing: lab stock"

    }
    LOCUS_TAG_PREFIX = "MPTK2"
else:
    common_metadata_file = "common.metadata.tsv"
    input_gff_file = "MpTak1_v7.1.ddbj.gff"
    input_fasta_file = "MpTak1_v7.1.ddbj.fa"

    output_ann_file = "SAMD00647143_marchantia.ann"
    output_seq_file = "SAMD00647143_marchantia.fasta"

    # Qualifiers for source features.
    # The values described here is commonly applied to each entry.
    souce_qualifiers = {
        "organism": "Marchantia polymorpha subsp. ruderalis",
        "sub_species": "ruderalis",
        "strain": "Tak-1",
        "mol_type": "genomic DNA",
        "sex": "male",
        "collection_date": "missing: lab stock",
        "country": "missing: lab stock"

    }
    LOCUS_TAG_PREFIX = "MPTK1"


print("input_file (metadata):", common_metadata_file)
print("input_file (gff):", input_gff_file)
print("input_file (fasta):", input_fasta_file)

print("output_file (ann):", output_ann_file)
print("output_file (fasta):", output_seq_file)
print("locus tag prefix:", LOCUS_TAG_PREFIX)



# constant values for assembly_gap feature
# see https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
MIN_GAP_LENGTH = 10  # Recommendation from the curation staff
GAP_TYPE = "within scaffold"  # see https://www.ddbj.nig.ac.jp/ddbj/qualifiers.html#gap_type
LINKAGE_EVIDENCE = "align genus"  # paired-ends, align genus, proximity ligation, etc. 
gap_pattern = re.compile('n{%d,}' % MIN_GAP_LENGTH)


TRANSL_TABLE = "1"


# codon table
codon_table = Bio.Data.CodonTable.standard_dna_table
start_codons = codon_table.start_codons  # ['TTG', 'CTG', 'ATG']
stop_codons = codon_table.stop_codons  # ['TAA', 'TAG', 'TGA']


# TODO
# featureのsortは頻繁に行うので関数化する

def read_fasta_and_gff(input_fasta_file, input_gff_file):
    seq_dict = SeqIO.to_dict(SeqIO.parse(input_fasta_file, "fasta"))
    return list(GFF.parse(open(input_gff_file), base_dict=seq_dict))

def get_locus_tag(gene_feature):
    # This function is specific to Marchantia. 
    # use ID as locus_tag
    # Mp1g00010 --> MpTak1_1g00010
    # Mp1g00010_R1 --> MpTak1_1g00010R1
    locus_tag = gene_feature.id.replace("_", "").replace("Mp", LOCUS_TAG_PREFIX + "_")
    return locus_tag

def add_note_describing_submitter_id(gene_id, transcript_id):
    return ["", "", "", "note", f"submitter_gene_id: {gene_id}, submitter_transcript_id: {transcript_id}"]

def create_source_feature(record, souce_qualifiers):
    ret = []
    for key, value in souce_qualifiers.items():
        ret.append(["", "", "", key, value])
    if record.id.startswith("chr"):
        # this must be adjusted for the sequence identifiers in FASTA
        chromosome_name = record.id.replace("chr", "")
        ret.append(["", "", "", "chromosome", chromosome_name])
        ret.append(["", "", "", "ff_definition", "@@[organism]@@ @@[strain]@@ DNA, chromosome: @@[chromosome]@@"])    
    else:
        ret.append(["", "", "", "submitter_seqid", record.id])
        ret.append(["", "", "", "ff_definition", "@@[organism]@@ @@[strain]@@ DNA, @@[entry]@@"])

    ret[0][0] = record.id
    ret[0][1] = "source"
    ret[0][2] = "1..{}".format(len(record))
    return ret


def create_gap_feature(seq, length_estimated=False):
    """
    Create assembly_feature
    If length_estimated is True, estimated_length is "known", otherwise "unknown"
    MIN_GAP_LENGTH, GAP_TYPE, LINKAGE_EVIDENCE are specific to Marchantia genome, and need adjustment for other species.
    """
    ret = []
    length_value = "known" if length_estimated else "unknown"

    for match in gap_pattern.finditer(seq.lower()):
        start = match.start() + 1  # +1 because INSDC coordinate is 1-based
        end = match.end()
        location = f"{start}..{end}"
        ret.append(["", "assembly_gap", location, "estimated_length", length_value])
        ret.append(["", "", "", "gap_type", GAP_TYPE])
        ret.append(["", "", "", "linkage_evidence", LINKAGE_EVIDENCE])

    return ret


def fix_partial_location(locations, left_partial, right_partial, strand):
    # locations must be sorted depending on the strand (1: ascending order of start pos, -1: descending order of end pos)
    first_location = locations[0]
    if strand == 1 and left_partial:
        new_location = FeatureLocation(BeforePosition(first_location.start), first_location.end, strand=first_location.strand)
        # sys.stderr.write(f"Location of 5' end fixed {first_location} --> {new_location}\n")
        locations[0] = new_location 
    elif strand == -1 and right_partial:
        new_location = FeatureLocation(first_location.start, AfterPosition(first_location.end), strand=first_location.strand)
        # sys.stderr.write(f"Location of 5' end fixed {right_location} --> {new_location}\n")
        locations[0] = new_location 

    last_location = locations[-1]
    if strand == 1 and right_partial:
        new_location = FeatureLocation(last_location.start, AfterPosition(last_location.end), strand=last_location.strand)
        # sys.stderr.write(f"Location of 5' end fixed {last_location} --> {new_location}\n")
        locations[-1] = new_location 
    elif strand == -1 and left_partial:
        new_location = FeatureLocation(BeforePosition(last_location.start), last_location.end, strand=last_location.strand)
        # sys.stderr.write(f"Location of 3' end fixed {last_location} --> {new_location}\n")
        locations[-1] = new_location



def create_mRNA_feature(mRNA_feature, gene_id, locus_tag, seq):
    """
    mRNAをMSSファイルに追加する。
    """

    def _check_partial_mRNA(mRNA_feature):
        # mRNAがpartialであるか判定。UTRの有無で判定を行う
        # ただし GFF ファイルにUTRが含まれない場合を想定し、エクソンとCDSの開始・終了位置で判断
        # エクソンの開始位置がCDSの開始位置より左なら left_partial = False
        # エクソンの終了位置がCDSの終了位置より右なら right_partial = False
        # trans_splicingには対応していない
        exon_locations = [sf.location for sf in mRNA_feature.sub_features if sf.type == "exon"]
        CDS_locations = [sf.location for sf in mRNA_feature.sub_features if sf.type == "CDS"]
        if len(exon_locations) == 0 or len(CDS_locations) == 0:
            sys.stderr.write(f"mRNA {mRNA_feature.id} has no exon or CDS feature. Could not determine if it is partial.\n")
            return False, False
        exon_start = min(exon_locations, key=lambda x: x.start).start
        cds_start = min(CDS_locations, key=lambda x: x.start).start
        assert exon_start <= cds_start
        left_partial = True if exon_start == cds_start else False
        exon_end = max(exon_locations, key=lambda x: x.end).end
        cds_end = max(CDS_locations, key=lambda x: x.end).end
        assert exon_end >= cds_end
        right_partial = True if exon_end == cds_end else False
        return left_partial, right_partial


    unwanted_qualifiers = ["ID", "Name", "Parent", "product", "inference", "note", "source", "db_xref"]
    qualifier_type_mapping = {
        "Dbxref": "db_xref",
        "Note": "note"
    }
    ret = []

    # mRNAのlocatioは exon を連結したものとする。ただし、exonがGFFに記載されていない場合 CDS で代用する場合がある
    exon_locations = [sf.location for sf in mRNA_feature.sub_features if sf.type == "exon"] 

    strand = mRNA_feature.location.strand
    assert all([location.strand == strand for location in exon_locations])

    # sort locations
    if strand == 1:
        exon_locations = sorted(exon_locations, key=lambda x: x.start)
    else:  # for minus strand
        exon_locations = sorted(exon_locations, key=lambda x: x.start, reverse=True)

    # fix partial mRNA
    left_partial, right_partial = _check_partial_mRNA(mRNA_feature)
    fix_partial_location(exon_locations, left_partial, right_partial, strand)

    compound_location = sum(exon_locations)
    location = _insdc_location_string(compound_location, len(seq))

    ret.append(["", "mRNA", location, "locus_tag", locus_tag])
    for qualifier, values in mRNA_feature.qualifiers.items():
        qualifier = qualifier_type_mapping.get(qualifier, qualifier)
        if qualifier in unwanted_qualifiers:
            continue
        for value in values:
            ret.append(["", "", "", qualifier, value])
    ret.append(add_note_describing_submitter_id(gene_id, mRNA_feature.id))
    return ret


def create_CDS_feature(mRNA_feature, gene_id, locus_tag, seq):
    unwanted_qualifiers = ["ID", "Name", "Parent", "source", "note", "db_xref"]
    qualifier_type_mapping = {
        "Dbxref": "db_xref",
        "Note": "note"
    }
    ret = []

    # set product name
    product = mRNA_feature.qualifiers.get("product", ["hypothetical protein"])[0]
    if "gene" in mRNA_feature.qualifiers:
        gene = mRNA_feature.qualifiers["gene"][0]
        if product == "hypothetical protein":
            product = f"protein {gene}"
        # else:
        #     product = f"{product} {gene}"
    mRNA_feature.qualifiers["product"] = [product]

    cds_locations = [sf.location for sf in mRNA_feature.sub_features if sf.type == "CDS"] 
    strand = mRNA_feature.location.strand
    assert all([location.strand == strand for location in cds_locations])

    # sort locations
    if strand == 1:
        cds_locations = sorted(cds_locations, key=lambda x: x.start)
    else:  # for minus strand
        cds_locations = sorted(cds_locations, key=lambda x: x.start, reverse=True)


    # check partial CDS
    compound_location = sum(cds_locations)
    cds_seq = compound_location.extract(seq)
    first_codon = str(cds_seq[:3]).upper()
    last_codon = str(cds_seq[-3:]).upper()
    five_prime_partial = True if first_codon not in start_codons else False
    three_prime_partial = True if last_codon not in stop_codons else False

    protein_seq = cds_seq.translate()

    # fix partial CDS
    if five_prime_partial or three_prime_partial:
        sys.stderr.write(f"CDS is partial: {mRNA_feature.id} 5' partial: {five_prime_partial}, 3' partial: {three_prime_partial}\n")
    if strand == 1:
        fix_partial_location(cds_locations, five_prime_partial, three_prime_partial, strand)
    else:
        fix_partial_location(cds_locations, three_prime_partial, five_prime_partial, strand)
    compound_location = sum(cds_locations)


    location = _insdc_location_string(compound_location, len(seq))

    ret.append(["", "CDS", location, "locus_tag", locus_tag])
    ret.append(["", "", "", "transl_table", TRANSL_TABLE])
    ret.append(["", "", "", "codon_start", "1"])  # TODO: codon_start should be set based on phase value of GFF
    for qualifier, values in mRNA_feature.qualifiers.items():
        qualifier = qualifier_type_mapping.get(qualifier, qualifier)
        if qualifier in unwanted_qualifiers:
            continue
        for value in values:
            ret.append(["", "", "", qualifier, value])
    ret.append(add_note_describing_submitter_id(gene_id, mRNA_feature.id))
    return ret

def create_UTR_feature(mRNA_feature, end, gene_id, locus_tag, seq):
    # end: "left" for 5' UTR of plus-strand CDS and 3' UTR of minus-strand CDS
    # "right" for 3' UTR of plus-strand CDS and 5' UTR of minus-strand CDS
    # TODO: GFFにが書かれていない場合がある。その場合には、CDSとexonからUTRを推定する必要がある
    strand = mRNA_feature.strand

    if strand == 1 and end == "left":
        target_features = ["five_prime_UTR", "5'UTR"]
        UTR_locations = [sf.location for sf in mRNA_feature.sub_features if sf.type in target_features]
        UTR_locations = sorted(UTR_locations, key=lambda x: x.start)
        feature_name = "5'UTR"
    elif strand == 1 and end == "right":
        target_features = ["three_prime_UTR", "3'UTR"]
        UTR_locations = [sf.location for sf in mRNA_feature.sub_features if sf.type in target_features]
        UTR_locations = sorted(UTR_locations, key=lambda x: x.start)
        feature_name = "3'UTR"
    elif strand == -1 and end == "left":
        target_features = ["three_prime_UTR", "3'UTR"]
        UTR_locations = [sf.location for sf in mRNA_feature.sub_features if sf.type in target_features]
        UTR_locations = sorted(UTR_locations, key=lambda x: x.start, reverse=True)
        feature_name = "3'UTR"
    elif strand == -1 and end == "right":
        target_features = ["five_prime_UTR", "5'UTR"]
        UTR_locations = [sf.location for sf in mRNA_feature.sub_features if sf.type in target_features]
        UTR_locations = sorted(UTR_locations, key=lambda x: x.start, reverse=True)
        feature_name = "5'UTR"
    else:
        raise ValueError(f"strand and end is value: {strand}")

    assert all([location.strand == strand for location in UTR_locations])

    compound_location = sum(UTR_locations)

    ret = []
    if compound_location:
        location = _insdc_location_string(compound_location, len(seq))
        ret.append(["", feature_name, location, "locus_tag", locus_tag])
        if "gene" in mRNA_feature.qualifiers:
            ret.append(["", "", "", "gene", mRNA_feature.qualifiers["gene"][0]])
        ret.append(add_note_describing_submitter_id(gene_id, mRNA_feature.id))
    return ret

def create_other_RNA_features(rna_feature, gene_id, locus_tag, seq):
    """
    for other biological features such as rRNA, tRNA, ncRNA, etc. (under development)
    """
    canonical_feature_names = ["misc_RNA", "ncRNA", "precursor_RNA", "rRNA", "tmRNA", "tRNA"]
    # tRNAとrRNAは別の扱いにした方が良さそう
    feature_name_mapping = {
        # RNAs not listed here will be treated as "misc_RNA". TO be updated
        "miRNA": "ncRNA",
        "pre_miRNA": "precursor_RNA",
        "ncRNA": "ncRNA",
        "snRNA": "ncRNA",
        "snoRNA": "ncRNA",
        "tRNA": "tRNA",
        "rRNA": "rRNA",
        "tmRNA": "tmRNA"
    }

    # unwanted_qualifiers = ["ID", "Name", "Parent", "product", "inference", "note", "Note", "source", "Dbxref", "db_xref"]
    allowed_qualifiers = ["note", "gene"]
    qualifier_type_mapping = {
        "Dbxref": "db_xref",
        "Note": "note"
    }

    # location は exon を連結したものとする。したがって、exonが child feature としてGFFに記載されていることを期待する 
    exon_locations = [sf.location for sf in rna_feature.sub_features if sf.type == "exon"] 

    strand = rna_feature.location.strand
    assert all([location.strand == strand for location in exon_locations])

    # sort locations
    if strand == 1:
        exon_locations = sorted(exon_locations, key=lambda x: x.start)
    else:  # for minus strand
        exon_locations = sorted(exon_locations, key=lambda x: x.start, reverse=True)


    compound_location = sum(exon_locations)
    location = _insdc_location_string(compound_location, len(seq))

    ret = []
    feature_name = feature_name_mapping.get(rna_feature.type, "misc_RNA")
    ret.append(["", feature_name, location, "locus_tag", locus_tag])
    if feature_name == "ncRNA":
        ret.append(["", "", "", "ncRNA_class", rna_feature.type])
    for qualifier, values in rna_feature.qualifiers.items():
        # if qualifier in unwanted_qualifiers:
        #     continue
        if qualifier in qualifier_type_mapping:
            qualifier = qualifier_type_mapping[qualifier]
        if qualifier not in allowed_qualifiers:
            continue

        for value in values:
            ret.append(["", "", "", qualifier, value])
    ret.append(add_note_describing_submitter_id(gene_id, rna_feature.id))
    return ret


seq_records = read_fasta_and_gff(input_fasta_file, input_gff_file)
ret = []
for record in seq_records:
    record.annotations = {}
    record.description = ""
    record.name = ""
    ret.extend(create_source_feature(record, souce_qualifiers))
    ret.extend(create_gap_feature(str(record.seq), length_estimated=False))
    for feature in record.features:  # gene level feature
        if feature.type != "gene":
            continue
        locus_tag = get_locus_tag(feature)
        gene_id = feature.id  # todo : GFF の　ID か Nameから取得するように変更した方が良いかも
        for sf in feature.sub_features:  # mRNA as well as other RNA level features
            if sf.type == "mRNA":
                ret.extend(create_mRNA_feature(sf, gene_id, locus_tag, record.seq))
                ret.extend(create_UTR_feature(sf, "left", gene_id, locus_tag, record.seq))
                ret.extend(create_CDS_feature(sf, gene_id, locus_tag, record.seq))
                ret.extend(create_UTR_feature(sf, "right", gene_id, locus_tag, record.seq))
                # print(create_mRNA_feature(sf, len(record)))
            else:
                ret.extend(create_other_RNA_features(sf, gene_id, locus_tag, record.seq))


# output
print()
with open(output_ann_file, "w") as f:
    common = open(common_metadata_file).read()
    f.write(common)
    for line in ret:
        f.write("\t".join(line) + "\n")

with open(output_seq_file, "w") as f:
    for record in seq_records:
        f.write(record.format("fasta"))
        f.write("//\n")

