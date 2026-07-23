# flatfile2gff Implementation Plan (nuclear 3-level)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert a DDBJ flatfile (INSDC feature table) into the project's canonical INSDC GFF3 for **nuclear eukaryotes (gene→mRNA→exon/CDS, 3-level)**, synthesizing the gene and exon features the flatfile lacks, so the biological features round-trip via `gff2mss`→`mss2ff`.

**Architecture:** New package `src/ddbj_gff/flatfile/` in `ddbj-gff`. Reads the flatfile with BioPython's genbank parser, detects molecule type (taxid/division/topology/compartment) to pick the hierarchy, synthesizes gene (grouped by locus_tag) + mRNA + exon, maps qualifiers/locations, assembles a `GffDocument`, runs the existing `normalize` (adds directives; no-op for the biology passes here), and writes with the existing `writer`. Round-trip is validated in `ddbj_mss_tools` via the existing `gff2mss`+`mss2ff`.

**Tech Stack:** Python 3.10+, BioPython 1.87 (`SeqIO.read(path,"genbank")`), pytest 9. All tests run inside the `ddbj-gff-dev` Docker container.

## Global Constraints

- **Two repos, one-way dependency:** the flatfile2gff module lives in `ddbj-gff` (`gff_submission`) and MUST NOT import `gff2mss`/`common`. The round-trip test (Task 5) lives in `ddbj_mss_tools` and may import `ddbj_gff`.
- **Output = project canonical INSDC GFF** (feature-mapping conformant, directives present, passes `ddbj_gff.validate`). **Do NOT emulate NCBI's GFF conventions** (NCBI's gene granularity differs — 33 genes vs 30 locus_tags for this record; the INSDC-GFF spec is authoritative).
- **Gene grouping key = `/locus_tag`.** One `gene` per distinct locus_tag.
- **CDS↔mRNA pairing = coordinate containment** (every CDS segment lies within an mRNA exon on the same strand, boundaries compatible). Ambiguous → best-boundary-match + diagnostic; none → synthesize mRNA = CDS location + diagnostic.
- **exon synthesis:** one exon per mRNA segment; if a transcript has no mRNA feature, synthesize `mRNA` = CDS location first, then exons.
- **codon_start → phase:** first biological CDS segment phase = `codon_start-1`; propagate cumulatively. Drop `/translation` (regenerated downstream) and `/codon_start` (→ phase). Drop `assembly_gap` features (regenerated downstream).
- **COMMON ignored** (submitter/reference/DBLINK/SRA). Extract only: taxid, organism, division, topology, source compartment.
- **Molecule type (this plan):** nuclear (source has no `/organelle`) ⇒ 3-level, `transl_table` from the CDS qualifier. Organelle/trans-splicing deferred.
- **Commit policy:** `gff_submission` (Tasks 1-4) = normal commits on branch `feat/flatfile2gff`. `ddbj_mss_tools` (Task 5) = commit ONLY `tests/**`; NEVER `git add -A` (~49 pre-existing uncommitted files); do NOT push.
- **Reference spec:** `docs/superpowers/specs/2026-07-05-flatfile2gff-design.md`.
- **Fixture (already created, in repo):** `tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk` — a 4000 bp re-coordinatized excerpt of *Citrus unshiu* `BDQV01000200.1` (DDBJ getentry). Contains: source (1..4000, `/db_xref=taxon:55188`, `/organism="Citrus unshiu"`, `/mol_type="genomic DNA"`, `/submitter_seqid="scaffold00200"`, `/cultivar`), 3 mRNA, 3 CDS over **2 locus_tags** — `CUMW_191330` (1 transcript) and `CUMW_191340` (**2 transcripts = alt-splicing**). division PLN, topology linear.

### Verified fixture facts (pin these in tests)

| locus_tag | transcript | mRNA (1-based) | CDS | protein (table 1) |
|---|---|---|---|---|
| CUMW_191330 | single | 1216..2597 (2 seg) | 1229..2241 (2 seg) | `MDLLINCILWLVFTL…` 314 aa, 0 internal stop |
| CUMW_191340 | A | 1229..2597 (3 seg) | 1229..2241 (3 seg) | `MDLLINCILWLVFTL…` 270 aa, 0 stop |
| CUMW_191340 | B | 2245..2762 (2 seg) | 2245..2762 (2 seg) | `MAELLHNPEALLKAK…` 170 aa, 0 stop |

All CDS: `transl_table=1`, `codon_start=1`, have `/protein_id`. Expected flatfile2gff output: **1 region + 2 gene + 3 mRNA + 7 exon (2+3+2) + 3 CDS**.

## Test Environment (`ddbj-gff-dev` container)

- `ddbj_gff` bind-mounted live at `/workspace/src`; tests at `/workspace/tests` (Tasks 1-4 need NO sync). Fixture at `/workspace/tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk`.
- `gff2mss`+`common` at `/opt/mss_src` (docker-cp'd — re-sync BOTH before Task 5, the container `common` goes stale). ddbj_mss_tools tests at `/opt/mss_tests`.
- venv `/opt/ddbj-venv` (Bio 1.87, pytest 9.1).

**Run ddbj-gff tests (Tasks 1-4):**
```bash
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev \
  bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/<file> -v'
```

**Sync + run Task 5 (ddbj_mss_tools):**
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/common/. ddbj-gff-dev:/opt/mss_src/common
docker cp /Users/tanizawa/projects/ddbj/gff_submission/tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk ddbj-gff-dev:/opt/mss_tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/<file> ddbj-gff-dev:/opt/mss_tests/<file>
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest <file> -v'
```

## File Structure

- `src/ddbj_gff/flatfile/__init__.py` — exports `flatfile_to_gff`, `detect_molecule`.
- `src/ddbj_gff/flatfile/molecule.py` — `MoleculeInfo` + `detect_molecule(rec)`.
- `src/ddbj_gff/flatfile/convert.py` — location/qualifier mappers, gene/transcript/exon synthesis, `flatfile_to_gff(rec, molecule)`.
- `src/ddbj_gff/flatfile/cli.py` — `flatfile2gff` CLI.
- Tests: `tests/test_flatfile_molecule.py`, `tests/test_flatfile_mapping.py`, `tests/test_flatfile_synthesis.py`, `tests/test_flatfile_convert.py` (gff_submission); `ddbj_mss_tools/tests/test_flatfile_roundtrip.py` (Task 5).

---

### Task 1: molecule-type detection

**Files:**
- Create: `src/ddbj_gff/flatfile/__init__.py` (start with `detect_molecule` export), `src/ddbj_gff/flatfile/molecule.py`
- Test: `tests/test_flatfile_molecule.py`

**Interfaces:**
- Produces: `MoleculeInfo` dataclass with fields `taxid:int|None, organism:str|None, division:str|None, topology:str, compartment:str, hierarchy:str, transl_table:int` and `detect_molecule(rec) -> MoleculeInfo`. `compartment` ∈ {`"nuclear"`,`"organelle"`}; `hierarchy` ∈ {`"three_level"`,`"two_level"`}. Consumed by Task 4.

- [ ] **Step 1: Write the failing test**

Create `tests/test_flatfile_molecule.py`:
```python
from Bio import SeqIO
from ddbj_gff.flatfile.molecule import detect_molecule

FIX = "tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk"


def test_detect_nuclear_plant():
    rec = SeqIO.read(FIX, "genbank")
    m = detect_molecule(rec)
    assert m.taxid == 55188
    assert m.organism == "Citrus unshiu"
    assert m.division == "PLN"
    assert m.topology == "linear"
    assert m.compartment == "nuclear"     # source has no /organelle
    assert m.hierarchy == "three_level"
    assert m.transl_table == 1            # from CDS /transl_table
```

- [ ] **Step 2: Run test to verify it fails**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_molecule.py -v'`
Expected: FAIL — `ModuleNotFoundError: No module named 'ddbj_gff.flatfile'`.

- [ ] **Step 3: Implement**

Create `src/ddbj_gff/flatfile/molecule.py`:
```python
from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass
class MoleculeInfo:
    taxid: int | None
    organism: str | None
    division: str | None
    topology: str          # "linear" | "circular"
    compartment: str       # "nuclear" | "organelle"
    hierarchy: str         # "three_level" | "two_level"
    transl_table: int


def _source_feature(rec):
    for f in rec.features:
        if f.type == "source":
            return f
    return None


def detect_molecule(rec) -> MoleculeInfo:
    src = _source_feature(rec)
    quals = src.qualifiers if src is not None else {}
    taxid = None
    for xref in quals.get("db_xref", []):
        m = re.match(r"taxon:(\d+)", xref)
        if m:
            taxid = int(m.group(1))
    organism = (quals.get("organism") or [rec.annotations.get("organism")])[0] \
        if (quals.get("organism") or rec.annotations.get("organism")) else None
    division = rec.annotations.get("data_file_division")
    topology = rec.annotations.get("topology") or "linear"
    compartment = "organelle" if quals.get("organelle") else "nuclear"
    hierarchy = "two_level" if compartment == "organelle" else "three_level"
    # transl_table: primary from a CDS qualifier, default 1
    tt = 1
    for f in rec.features:
        if f.type == "CDS" and f.qualifiers.get("transl_table"):
            tt = int(f.qualifiers["transl_table"][0])
            break
    return MoleculeInfo(taxid=taxid, organism=organism, division=division,
                        topology=topology, compartment=compartment,
                        hierarchy=hierarchy, transl_table=tt)
```
Create `src/ddbj_gff/flatfile/__init__.py`:
```python
"""flatfile2gff: DDBJ flatfile -> canonical INSDC GFF3."""
from .molecule import MoleculeInfo, detect_molecule

__all__ = ["MoleculeInfo", "detect_molecule"]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_molecule.py -v'`
Expected: PASS (1 passed).

- [ ] **Step 5: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add src/ddbj_gff/flatfile/__init__.py src/ddbj_gff/flatfile/molecule.py tests/test_flatfile_molecule.py
git commit -m "feat(flatfile): detect_molecule (taxid/division/topology/compartment/hierarchy)"
```

---

### Task 2: location + qualifier mapping helpers

**Files:**
- Create: `src/ddbj_gff/flatfile/convert.py`
- Test: `tests/test_flatfile_mapping.py`

**Interfaces:**
- Consumes: `ddbj_gff.model.Span`.
- Produces: `bio_location_to_spans(location, seqid, *, is_cds, codon_start=1) -> list[Span]` (BioPython location parts → 1-based `Span`s in biological order; per-segment `phase` for CDS, else `None`); `qualifiers_to_attrs(feature) -> dict[str,list[str]]` (INSDC qualifier → GFF attribute; drops `translation`/`codon_start`). Consumed by Task 3.

- [ ] **Step 1: Write the failing test**

Create `tests/test_flatfile_mapping.py`:
```python
from Bio import SeqIO
from ddbj_gff.flatfile.convert import bio_location_to_spans, qualifiers_to_attrs

FIX = "tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk"


def _feat(rec, ftype, lt, seg):
    for f in rec.features:
        if f.type == ftype and f.qualifiers.get("locus_tag", [None])[0] == lt \
                and len(f.location.parts) == seg:
            return f
    raise AssertionError("feature not found")


def test_cds_location_to_spans_1based_phase():
    rec = SeqIO.read(FIX, "genbank")
    cds = _feat(rec, "CDS", "CUMW_191330", 2)          # 1229..2241, 2 segments, + strand
    spans = bio_location_to_spans(cds.location, "S", is_cds=True, codon_start=1)
    assert len(spans) == 2
    assert (spans[0].start, spans[0].strand, spans[0].phase) == (1229, "+", 0)  # codon_start 1 -> phase 0
    assert all(s.strand == "+" for s in spans)


def test_mrna_location_no_phase():
    rec = SeqIO.read(FIX, "genbank")
    mrna = _feat(rec, "mRNA", "CUMW_191340", 3)        # 3 exons
    spans = bio_location_to_spans(mrna.location, "S", is_cds=False)
    assert len(spans) == 3
    assert all(s.phase is None for s in spans)


def test_qualifiers_to_attrs_maps_and_drops():
    rec = SeqIO.read(FIX, "genbank")
    cds = _feat(rec, "CDS", "CUMW_191330", 2)
    attrs = qualifiers_to_attrs(cds)
    assert attrs["locus_tag"] == ["CUMW_191330"]
    assert attrs["product"] == ["hypothetical protein"]
    assert "protein_id" in attrs and "transl_table" in attrs
    assert "translation" not in attrs and "codon_start" not in attrs   # dropped
    assert "Note" in attrs or "note" not in cds.qualifiers             # /note -> Note
```

- [ ] **Step 2: Run test to verify it fails**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_mapping.py -v'`
Expected: FAIL — `ImportError` (convert.py absent).

- [ ] **Step 3: Implement**

Create `src/ddbj_gff/flatfile/convert.py`:
```python
from __future__ import annotations

from ..model import Span

_STRAND = {1: "+", -1: "-", 0: "?", None: "."}

# INSDC qualifier -> canonical GFF attribute key
_QUAL_MAP = {
    "gene": "gene", "locus_tag": "locus_tag", "product": "product",
    "note": "Note", "protein_id": "protein_id", "gene_synonym": "gene_synonym",
    "transl_table": "transl_table", "db_xref": "Dbxref", "pseudo": "pseudo",
    "ncRNA_class": "ncRNA_class",
}
_DROP_QUALS = {"translation", "codon_start"}


def bio_location_to_spans(location, seqid, *, is_cds, codon_start=1):
    """BioPython location -> 1-based Spans in biological (5'->3') order.
    parts come from BioPython in transcription order. For CDS, per-segment phase
    is derived from codon_start; else phase is None."""
    parts = list(location.parts)
    spans = []
    phase = (codon_start - 1) if is_cds else None
    for p in parts:
        start = int(p.start) + 1
        end = int(p.end)
        strand = _STRAND.get(p.strand, ".")
        spans.append(Span(seqid, start, end, strand, phase=phase))
        if is_cds:
            seg_len = end - start + 1
            phase = (3 - ((seg_len - phase) % 3)) % 3
    return spans


def qualifiers_to_attrs(feature) -> dict:
    attrs: dict[str, list[str]] = {}
    for k, vals in feature.qualifiers.items():
        if k in _DROP_QUALS:
            continue
        gk = _QUAL_MAP.get(k, k)
        attrs[gk] = list(vals)
    return attrs
```

- [ ] **Step 4: Run test to verify it passes**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_mapping.py -v'`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add src/ddbj_gff/flatfile/convert.py tests/test_flatfile_mapping.py
git commit -m "feat(flatfile): bio_location_to_spans (phase from codon_start) + qualifiers_to_attrs"
```

---

### Task 3: gene / transcript / exon synthesis

**Files:**
- Modify: `src/ddbj_gff/flatfile/convert.py` (add synthesis functions)
- Test: `tests/test_flatfile_synthesis.py`

**Interfaces:**
- Consumes: `bio_location_to_spans`, `qualifiers_to_attrs` (Task 2), `ddbj_gff.model.Feature`, `Span`.
- Produces: `synthesize_features(rec, seqid) -> list[Feature]` — returns canonical Features: one `gene` per locus_tag (Parent-less), each mRNA (Parent=gene), synthesized `exon`s (Parent=mRNA), each CDS (Parent=mRNA), with `ID`/`Parent` wired. CDS↔mRNA paired by containment. Consumed by Task 4.

- [ ] **Step 1: Write the failing test**

Create `tests/test_flatfile_synthesis.py`:
```python
from Bio import SeqIO
from collections import Counter
from ddbj_gff.flatfile.convert import synthesize_features

FIX = "tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk"


def test_synthesis_counts_and_parentage():
    rec = SeqIO.read(FIX, "genbank")
    feats = synthesize_features(rec, "BDQV01000200.1")
    counts = Counter(f.type for f in feats)
    assert counts["gene"] == 2
    assert counts["mRNA"] == 3
    assert counts["exon"] == 7      # 2 + 3 + 2
    assert counts["CDS"] == 3

    by_id = {f.id: f for f in feats}
    genes = [f for f in feats if f.type == "gene"]
    gene_ids = {g.id for g in genes}
    # every mRNA parents a gene; every exon/CDS parents an mRNA
    mrnas = [f for f in feats if f.type == "mRNA"]
    for m in mrnas:
        assert m.parent_ids and m.parent_ids[0] in gene_ids
    mrna_ids = {m.id for m in mrnas}
    for f in feats:
        if f.type in ("exon", "CDS"):
            assert f.parent_ids and f.parent_ids[0] in mrna_ids


def test_altsplice_cds_paired_by_containment():
    rec = SeqIO.read(FIX, "genbank")
    feats = synthesize_features(rec, "BDQV01000200.1")
    mrnas = {f.id: f for f in feats if f.type == "mRNA"}
    # CUMW_191340 has two transcripts: CDS-A (ends 2241) pairs with mRNA spanning ..2597,
    # CDS-B (2245..2762) pairs with the mRNA spanning 2245..2762
    cds = [f for f in feats if f.type == "CDS" and f._first("locus_tag") == "CUMW_191340"]
    assert len(cds) == 2
    for c in cds:
        parent = mrnas[c.parent_ids[0]]
        cs = c.ordered_spans()
        ms = parent.ordered_spans()
        lo_m, hi_m = min(s.start for s in ms), max(s.end for s in ms)
        # CDS fully within its paired mRNA
        assert lo_m <= min(s.start for s in cs) and max(s.end for s in cs) <= hi_m
    # the two CDS are paired to DIFFERENT mRNAs
    assert cds[0].parent_ids[0] != cds[1].parent_ids[0]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_synthesis.py -v'`
Expected: FAIL — `ImportError: cannot import name 'synthesize_features'`.

- [ ] **Step 3: Implement**

Append to `src/ddbj_gff/flatfile/convert.py`:
```python
from collections import OrderedDict
from ..model import Feature

_RNA_TYPES = {"mRNA", "tRNA", "rRNA", "ncRNA", "misc_RNA"}
_BIOTYPE = {"CDS": "protein_coding", "mRNA": "protein_coding",
            "tRNA": "tRNA", "rRNA": "rRNA", "ncRNA": "ncRNA"}


def _locus_tag(f):
    return f.qualifiers.get("locus_tag", [None])[0]


def _cds_within(cds_spans, mrna_spans) -> bool:
    for cs in cds_spans:
        if not any(ms.strand == cs.strand and ms.start <= cs.start and cs.end <= ms.end
                   for ms in mrna_spans):
            return False
    return True


def _shared_boundaries(cds_spans, mrna_spans) -> int:
    edges = {(s.start) for s in mrna_spans} | {(s.end) for s in mrna_spans}
    return sum((s.start in edges) + (s.end in edges) for s in cds_spans)


def synthesize_features(rec, seqid) -> list:
    """Group flatfile mRNA/CDS/RNA by locus_tag into canonical gene->mRNA->exon/CDS."""
    bio = [f for f in rec.features if f.type in _RNA_TYPES or f.type == "CDS"]
    groups = OrderedDict()
    for f in bio:
        groups.setdefault(_locus_tag(f), []).append(f)

    out: list = []
    for lt, members in groups.items():
        gene_id = f"gene-{lt}" if lt else f"gene-{len(out)}"
        mrnas = [f for f in members if f.type == "mRNA"]
        cdss = [f for f in members if f.type == "CDS"]
        rnas = [f for f in members if f.type in _RNA_TYPES and f.type != "mRNA"]

        # pair each CDS to a containing mRNA; synthesize an mRNA if none
        transcripts = []  # list of (mrna_feature_or_None_synth, [cds])
        used = {id(m): [] for m in mrnas}
        for c in cdss:
            cspans = bio_location_to_spans(c.location, seqid, is_cds=True)
            cand = [m for m in mrnas
                    if _cds_within(cspans, bio_location_to_spans(m.location, seqid, is_cds=False))]
            if cand:
                best = max(cand, key=lambda m: _shared_boundaries(
                    cspans, bio_location_to_spans(m.location, seqid, is_cds=False)))
                used[id(best)].append(c)
            else:
                transcripts.append((None, [c]))   # synth mRNA = CDS
        for m in mrnas:
            transcripts.append((m, used[id(m)]))
        for r in rnas:                              # tRNA/rRNA: transcript is the RNA itself
            transcripts.append((r, []))

        # gene span = union of all member spans
        all_spans = []
        for m in members:
            all_spans += bio_location_to_spans(m.location, seqid, is_cds=(m.type == "CDS"))
        g_lo, g_hi = min(s.start for s in all_spans), max(s.end for s in all_spans)
        g_strand = all_spans[0].strand
        biotype = _BIOTYPE.get(members[0].type, "other")
        gene_attrs = {"ID": [gene_id]}
        if _locus_tag(members[0]):
            gene_attrs["locus_tag"] = [_locus_tag(members[0])]
        if members[0].qualifiers.get("gene"):
            gene_attrs["gene"] = list(members[0].qualifiers["gene"])
            gene_attrs["Name"] = list(members[0].qualifiers["gene"])
        if members[0].qualifiers.get("gene_synonym"):
            gene_attrs["gene_synonym"] = list(members[0].qualifiers["gene_synonym"])
        gene_attrs["gene_biotype"] = [biotype]
        out.append(Feature(gene_id, "DDBJ", "gene", [Span(seqid, g_lo, g_hi, g_strand)],
                           gene_attrs, []))

        for i, (mfeat, member_cds) in enumerate(transcripts, 1):
            is_rna = mfeat is not None and mfeat.type in _RNA_TYPES and mfeat.type != "mRNA"
            tx_type = mfeat.type if is_rna else "mRNA"
            tx_id = f"{('rna' if is_rna else 'mrna')}-{lt}-{i}"
            if mfeat is not None:
                tx_spans = bio_location_to_spans(mfeat.location, seqid, is_cds=False)
                tx_attrs = qualifiers_to_attrs(mfeat)
            else:                                   # synth mRNA from the CDS
                cspans = bio_location_to_spans(member_cds[0].location, seqid, is_cds=True)
                tx_spans = [Span(seqid, s.start, s.end, s.strand) for s in cspans]
                tx_attrs = {k: v for k, v in qualifiers_to_attrs(member_cds[0]).items()
                            if k in ("locus_tag", "gene", "product", "Note")}
            tx_attrs["ID"] = [tx_id]
            tx_attrs["Parent"] = [gene_id]
            out.append(Feature(tx_id, "DDBJ", tx_type, tx_spans, tx_attrs, [gene_id]))

            if not is_rna:                          # exons for mRNA transcripts
                for j, sp in enumerate(tx_spans, 1):
                    ex_id = f"exon-{lt}-{i}-{j}"
                    out.append(Feature(ex_id, "DDBJ", "exon",
                                       [Span(seqid, sp.start, sp.end, sp.strand)],
                                       {"ID": [ex_id], "Parent": [tx_id]}, [tx_id]))
            for c in member_cds:
                cspans = bio_location_to_spans(c.location, seqid, is_cds=True,
                    codon_start=int(c.qualifiers.get("codon_start", ["1"])[0]))
                c_attrs = qualifiers_to_attrs(c)
                c_id = f"cds-{c.qualifiers.get('protein_id', [tx_id])[0]}"
                c_attrs["ID"] = [c_id]
                c_attrs["Parent"] = [tx_id]
                out.append(Feature(c_id, "DDBJ", "CDS", cspans, c_attrs, [tx_id]))
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_synthesis.py -v'`
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add src/ddbj_gff/flatfile/convert.py tests/test_flatfile_synthesis.py
git commit -m "feat(flatfile): synthesize gene/mRNA/exon; pair CDS<->mRNA by containment"
```

---

### Task 4: assemble GffDocument + directives + CLI (validate-clean)

**Files:**
- Modify: `src/ddbj_gff/flatfile/convert.py` (add `flatfile_to_gff`), `src/ddbj_gff/flatfile/__init__.py` (export it)
- Create: `src/ddbj_gff/flatfile/cli.py`
- Test: `tests/test_flatfile_convert.py`

**Interfaces:**
- Consumes: `synthesize_features`, `detect_molecule`, `ddbj_gff.normalize.normalize`, `ddbj_gff.normalize.config.NormalizeConfig`, `ddbj_gff.writer.write`, `ddbj_gff.validate.validate`, `ddbj_gff.model.{GffDocument,Feature,Span}`.
- Produces: `flatfile_to_gff(rec) -> GffDocument` (region feature + synthesized features, run through `normalize` to add directives); `cli.main()` reads a `.gbk` and writes GFF3.

- [ ] **Step 1: Write the failing test**

Create `tests/test_flatfile_convert.py`:
```python
from Bio import SeqIO
from ddbj_gff.flatfile.convert import flatfile_to_gff
from ddbj_gff.writer import write
from ddbj_gff.validate import validate

FIX = "tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk"


def test_flatfile_to_gff_validates_and_has_hierarchy():
    rec = SeqIO.read(FIX, "genbank")
    doc = flatfile_to_gff(rec)
    types = [f.type for f in doc.features]
    assert types.count("gene") == 2 and types.count("mRNA") == 3
    assert types.count("exon") == 7 and types.count("CDS") == 3
    assert any(f.type == "region" for f in doc.features)
    # directives added by normalize
    text = write(doc)
    assert "##sequence-region" in text
    assert "id=55188" in text            # ##species taxid
    assert "##gff-version 3" in text
    # canonical GFF validates with no ERROR-level diagnostics
    diags = validate(doc)
    errors = [d for d in diags if getattr(d, "severity", None) and d.severity.name == "ERROR"]
    assert errors == [], f"unexpected validate errors: {[d.code for d in errors]}"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_convert.py -v'`
Expected: FAIL — `ImportError: cannot import name 'flatfile_to_gff'`.

- [ ] **Step 3: Implement `flatfile_to_gff`**

Append to `src/ddbj_gff/flatfile/convert.py`:
```python
from ..model import GffDocument
from ..normalize.normalize import normalize
from ..normalize.config import NormalizeConfig
from .molecule import detect_molecule


def _region_feature(rec, mol, seqid):
    src = next((f for f in rec.features if f.type == "source"), None)
    attrs = {"ID": [f"{seqid}:1..{len(rec.seq)}"]}
    if mol.taxid:
        attrs["Dbxref"] = [f"taxon:{mol.taxid}"]
    if mol.topology == "circular":
        attrs["Is_circular"] = ["true"]
    if src is not None:
        for k in ("mol_type", "organism", "submitter_seqid", "chromosome", "organelle"):
            if src.qualifiers.get(k):
                attrs[k] = list(src.qualifiers[k])
    return Feature(attrs["ID"][0], "DDBJ", "region",
                   [Span(seqid, 1, len(rec.seq), "+")], attrs, [])


def flatfile_to_gff(rec) -> GffDocument:
    mol = detect_molecule(rec)
    seqid = rec.id
    feats = [_region_feature(rec, mol, seqid)] + synthesize_features(rec, seqid)
    doc = GffDocument(directives=[], features=feats,
                      fasta={seqid: str(rec.seq)})
    norm, _report = normalize(doc, seq_lengths={seqid: len(rec.seq)},
                              config=NormalizeConfig(taxid=mol.taxid,
                                                     transl_table=mol.transl_table))
    return norm
```
Update `src/ddbj_gff/flatfile/__init__.py`:
```python
"""flatfile2gff: DDBJ flatfile -> canonical INSDC GFF3."""
from .molecule import MoleculeInfo, detect_molecule
from .convert import flatfile_to_gff

__all__ = ["MoleculeInfo", "detect_molecule", "flatfile_to_gff"]
```
Create `src/ddbj_gff/flatfile/cli.py`:
```python
from __future__ import annotations

import argparse

from Bio import SeqIO

from .convert import flatfile_to_gff
from ..writer import write


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="flatfile2gff",
                                 description="DDBJ flatfile -> canonical INSDC GFF3")
    ap.add_argument("--in", dest="infile", required=True, help="input DDBJ flatfile (.gbk)")
    ap.add_argument("--out", dest="outfile", required=True, help="output GFF3")
    args = ap.parse_args(argv)
    rec = SeqIO.read(args.infile, "genbank")
    doc = flatfile_to_gff(rec)
    with open(args.outfile, "w", encoding="utf-8") as fh:
        fh.write(write(doc))
    print(f"[flatfile2gff] -> {args.outfile}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests/test_flatfile_convert.py -v'`
Expected: PASS. Note: the `region` landmark type is NOT in the SO-term vocabulary, so validate
emits `feature-type-not-insdc` for it — but that code is **WARNING**, not ERROR (verified), and
the established canonical fixtures (chloroplast.gff3, cp187952) use `region` the same way, so the
test (which gates on ERROR only) passes. If an actual ERROR appears, read the code: `undefined-seqid`
means the `##sequence-region` was not added (check `normalize` got `seq_lengths`); `cds-invalid-phase`
means `bio_location_to_spans` produced a bad phase; `dangling-parent` means an `ID`/`Parent` wiring
bug in Task 3. Fix the mapping, not the test.

- [ ] **Step 5: Full ddbj-gff regression (not-slow)**

Run: `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && PYTHONPATH=/workspace/src /opt/ddbj-venv/bin/python -m pytest tests -m "not slow" -q'`
Expected: PASS (prior 158 + the 4 new files' tests). Zero failures.

- [ ] **Step 6: Commit**

```bash
cd /Users/tanizawa/projects/ddbj/gff_submission
git add src/ddbj_gff/flatfile/convert.py src/ddbj_gff/flatfile/__init__.py src/ddbj_gff/flatfile/cli.py tests/test_flatfile_convert.py
git commit -m "feat(flatfile): flatfile_to_gff (region + normalize directives) + flatfile2gff CLI"
```

---

### Task 5: round-trip verification (ddbj_mss_tools)

**Files:**
- Create: `ddbj_mss_tools/tests/test_flatfile_roundtrip.py`
- Create: `ddbj_mss_tools/tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk` (copy of the ddbj-gff fixture)

**Interfaces:**
- Consumes: `ddbj_gff.flatfile.flatfile_to_gff`, `ddbj_gff.writer.write`, `ddbj_gff.parse`, `gff2mss.convert.convert`, `gff2mss.emit.emit_ann`, `gff2mss.config.MssConfig`, `mss2ff` (`ann_parser.parse_ann` + `ff_writer`, or reuse `build_ann_text`→file→`mss2ff.cli`). Simplest: assert on the gff2mss `.ann` text (which mss2ff renders verbatim), plus verify translation via BioPython.
- Produces: proof that CDS features survive `flatfile → GFF → gff2mss` with matching coordinates + translation.

- [ ] **Step 1: Copy the fixture into ddbj_mss_tools**

```bash
mkdir -p /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/flatfile_fixtures
cp /Users/tanizawa/projects/ddbj/gff_submission/tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk \
   /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk
```

- [ ] **Step 2: Write the round-trip test**

Create `ddbj_mss_tools/tests/test_flatfile_roundtrip.py`:
```python
import os
from Bio import SeqIO
from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.flatfile import flatfile_to_gff
from ddbj_gff.writer import write
from gff2mss.convert import convert
from gff2mss.emit import emit_ann
from gff2mss.config import MssConfig

FIX = os.path.join(os.path.dirname(__file__), "flatfile_fixtures", "citrus_unshiu_excerpt.gbk")


def _translations_from_flatfile(rec):
    out = []
    for f in rec.features:
        if f.type == "CDS":
            tt = int(f.qualifiers.get("transl_table", ["1"])[0])
            prot = str(f.extract(rec.seq).translate(table=tt))
            out.append(prot[:-1] if prot.endswith("*") else prot)
    return sorted(out)


def test_flatfile_to_gff_roundtrip_cds_translation():
    rec = SeqIO.read(FIX, "genbank")
    # forward: flatfile -> canonical GFF
    gff_text = write(flatfile_to_gff(rec))
    doc = parse(gff_text)
    seqs = {rec.id: rec.seq}
    cfg = MssConfig(source={}, transl_table=1, product_default="hypothetical protein")
    cfg.emit_mrna = True                                   # nuclear 3-level
    mss_doc, _ = convert(doc, seqs, cfg, common_rows=[])
    ann = emit_ann(mss_doc)

    # the .ann (which mss2ff renders verbatim to a flatfile) carries mRNA + CDS
    assert "\tCDS\t" in ann and "\tmRNA\t" in ann

    # CDS translations survive the loop: recompute from the round-tripped CDS locations
    # by re-parsing the GFF CDS spans and translating against the original sequence.
    from ddbj_gff.model import Feature
    cds_feats = [f for f in doc.features if f.type == "CDS"]
    got = []
    for c in cds_feats:
        loc = c.to_biopython_location()
        prot = str(loc.extract(rec.seq).translate(table=1))
        got.append(prot[:-1] if prot.endswith("*") else prot)
    assert sorted(got) == _translations_from_flatfile(rec)   # 3 CDS, translations match
```

- [ ] **Step 3: Sync + run**

Run:
```bash
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/gff2mss/. ddbj-gff-dev:/opt/mss_src/gff2mss
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/src/common/. ddbj-gff-dev:/opt/mss_src/common
docker exec ddbj-gff-dev bash -lc 'mkdir -p /opt/mss_tests/flatfile_fixtures'
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk ddbj-gff-dev:/opt/mss_tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk
docker cp /Users/tanizawa/projects/ddbj/ddbj_mss_tools/tests/test_flatfile_roundtrip.py ddbj-gff-dev:/opt/mss_tests/test_flatfile_roundtrip.py
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest test_flatfile_roundtrip.py -v'
```
Expected: PASS. The 3 CDS translations (`MDLLINCILWLVFTL…` ×2 from the two 191330/191340-A CDS and `MAELLHNPEALLKAK…` for 191340-B) match between the original flatfile and the round-tripped GFF CDS. If it fails on translation, the phase/codon_start or the multi-exon span order in `bio_location_to_spans` is wrong (Task 2/3) — fix there.

- [ ] **Step 4: Full gff2mss regression**

Run:
```bash
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev \
  bash -lc 'cd /opt/mss_tests && /opt/ddbj-venv/bin/python -m pytest . -q'
```
Expected: PASS (prior suite + new test; the pre-existing `@slow` skip). Zero failures.

- [ ] **Step 5: Commit (ddbj_mss_tools test + fixture only)**

```bash
cd /Users/tanizawa/projects/ddbj/ddbj_mss_tools
git add tests/test_flatfile_roundtrip.py tests/flatfile_fixtures/citrus_unshiu_excerpt.gbk
git commit -m "test(flatfile): round-trip DDBJ flatfile -> GFF -> gff2mss CDS translation match"
```

---

## Self-Review

**Spec coverage:** molecule detection (taxid/organism/division/topology/compartment/hierarchy/transl_table) → Task 1. Location→spans+phase & qualifier→attr mapping (drop translation/codon_start) → Task 2. gene synthesis by locus_tag + CDS↔mRNA containment pairing + exon synthesis (+ synth mRNA when absent) → Task 3. region feature + directives + canonical output + validate-clean + CLI → Task 4. Round-trip biological-feature (CDS translation) verification → Task 5. assembly_gap dropped (not in `_RNA_TYPES`/CDS collection, so never emitted) — Task 3. NCBI-not-emulated: grouping by locus_tag (Task 3) gives 2 genes for this fixture, independent of NCBI. Nuclear-3-level only; organelle deferred (spec Out of scope).

**Placeholder scan:** No TBD/TODO. Every code step is complete. The Task 4/5 contingency notes point at concrete diagnostics/functions to check, not vague instructions.

**Type consistency:** `MoleculeInfo` fields (Task 1) consumed in Task 4. `bio_location_to_spans(location, seqid, *, is_cds, codon_start=1) -> list[Span]` and `qualifiers_to_attrs(feature) -> dict` (Task 2) used identically in Task 3. `synthesize_features(rec, seqid) -> list[Feature]` (Task 3) used in Task 4. `flatfile_to_gff(rec) -> GffDocument` (Task 4) used in Task 5. `Span(seqid, start, end, strand, phase=…)` and `Feature(id, source, type, spans, attributes, parent_ids)` match `model.py` constructors. `MssConfig(source={}, transl_table=1, product_default=…)` + `cfg.emit_mrna=True` and `convert(doc, seqs, cfg, common_rows=[]) -> (MssDocument, diags)` + `emit_ann` match the gff2mss usage established in prior plans.

**Known caveat (documented, non-blocking):** partial-location markers (`<`/`>`) on the alt-spliced mRNA-B are not preserved through the GFF (the canonical model has no explicit partial field; gff2mss recomputes partiality from start/stop codons). The round-trip is therefore asserted on **CDS translation equality** (complete CDS), not on mRNA partial markers.
