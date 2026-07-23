# Heterosigma DDBJ 登録パイプライン Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `dev/heterosigma/` のヘテロシグマゲノムを、核 WGS と organelle(MT+CP)の 2 本の DDBJ MSS 登録ファイル(.ann/.fasta)に変換する。

**Architecture:** 4 ステップ(Step1 AGAT 標準化 → Step2 heterosigma 固有スクリプト → Step3 normalize → Step4 feature=本プロジェクト mss / COMMON・source=ddbj_mss_tools の共有 `common` を再利用する adapter)。共有コード(`src/ddbj_gff`)は feature 生成と document 組み立てを分離し、product ルール等を仕様へ是正する。

**Tech Stack:** Python 3.11+, BioPython, pytest(TDD)。AGAT(amd64 コンテナ)。ddbj_mss_tools の `common` パッケージ(pydantic)を実コード再利用。

## Global Constraints

- Python `>=3.11`。依存は BioPython のみ追加禁止(pydantic 等は `ddbj_mss_tools` 側の依存で、adapter 実行時のみ `sys.path` 経由)。
- TDD 厳守。テストは `-m "not slow"` が既定(`pyproject.toml`)。大規模実データは `@pytest.mark.slow`。
- 既存 mss テストを壊さない(`convert()`/`emit_ann()` の出力を後方互換で維持)。
- product 基本ルール: `(1)` id→product 表 → `(2)` col-9 `product` → `(3)` `"hypothetical protein"`。`"protein {gene名}"` フォールバックは**廃止**。
- 生物種: Heterosigma akashiwo, taxid `2829`。transl_table: 核=1 / MT=1 / CP=11(CP は CDS 属性優先で自動)。
- pseudogene(主に tRNAscan tRNA 擬遺伝子)は登録から**除外**。
- product クリーニングは**末尾 ` [EC:...]` の除去のみ**(`/EC_number` 化や `---` 正規化はしない)。
- CDS 内部 stop 検出時は CDS/mRNA を出さず `misc_feature`(translate/product なし、note に概要)。末端の開始/終止 codon 欠如のみは partial `<`,`>`。
- annotation_kaas.tsv は **CRLF**。キーは `anno1.` 除去で GFF transcript ID と一致。
- heterosigma 固有物は `dev/heterosigma/` に置く(公開レポジトリに追加不要)。
- ddbj_mss_tools は `../../../../ddbj_mss_tools`(env `DDBJ_MSS_TOOLS_SRC` で上書き可)。

---

## File Structure

**共有コード改修(`src/ddbj_gff/`)**
- `mss/config.py` — `MssConfig` に `product_map_path` / `product_map` 追加、`[product] map` パース(modify)
- `mss/product_map.py` — `load_product_map(path)->dict[str,str]`(create)
- `mss/convert.py` — `_product` 是正、`build_entry_features` 新設、RNA ヘルパ、parentless RNA、pseudogene skip、内部stop→misc_feature(modify)
- `mss/emit.py` — `feature_rows(feat)` 抽出(modify)
- `mss/cli.py` — gzip FASTA 対応 + product_map ロード(modify)
- `normalize/passes.py` — `pass_coerce_transcript_to_mrna`(modify)
- `normalize/normalize.py` / `normalize/config.py` — 新パス登録 + gate(modify)
- `normalize/cli.py` — gzip FASTA 対応(modify)

**heterosigma 固有(`dev/heterosigma/scripts/`)**
- `build_product_map.py`(create)
- `split_by_compartment.py`(create)
- `verify_agat.py`(create)
- `make_ann.py`(create, adapter)

**設定(`dev/heterosigma/`)**
- `nuclear.mss.toml` / `organelle.mss.toml` / `common_nuclear.json` / `common_organelle.json` / `sequence_roles.tsv`(create)

**テスト(`tests/`)**
- `test_mss_product.py` / `test_mss_parentless_rna.py` / `test_mss_pseudogene.py` / `test_mss_misc_feature.py` / `test_mss_entry_features.py`(create)
- `test_normalize_transcript_mrna.py`(create)
- 既存更新: `test_mss_cds.py`
- heterosigma スクリプト単体: `dev/heterosigma/scripts/tests/`(create)

---

## Task 1: product ルール是正 + product_map プラミング

**Files:**
- Create: `src/ddbj_gff/mss/product_map.py`
- Modify: `src/ddbj_gff/mss/config.py`, `src/ddbj_gff/mss/convert.py`
- Test: `tests/test_mss_product.py`（create）, `tests/test_mss_cds.py`（modify）

**Interfaces:**
- Produces: `load_product_map(path: str) -> dict[str, str]`; `MssConfig.product_map: dict[str,str]`（既定空）, `MssConfig.product_map_path: str | None`; `_product(mrna, gene, cfg) -> str`（ルール: map[mrna.id] → map[gene.id] → col-9 product → cfg.product_default）

- [ ] **Step 1: `test_mss_product.py` に失敗するテストを書く**

```python
from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import _product
from ddbj_gff.mss.product_map import load_product_map


def _mrna(mid, gid, product=None):
    attr = {"product": [product]} if product else {}
    mrna = Feature(mid, "S", "mRNA", [Span("c", 1, 9, "+")], attr, [])
    gene = Feature(gid, "S", "gene", [Span("c", 1, 9, "+")], {}, [])
    return mrna, gene


def test_product_map_hit_by_transcript_id():
    mrna, gene = _mrna("g1.t1", "g1")
    cfg = MssConfig(source={}, product_map={"g1.t1": "tubulin-tyrosine ligase"})
    assert _product(mrna, gene, cfg) == "tubulin-tyrosine ligase"


def test_product_map_hit_by_gene_id_when_transcript_misses():
    mrna, gene = _mrna("g1.t2", "g1")
    cfg = MssConfig(source={}, product_map={"g1": "some kinase"})
    assert _product(mrna, gene, cfg) == "some kinase"


def test_product_falls_back_to_col9():
    mrna, gene = _mrna("g1.t1", "g1", product="50S ribosomal protein L5")
    cfg = MssConfig(source={})
    assert _product(mrna, gene, cfg) == "50S ribosomal protein L5"


def test_product_defaults_to_hypothetical_no_gene_name_fallback():
    mrna, gene = _mrna("g1.t1", "g1")
    mrna.attributes["gene"] = ["MpX"]  # 旧仕様なら "protein MpX" だったが廃止
    cfg = MssConfig(source={}, product_default="hypothetical protein")
    assert _product(mrna, gene, cfg) == "hypothetical protein"


def test_load_product_map_reads_tsv(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("g1.t1\ttubulin\ng1\ttubulin\n", encoding="utf-8")
    m = load_product_map(str(p))
    assert m == {"g1.t1": "tubulin", "g1": "tubulin"}
```

- [ ] **Step 2: テストが失敗することを確認**

Run: `uv run pytest tests/test_mss_product.py -v`
Expected: FAIL（`load_product_map` 未定義 / `MssConfig` に `product_map` なし / `_product` が旧仕様）

- [ ] **Step 3: `product_map.py` を実装**

```python
from __future__ import annotations


def load_product_map(path: str) -> dict[str, str]:
    """Read a 2-column TSV (id<TAB>product) into a dict. Blank lines skipped."""
    result: dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            key, product = parts[0].strip(), parts[1].strip()
            if key and product:
                result[key] = product
    return result
```

- [ ] **Step 4: `MssConfig` に product_map を追加**

`src/ddbj_gff/mss/config.py` の `MssConfig` dataclass に追加:

```python
    product_default: str = "hypothetical protein"
    product_map_path: str | None = None
    product_map: dict = field(default_factory=dict)
    transcript_mode: str = "nonredundant"
```

`from dataclasses import dataclass, field` に `field` を含める。`load_config` の `MssConfig(...)` 構築に `product_map_path=product.get("map")` を追加。

- [ ] **Step 5: `_product` を是正**

`src/ddbj_gff/mss/convert.py` の `_product` を置換:

```python
def _product(mrna, gene, cfg: MssConfig) -> str:
    pmap = cfg.product_map or {}
    hit = pmap.get(mrna.id) or (pmap.get(gene.id) if gene and gene.id else None)
    if hit:
        return hit
    vals = mrna.attributes.get("product")
    if vals and vals[0]:
        return vals[0]
    return cfg.product_default
```

- [ ] **Step 6: 既存テストを新仕様へ更新**

`tests/test_mss_cds.py::test_product_protein_gene_default` を置換:

```python
def test_product_defaults_to_hypothetical_when_only_gene_name():
    genome = Seq("ATGAAATAA")
    gene, mrna = mrna_with_cds([(1, 9)], strand="+", phase0=0, mrna_attr={"gene": ["MpX"]})
    f = build_cds_feature(mrna, gene, "PFX_000010", genome, cfg(), [])
    q = {x.key: x.value for x in f.qualifiers}
    assert q["product"] == "hypothetical protein"  # "protein MpX" フォールバックは廃止
    assert q["gene"] == "MpX"
```

- [ ] **Step 7: テストが通ることを確認**

Run: `uv run pytest tests/test_mss_product.py tests/test_mss_cds.py tests/test_mss_config.py -v`
Expected: PASS

- [ ] **Step 8: Commit**

```bash
git add src/ddbj_gff/mss/product_map.py src/ddbj_gff/mss/config.py src/ddbj_gff/mss/convert.py tests/test_mss_product.py tests/test_mss_cds.py
git commit -m "feat(mss): product rule (map→col9→hypothetical), drop 'protein {gene}' fallback"
```

---

## Task 2: `emit.py` の feature_rows ヘルパ抽出

**Files:**
- Modify: `src/ddbj_gff/mss/emit.py`
- Test: `tests/test_mss_emit.py`（既存が通ること）

**Interfaces:**
- Produces: `feature_rows(feat: MssFeature) -> list[list[str]]`（col1 は常に空。行 = `["", key_or_blank, loc_or_blank, qual_key, qual_val]`）。adapter が再利用する。

- [ ] **Step 1: 失敗テストを書く**（`tests/test_mss_emit.py` に追記）

```python
from ddbj_gff.mss.emit import feature_rows


def test_feature_rows_first_row_carries_key_and_location():
    feat = MssFeature("CDS", "1..9", [MssQualifier("locus_tag", "L_1"),
                                       MssQualifier("product", "x")])
    rows = feature_rows(feat)
    assert rows[0] == ["", "CDS", "1..9", "locus_tag", "L_1"]
    assert rows[1] == ["", "", "", "product", "x"]
```

- [ ] **Step 2: 失敗確認**

Run: `uv run pytest tests/test_mss_emit.py::test_feature_rows_first_row_carries_key_and_location -v`
Expected: FAIL（`feature_rows` 未定義）

- [ ] **Step 3: `emit.py` を実装(抽出 + emit_ann 再利用)**

```python
from __future__ import annotations

from .model import MssDocument, MssFeature, MssQualifier


def feature_rows(feat: MssFeature) -> list[list[str]]:
    quals = feat.qualifiers or [MssQualifier("", "")]
    rows: list[list[str]] = []
    for i, q in enumerate(quals):
        col2 = feat.key if i == 0 else ""
        col3 = feat.location if i == 0 else ""
        rows.append(["", col2, col3, q.key, q.value])
    return rows


def emit_ann(doc: MssDocument) -> str:
    lines: list[str] = list(doc.common_rows)
    for entry in doc.entries:
        rows: list[list[str]] = []
        for feat in entry.features:
            rows.extend(feature_rows(feat))
        if rows:
            rows[0][0] = entry.name
        for r in rows:
            lines.append("\t".join(r))
    return "\n".join(lines) + "\n"
```

`emit_fasta` は変更しない。

- [ ] **Step 4: テストが通ることを確認(既存の emit テスト含む)**

Run: `uv run pytest tests/test_mss_emit.py tests/test_mss_snapshot.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/mss/emit.py tests/test_mss_emit.py
git commit -m "refactor(mss): extract feature_rows() from emit_ann for adapter reuse"
```

---

## Task 3: feature-only API（`build_entry_features`）+ convert 委譲 + 位置ソート

**Files:**
- Modify: `src/ddbj_gff/mss/convert.py`
- Test: `tests/test_mss_entry_features.py`（create）, 既存 `tests/test_mss_convert.py` が通ること

**Interfaces:**
- Produces: `build_entry_features(doc, seqs, cfg, diagnostics) -> dict[str, list[MssFeature]]`（source/gap を含まない。entry 内 feature は開始位置でソート。locus_tag は run 全体で連番）
- Consumes: `LocusTagAssigner`, `build_gene_features`（Task 1 の `_product`）

- [ ] **Step 1: 失敗テストを書く**

```python
from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import build_entry_features, convert

GFF = """##gff-version 3
c1\tS\tgene\t1\t9\t.\t+\t.\tID=g1
c1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1
c1\tS\texon\t1\t9\t.\t+\t.\tID=e1;Parent=g1.t1
c1\tS\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=g1.t1
"""


def cfg():
    return MssConfig(source={"organism": "x", "mol_type": "genomic DNA"},
                     locus_tag_prefix="PFX")


def test_build_entry_features_has_no_source():
    doc = parse(GFF)
    per = build_entry_features(doc, {"c1": Seq("ATGAAATAA")}, cfg(), [])
    assert set(per) == {"c1"}
    assert all(f.key != "source" for f in per["c1"])
    assert any(f.key == "CDS" for f in per["c1"])


def test_convert_still_emits_source_first():
    doc = parse(GFF)
    mss, _ = convert(doc, {"c1": Seq("ATGAAATAA")}, cfg(), ["COMMON"])
    assert mss.entries[0].features[0].key == "source"
```

- [ ] **Step 2: 失敗確認**

Run: `uv run pytest tests/test_mss_entry_features.py -v`
Expected: FAIL（`build_entry_features` 未定義）

- [ ] **Step 3: `build_entry_features` を実装し、`convert` を委譲に書き換える**

`convert.py` に追加 + 既存 `convert` を置換:

```python
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
        items = [(_span_start(g), g) for g in genes]
        items.sort(key=lambda t: t[0])
        feats: list = []
        for _, gene in items:
            feats.extend(build_gene_features(gene, cfg.transcript_mode, assigner,
                                             genome_seq, cfg, diagnostics))
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
```

（旧 `convert` 内のインラインループは削除。`_span_start` は既存を使用。）

- [ ] **Step 4: テストが通ることを確認(既存 convert / snapshot / integration 含む)**

Run: `uv run pytest tests/test_mss_entry_features.py tests/test_mss_convert.py tests/test_mss_snapshot.py tests/test_mss_integration.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/mss/convert.py tests/test_mss_entry_features.py
git commit -m "refactor(mss): add build_entry_features(); convert() delegates (source/gap split out)"
```

---

## Task 4: 親 gene のない ncRNA/tRNA/rRNA のトップレベル出力

**Files:**
- Modify: `src/ddbj_gff/mss/convert.py`
- Test: `tests/test_mss_parentless_rna.py`（create）, 既存 `tests/test_mss_noncoding.py` が通ること

**Interfaces:**
- Produces: `build_rna_feature(rna, locus_tag, seqlen, gene_id, tx_id) -> MssFeature`; `build_gene_features`/`build_noncoding_features` はこれを再利用。`build_entry_features` は親なし RNA(type ∈ ncRNA/tRNA/rRNA/miRNA/snRNA/snoRNA/tmRNA/pre_miRNA + misc)を単独 emit。

- [ ] **Step 1: 失敗テストを書く**

```python
from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import build_entry_features

GFF = """##gff-version 3
c1\tInfernal\tncRNA\t10\t50\t.\t+\t.\tID=n1;Name=U1;Dbxref=RFAM:RF00003;note=hit
c1\ttRNAscan-SE\ttRNA\t60\t90\t.\t-\t.\tID=t1;isotype=Thr;anticodon=CGT
c1\tpybarrnap\trRNA\t100\t200\t.\t+\t.\tID=r1;product=18S ribosomal RNA
"""


def cfg():
    return MssConfig(source={}, locus_tag_prefix="PFX")


def test_parentless_rna_emitted_as_toplevel_with_locus_tag():
    doc = parse(GFF)
    per = build_entry_features(doc, {"c1": Seq("A" * 300)}, cfg(), [])
    keys = [f.key for f in per["c1"]]
    assert keys == ["ncRNA", "tRNA", "rRNA"]  # 位置順
    nc = {q.key: q.value for q in per["c1"][0].qualifiers}
    assert nc["locus_tag"] == "PFX_000010"
    assert nc["ncRNA_class"] == "other"  # RFAM 由来の未知クラスは other
    assert any(q.key == "db_xref" and q.value == "RFAM:RF00003" for q in per["c1"][0].qualifiers)
    rr = {q.key: q.value for q in per["c1"][2].qualifiers}
    assert rr["product"] == "18S ribosomal RNA"
    tr = {q.key: q.value for q in per["c1"][1].qualifiers}
    assert tr["product"] == "tRNA-Thr"  # isotype から導出
```

- [ ] **Step 2: 失敗確認**

Run: `uv run pytest tests/test_mss_parentless_rna.py -v`
Expected: FAIL

- [ ] **Step 3: `build_rna_feature` を実装し、既存 noncoding と parentless で共有**

`convert.py` に追加。`_RNA_MAP`/`_STRUCTURAL` は既存を使用。`_ncRNA_KNOWN_CLASSES` を定義（既知クラスはそのまま、未知は `other`）。

```python
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
    quals.append(_submitter_note_ids(gene_id, tx_id))
    return MssFeature(feat_key, location, quals)
```

`build_noncoding_features` を `build_rna_feature` 使用へ書き換え（挙動を維持）:

```python
def build_noncoding_features(gene, locus_tag: str, seqlen: int, cfg) -> list:
    features = []
    for rna in gene.children:
        if rna.type in _STRUCTURAL:
            continue
        features.append(build_rna_feature(rna, locus_tag, seqlen, gene.id, rna.id,
                                           gene.gene or rna.gene))
    return features
```

`build_entry_features` の `items` に親なし RNA を追加:

```python
        parentless = [f for f in doc.roots
                      if f.type in _PARENTLESS_RNA_TYPES
                      and any(s.seqid == seqid for s in f.spans)]
        items = [(_span_start(g), g) for g in genes] + [(_span_start(r), r) for r in parentless]
        items.sort(key=lambda t: t[0])
        feats = []
        for _, feat in items:
            if feat.type == "gene":
                feats.extend(build_gene_features(feat, cfg.transcript_mode, assigner,
                                                 genome_seq, cfg, diagnostics))
            else:
                feats.append(build_rna_feature(feat, assigner.assign(feat),
                                                len(genome_seq), feat.id, feat.id, feat.gene))
```

- [ ] **Step 4: テストが通ることを確認(既存 noncoding 含む)**

Run: `uv run pytest tests/test_mss_parentless_rna.py tests/test_mss_noncoding.py tests/test_mss_convert.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/mss/convert.py tests/test_mss_parentless_rna.py
git commit -m "feat(mss): emit parentless ncRNA/tRNA/rRNA as top-level features (shared build_rna_feature)"
```

---

## Task 5: pseudogene の除外(変換器側の二重防御)

**Files:**
- Modify: `src/ddbj_gff/mss/convert.py`
- Test: `tests/test_mss_pseudogene.py`（create）

**Interfaces:**
- Consumes: `build_entry_features`
- Behavior: `type == "pseudogene"`、または gene/root で属性 `gene_biotype == "pseudogene"` の feature を emit しない。除外件数を `Diagnostic("pseudogene-skipped")` で報告。

- [ ] **Step 1: 失敗テストを書く**

```python
from Bio.Seq import Seq
from ddbj_gff import parse
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import build_entry_features

GFF = """##gff-version 3
c1\ttRNAscan-SE\tpseudogene\t10\t80\t.\t-\t.\tID=p1;gene_biotype=pseudogene
c1\ttRNAscan-SE\ttRNA\t100\t170\t.\t+\t.\tID=t1;isotype=Ala
"""


def test_pseudogene_excluded():
    doc = parse(GFF)
    diags = []
    per = build_entry_features(doc, {"c1": Seq("A" * 300)}, MssConfig(source={}, locus_tag_prefix="P"), diags)
    keys = [f.key for f in per["c1"]]
    assert "pseudogene" not in keys
    assert keys == ["tRNA"]
    assert any(d.code == "pseudogene-skipped" for d in diags)
```

- [ ] **Step 2: 失敗確認**

Run: `uv run pytest tests/test_mss_pseudogene.py -v`
Expected: FAIL（pseudogene が misc_RNA として出る or カウントされない）

- [ ] **Step 3: `build_entry_features` に除外を実装**

`items` 構築前にフィルタを追加:

```python
        def _is_pseudogene(f):
            return f.type == "pseudogene" or f._first("gene_biotype") == "pseudogene"

        skipped = 0
        genes2, parentless2 = [], []
        for g in genes:
            if _is_pseudogene(g):
                skipped += 1
            else:
                genes2.append(g)
        for r in parentless:
            if _is_pseudogene(r):
                skipped += 1
            else:
                parentless2.append(r)
        if skipped:
            diagnostics.append(Diagnostic(Severity.WARNING, None, "pseudogene-skipped",
                                          f"{seqid}: skipped {skipped} pseudogene feature(s)"))
        genes, parentless = genes2, parentless2
```

（`parentless` を先に定義しておくこと。Task 4 の `parentless` 定義行の直後にこのブロックを置き、`items` はフィルタ後の `genes`/`parentless` から作る。）

- [ ] **Step 4: テストが通ることを確認**

Run: `uv run pytest tests/test_mss_pseudogene.py tests/test_mss_parentless_rna.py tests/test_mss_convert.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/ddbj_gff/mss/convert.py tests/test_mss_pseudogene.py
git commit -m "feat(mss): exclude pseudogene features from MSS output (with report)"
```

---

## Task 6: CDS 内部 stop → misc_feature

**Files:**
- Modify: `src/ddbj_gff/mss/convert.py`
- Test: `tests/test_mss_misc_feature.py`（create）, 既存 `tests/test_mss_cds.py` の内部stopテストを更新

**Interfaces:**
- Behavior: `build_cds_feature` は body に内部 `*` があるとき、CDS ではなく `misc_feature`（location=CDS spans、qual=locus_tag + note、product/translate なし）を返す。`build_gene_features` は返りが `misc_feature` のとき mRNA を出さずそれのみ append。

- [ ] **Step 1: 失敗テスト（`tests/test_mss_misc_feature.py`）**

```python
from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import build_cds_feature, build_gene_features
from ddbj_gff.mss.locus_tag import LocusTagAssigner


def cfg():
    return MssConfig(source={}, transl_table=1, product_default="hypothetical protein")


def _gene_with_internal_stop():
    genome = Seq("ATGTAAAAA")  # ATG TAA AAA -> M * K  内部stop
    gene = Feature("g", "S", "gene", [Span("c", 1, 9, "+")], {"ID": ["g"]}, [])
    mrna = Feature("g.t1", "S", "mRNA", [Span("c", 1, 9, "+")], {"ID": ["g.t1"]}, ["g"])
    cds = Feature("cds", "S", "CDS", [Span("c", 1, 9, "+", 0)], {"ID": ["cds"]}, ["g.t1"])
    mrna.children = [cds]
    gene.children = [mrna]
    return gene, mrna, genome


def test_internal_stop_returns_misc_feature():
    gene, mrna, genome = _gene_with_internal_stop()
    f = build_cds_feature(mrna, gene, "L_1", genome, cfg(), [])
    assert f.key == "misc_feature"
    q = {x.key: x.value for x in f.qualifiers}
    assert "product" not in q
    assert q["locus_tag"] == "L_1"
    assert any(x.key == "note" and "internal stop" in x.value.lower() for x in f.qualifiers)


def test_gene_with_internal_stop_emits_only_misc_feature():
    gene, mrna, genome = _gene_with_internal_stop()
    assigner = LocusTagAssigner("L", 6, 10, 10)
    feats = build_gene_features(gene, "nonredundant", assigner, genome, cfg(), [])
    keys = [f.key for f in feats]
    assert keys == ["misc_feature"]  # mRNA/CDS は出さない
```

- [ ] **Step 2: 失敗確認**

Run: `uv run pytest tests/test_mss_misc_feature.py -v`
Expected: FAIL

- [ ] **Step 3: `build_cds_feature` を修正（内部stop→misc_feature）**

`build_cds_feature` 内、`body` 計算後の内部stopチェックを置換。現状:

```python
    body = protein[:-1] if protein.endswith("*") else protein
    if "*" in body:
        diagnostics.append(Diagnostic(Severity.WARNING, None, "translation-internal-stop",
                                      f"CDS {mrna.id!r} has an internal stop codon"))
```

を次に変更（misc_feature を返す）:

```python
    body = protein[:-1] if protein.endswith("*") else protein
    if "*" in body:
        diagnostics.append(Diagnostic(Severity.WARNING, None, "translation-internal-stop",
                                      f"CDS {mrna.id!r} has an internal stop codon"))
        loc = build_insdc_location(spans, len(genome_seq))
        note = (f"internal stop codon(s) detected in CDS {mrna.id}; "
                f"not translated")
        quals = [MssQualifier("locus_tag", locus_tag), MssQualifier("note", note)]
        return MssFeature("misc_feature", loc, quals)
```

（`build_insdc_location` は import 済み。以降の CDS 組み立てはそのまま。）

- [ ] **Step 4: `build_gene_features` を修正（misc_feature 時は mRNA を出さない）**

`build_gene_features` の per-transcript ループを、CDS を先に計算して分岐する形へ:

```python
    for mrna in transcripts:
        if not collect_spans(mrna, "exon") and not collect_spans(mrna, "CDS"):
            diagnostics.append(Diagnostic(Severity.WARNING, None, "no-exon",
                                          f"mRNA {mrna.id!r} has no exon or CDS; skipped"))
            continue
        if locus_tag is None:
            locus_tag = assigner.assign(gene)
        cds = build_cds_feature(mrna, gene, locus_tag, genome_seq, cfg, diagnostics)
        if cds is not None and cds.key == "misc_feature":
            features.append(cds)
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
```

- [ ] **Step 5: 既存の内部stopテストを更新**

`tests/test_mss_cds.py::test_internal_stop_diagnostic` を新挙動へ:

```python
def test_internal_stop_returns_misc_feature():
    genome = Seq("ATGTAAAAA")  # M*K
    gene, mrna = mrna_with_cds([(1, 9)], strand="+", phase0=0)
    diags = []
    f = build_cds_feature(mrna, gene, "PFX_000010", genome, cfg(), diags)
    assert f.key == "misc_feature"
    assert any(d.code == "translation-internal-stop" for d in diags)
```

- [ ] **Step 6: テストが通ることを確認（recoded_codon が CDS のままなことも）**

Run: `uv run pytest tests/test_mss_misc_feature.py tests/test_mss_cds.py tests/test_mss_gene_features.py -v`
Expected: PASS（`test_recoded_codon_child_avoids_internal_stop_warning` は CDS を返す）

- [ ] **Step 7: Commit**

```bash
git add src/ddbj_gff/mss/convert.py tests/test_mss_misc_feature.py tests/test_mss_cds.py
git commit -m "feat(mss): emit misc_feature for CDS with internal stop (no translate/product)"
```

---

## Task 7: gzip FASTA 対応 + normalize transcript→mRNA パス + product_map ロード

**Files:**
- Create: `src/ddbj_gff/io.py`
- Modify: `src/ddbj_gff/mss/cli.py`, `src/ddbj_gff/normalize/cli.py`, `src/ddbj_gff/normalize/passes.py`, `src/ddbj_gff/normalize/normalize.py`, `src/ddbj_gff/normalize/config.py`
- Test: `tests/test_normalize_transcript_mrna.py`（create）, `tests/test_io.py`（create）

**Interfaces:**
- Produces: `open_text(path)`（`.gz` は gzip 展開して text mode で開く context manager）; `pass_coerce_transcript_to_mrna(doc, ctx) -> list[Change]`（CDS を子に持つ `transcript` を `mRNA` へ）; `NormalizeConfig.coerce_transcript_to_mrna: bool = True`

- [ ] **Step 1: 失敗テストを書く**

`tests/test_io.py`:

```python
import gzip
from ddbj_gff.io import open_text


def test_open_text_reads_gzip(tmp_path):
    p = tmp_path / "x.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">a\nACGT\n")
    with open_text(str(p)) as fh:
        assert fh.read() == ">a\nACGT\n"
```

`tests/test_normalize_transcript_mrna.py`:

```python
from ddbj_gff import parse
from ddbj_gff.normalize.passes import pass_coerce_transcript_to_mrna, NormalizeContext
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.validate.vocab import Vocab

GFF = """##gff-version 3
c1\tS\tgene\t1\t9\t.\t+\t.\tID=g1
c1\tS\ttranscript\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1
c1\tS\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=g1.t1
c1\tS\ttranscript\t20\t30\t.\t+\t.\tID=nc1;Parent=g2
"""


def test_coding_transcript_becomes_mrna_noncoding_untouched():
    doc = parse(GFF)
    ctx = NormalizeContext(vocab=Vocab.load(), seq_lengths=None, config=NormalizeConfig())
    pass_coerce_transcript_to_mrna(doc, ctx)
    types = {f.id: f.type for f in doc.features}
    assert types["g1.t1"] == "mRNA"   # CDS を持つ transcript
    assert types["nc1"] == "transcript"  # CDS なしは不変
```

- [ ] **Step 2: 失敗確認**

Run: `uv run pytest tests/test_io.py tests/test_normalize_transcript_mrna.py -v`
Expected: FAIL

- [ ] **Step 3: `io.py` を実装**

```python
from __future__ import annotations

import gzip
from contextlib import contextmanager


@contextmanager
def open_text(path: str, encoding: str = "utf-8", errors: str = "strict"):
    if path.endswith(".gz"):
        fh = gzip.open(path, "rt", encoding=encoding, errors=errors)
    else:
        fh = open(path, "r", encoding=encoding, errors=errors)
    try:
        yield fh
    finally:
        fh.close()
```

- [ ] **Step 4: normalize パスを実装**

`normalize/passes.py` に追加:

```python
def pass_coerce_transcript_to_mrna(doc, ctx) -> list:
    changes: list = []
    if not getattr(ctx.config, "coerce_transcript_to_mrna", True):
        return changes
    for f in doc.features:
        if f.type != "transcript":
            continue
        if any(c.type == "CDS" for c in f.children):
            f.type = "mRNA"
            changes.append(Change("rename-type", f.id or "?", "transcript -> mRNA (has CDS)"))
    return changes
```

`normalize/config.py` の `NormalizeConfig` に `coerce_transcript_to_mrna: bool = True` を追加し、`load_normalize_config` で `n.get("coerce_transcript_to_mrna", True)` を読む。`normalize/normalize.py` のパス実行リストに `pass_coerce_transcript_to_mrna` を（`pass_so_terms` の前に）登録する。

- [ ] **Step 5: cli の FASTA を gzip 対応にする**

`mss/cli.py` と `normalize/cli.py` の FASTA 読み込みを `SeqIO.parse` へ渡す前に `open_text` 経由にする。例（`mss/cli.py`）:

```python
from ..io import open_text
...
    with open_text(args.fasta) as fh:
        seqs = {rec.id: rec.seq for rec in SeqIO.parse(fh, "fasta")}
```

`mss/cli.py` に product_map ロードも追加:

```python
from .product_map import load_product_map
...
    if cfg.product_map_path:
        cfg.product_map = load_product_map(cfg.product_map_path)
```

- [ ] **Step 6: テストが通ることを確認**

Run: `uv run pytest tests/test_io.py tests/test_normalize_transcript_mrna.py tests/test_mss_cli.py -v`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add src/ddbj_gff/io.py src/ddbj_gff/mss/cli.py src/ddbj_gff/normalize/ tests/test_io.py tests/test_normalize_transcript_mrna.py
git commit -m "feat: gzip FASTA support, coding transcript->mRNA normalize pass, product_map loading"
```

---

## Task 8: `build_product_map.py`（heterosigma 固有）

**Files:**
- Create: `dev/heterosigma/scripts/build_product_map.py`
- Test: `dev/heterosigma/scripts/tests/test_build_product_map.py`

**Interfaces:**
- CLI: `python build_product_map.py --kaas annotation_kaas.tsv --out product_map.tsv`
- Produces: `build_map(rows: list[str]) -> dict[str,str]`（CRLF 済みの各行文字列を受け、`anno1.` 除去・末尾 ` [EC:...]` 除去・空スキップ・transcript+gene 両キー）

- [ ] **Step 1: 失敗テストを書く**

```python
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from build_product_map import build_map, clean_product


def test_clean_product_strips_trailing_ec():
    assert clean_product("tubulin---tyrosine ligase [EC:6.3.2.25]") == "tubulin---tyrosine ligase"
    assert clean_product("serine/threonine protein kinase SCH9 [EC:2.7.11.1]") == "serine/threonine protein kinase SCH9"
    assert clean_product("acyl carrier protein") == "acyl carrier protein"


def test_build_map_strips_anno1_and_adds_gene_key():
    rows = ["protein_acc\tDescription_2",
            "anno1.g7.t1\ttubulin [EC:6.3.2.25]",
            "anno1.g8.t1\t",  # 空 -> スキップ
            "anno1.both_agree_g273.t2\tsome kinase"]
    m = build_map(rows)
    assert m["g7.t1"] == "tubulin"
    assert m["g7"] == "tubulin"                     # gene キーも
    assert "g8.t1" not in m and "g8" not in m       # 空はスキップ
    assert m["both_agree_g273.t2"] == "some kinase"
    assert m["both_agree_g273"] == "some kinase"
```

- [ ] **Step 2: 失敗確認**

Run: `cd dev/heterosigma/scripts && uv run pytest tests/test_build_product_map.py -v`
Expected: FAIL

- [ ] **Step 3: 実装**

```python
#!/usr/bin/env python3
"""Build id->product TSV from KAAS annotation (CRLF, anno1. prefix, trailing [EC:...])."""
from __future__ import annotations

import argparse
import re

_EC_RE = re.compile(r"\s*\[EC:[^\]]*\]\s*$")


def clean_product(desc: str) -> str:
    return _EC_RE.sub("", desc).strip()


def build_map(rows: list[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for i, raw in enumerate(rows):
        line = raw.rstrip("\r\n")
        if i == 0 and line.lower().startswith("protein_acc"):
            continue
        if not line:
            continue
        parts = line.split("\t")
        key = parts[0].strip()
        desc = parts[1].strip() if len(parts) > 1 else ""
        product = clean_product(desc)
        if not key or not product:
            continue
        tx = key[len("anno1."):] if key.startswith("anno1.") else key
        gene = re.sub(r"\.t\d+$", "", tx)
        result[tx] = product
        result[gene] = product
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kaas", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    with open(args.kaas, encoding="utf-8", errors="replace") as fh:
        rows = fh.readlines()
    m = build_map(rows)
    with open(args.out, "w", encoding="utf-8") as fh:
        for k in sorted(m):
            fh.write(f"{k}\t{m[k]}\n")
    print(f"[build_product_map] {len(m)} keys -> {args.out}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: テストが通ることを確認**

Run: `cd dev/heterosigma/scripts && uv run pytest tests/test_build_product_map.py -v`
Expected: PASS

- [ ] **Step 5: 実データで生成して健全性確認**

Run:
```bash
cd dev/heterosigma
uv run python scripts/build_product_map.py --kaas annotation_kaas.tsv --out product_map.tsv
wc -l product_map.tsv   # 期待: 4446(transcript) + gene キー分（空は除外済み）
```
Expected: 空 description が除外され、transcript/gene 両キーが出力される

- [ ] **Step 6: Commit**（テストとスクリプトのみ。product_map.tsv/実データはコミットしない）

```bash
git add dev/heterosigma/scripts/build_product_map.py dev/heterosigma/scripts/tests/test_build_product_map.py
git commit -m "feat(heterosigma): build_product_map.py (CRLF, anno1., EC strip, tx+gene keys)"
```

---

## Task 9: `split_by_compartment.py`（分割 + pseudogene 除外）

**Files:**
- Create: `dev/heterosigma/scripts/split_by_compartment.py`
- Test: `dev/heterosigma/scripts/tests/test_split_by_compartment.py`

**Interfaces:**
- CLI: `python split_by_compartment.py --gff standardized.agat.gff3 --nuclear nuclear.gff3 --organelle organelle.gff3 [--organelle-seqids MT,CP] [--drop-pseudogene]`
- Produces: `classify(seqid, organelle_seqids) -> "organelle"|"nuclear"`; `split_lines(lines, organelle_seqids, drop_pseudogene) -> tuple[list,list]`（Parent 鎖ごと分割。pseudogene とその子孫を除外）

- [ ] **Step 1: 失敗テストを書く**

```python
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from split_by_compartment import split_lines

LINES = [
    "##gff-version 3\n",
    "scaffold_1\tS\tgene\t1\t9\t.\t+\t.\tID=g1\n",
    "scaffold_1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1\n",
    "scaffold_1\ttRNAscan\tpseudogene\t20\t30\t.\t+\t.\tID=p1;gene_biotype=pseudogene\n",
    "MT\tLiftoff\tgene\t1\t9\t.\t+\t.\tID=gm1\n",
    "MT\tLiftoff\tCDS\t1\t9\t.\t+\t.\tID=cm1;Parent=gm1\n",
]


def test_split_and_drop_pseudogene():
    nuc, org = split_lines(LINES, {"MT", "CP"}, drop_pseudogene=True)
    assert any("ID=g1;" in l or "ID=g1\n" in l for l in nuc)
    assert not any("pseudogene" in l for l in nuc)   # p1 除外
    assert any("ID=gm1" in l for l in org)
    assert all(not l.startswith("scaffold_1") for l in org)
    assert all(not l.startswith("MT") for l in nuc)
```

- [ ] **Step 2: 失敗確認**

Run: `cd dev/heterosigma/scripts && uv run pytest tests/test_split_by_compartment.py -v`
Expected: FAIL

- [ ] **Step 3: 実装**

```python
#!/usr/bin/env python3
"""Split a standardized GFF3 into nuclear vs organelle by seqid; drop pseudogenes."""
from __future__ import annotations

import argparse
import re


def _attr(col9: str, key: str) -> str | None:
    m = re.search(rf"(?:^|;){key}=([^;]+)", col9)
    return m.group(1) if m else None


def split_lines(lines: list[str], organelle_seqids: set[str], drop_pseudogene: bool):
    # 1st pass: collect IDs of pseudogene features to drop (with descendants).
    drop_ids: set[str] = set()
    if drop_pseudogene:
        for line in lines:
            if line.startswith("#") or "\t" not in line:
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            if cols[2] == "pseudogene" or _attr(cols[8], "gene_biotype") == "pseudogene":
                fid = _attr(cols[8], "ID")
                if fid:
                    drop_ids.add(fid)
        # propagate to descendants (Parent chain), iterate to fixpoint
        changed = True
        while changed:
            changed = False
            for line in lines:
                if line.startswith("#") or "\t" not in line:
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                parent = _attr(cols[8], "Parent")
                fid = _attr(cols[8], "ID")
                if parent and any(p in drop_ids for p in parent.split(",")) and fid and fid not in drop_ids:
                    drop_ids.add(fid)
                    changed = True

    header = [l for l in lines if l.startswith("#")]
    nuclear = list(header)
    organelle = list(header)
    for line in lines:
        if line.startswith("#") or "\t" not in line:
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        fid = _attr(cols[8], "ID")
        parent = _attr(cols[8], "Parent")
        if drop_pseudogene and (fid in drop_ids or (parent and any(p in drop_ids for p in parent.split(",")))):
            continue
        target = organelle if cols[0] in organelle_seqids else nuclear
        target.append(line)
    return nuclear, organelle


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--nuclear", required=True)
    ap.add_argument("--organelle", required=True)
    ap.add_argument("--organelle-seqids", default="MT,CP")
    ap.add_argument("--drop-pseudogene", action="store_true")
    args = ap.parse_args()
    org_ids = {s.strip() for s in args.organelle_seqids.split(",") if s.strip()}
    with open(args.gff, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    nuc, org = split_lines(lines, org_ids, args.drop_pseudogene)
    with open(args.nuclear, "w", encoding="utf-8") as fh:
        fh.writelines(nuc)
    with open(args.organelle, "w", encoding="utf-8") as fh:
        fh.writelines(org)
    print(f"[split] nuclear={len(nuc)} organelle={len(org)} lines")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: テストが通ることを確認**

Run: `cd dev/heterosigma/scripts && uv run pytest tests/test_split_by_compartment.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dev/heterosigma/scripts/split_by_compartment.py dev/heterosigma/scripts/tests/test_split_by_compartment.py
git commit -m "feat(heterosigma): split_by_compartment.py (nuclear/organelle split + pseudogene drop)"
```

---

## Task 10: `make_ann.py` adapter（feature=mss / COMMON・source=ddbj_mss_tools.common）

**Files:**
- Create: `dev/heterosigma/scripts/make_ann.py`
- Test: `dev/heterosigma/scripts/tests/test_make_ann.py`

**Interfaces:**
- CLI: `python make_ann.py --gff G.gff3 --fasta F.fa[.gz] --mss-config C.toml --common common.json --sequence-roles roles.tsv --submission-category WGS|GNM --out OUT`（`OUT.ann` と `OUT.fasta` を出力）
- Consumes: `ddbj_gff.mss`（`build_entry_features`, `feature_rows`, `load_config`, `load_product_map`）, `ddbj_gff.io.open_text`, `ddbj_gff.parse`; `common.models.load_common_json`, `common.common_builder.create_common`, `common.source_builder.load_sequence_roles/source_qualifier/ff_definition`, `common.gap_annotator.GapAnnotator/annotate_gaps`, `common.submission_category.get_category_rules`
- Behavior: egapx2mss の `write_ddbj_ann` を踏襲。COMMON → 各 entry（circular なら `TOPOLOGY` 行 → source 行 → feature 行 → gap 行）。`ddbj_mss_tools/src` を `sys.path` に追加（env `DDBJ_MSS_TOOLS_SRC`、既定は 4 階層上の `ddbj_mss_tools/src`）。

- [ ] **Step 1: 失敗テストを書く**（tiny fixture、organelle circular 経路）

```python
import sys, os, importlib.util
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

# ddbj_mss_tools が import できない環境では skip
import pytest
_MSS_TOOLS = os.environ.get("DDBJ_MSS_TOOLS_SRC",
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../ddbj_mss_tools/src")))
if not os.path.isdir(os.path.join(_MSS_TOOLS, "common")):
    pytest.skip("ddbj_mss_tools not available", allow_module_level=True)

from make_ann import build_ann_text


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return str(p)


def test_organelle_circular_source_and_topology(tmp_path):
    gff = _write(tmp_path, "o.gff3",
        "##gff-version 3\n"
        "CP\tLiftoff\tgene\t1\t9\t.\t+\t.\tID=g1;gene=rpl5\n"
        "CP\tLiftoff\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        "CP\tLiftoff\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g1.t1;product=50S ribosomal protein L5;transl_table=11\n")
    fasta = _write(tmp_path, "o.fa", ">CP\nATGAAATAA\n")
    mss_cfg = _write(tmp_path, "o.toml",
        '[source]\norganism="x"\nmol_type="genomic DNA"\n[locus_tag]\nprefix="HAKA"\n[cds]\ntransl_table=1\n')
    common = _write(tmp_path, "c.json",
        '{"SOURCE":{"organism":"Heterosigma akashiwo","mol_type":"genomic DNA"},"SOURCE_IDENTIFIER":"strain"}')
    roles = _write(tmp_path, "r.tsv", "CP\torganelle\tplastid:chloroplast\tcomplete\tcircular\n")
    text = build_ann_text(gff, fasta, mss_cfg, common, roles, "GNM")
    assert "TOPOLOGY" in text and "circular" in text
    assert "\torganelle\tplastid:chloroplast" in text
    assert "50S ribosomal protein L5" in text
    assert "COMMON" in text
```

- [ ] **Step 2: 失敗確認**

Run: `cd dev/heterosigma/scripts && uv run pytest tests/test_make_ann.py -v`
Expected: FAIL（`make_ann` 未実装。ddbj_mss_tools が無ければ skip）

- [ ] **Step 3: 実装**

```python
#!/usr/bin/env python3
"""Assemble a DDBJ MSS .ann: features from ddbj_gff.mss, COMMON/source from ddbj_mss_tools.common."""
from __future__ import annotations

import argparse
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_MSS_TOOLS = os.environ.get(
    "DDBJ_MSS_TOOLS_SRC", os.path.abspath(os.path.join(_HERE, "../../../../ddbj_mss_tools/src")))
sys.path.insert(0, _MSS_TOOLS)
# gff_submission の src も import 可能に
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "../../../src")))

from Bio import SeqIO

from ddbj_gff import parse
from ddbj_gff.io import open_text
from ddbj_gff.mss.config import load_config
from ddbj_gff.mss.convert import build_entry_features
from ddbj_gff.mss.emit import feature_rows, emit_fasta
from ddbj_gff.mss.product_map import load_product_map

from common.models import load_common_json
from common.common_builder import create_common
from common.source_builder import load_sequence_roles, source_qualifier, ff_definition
from common.gap_annotator import GapAnnotator, annotate_gaps
from common.submission_category import get_category_rules


def build_ann_text(gff_path, fasta_path, mss_config_path, common_path,
                   sequence_roles_path, submission_category) -> tuple[str, dict]:
    cfg, _ = load_config(mss_config_path)
    if cfg.product_map_path:
        cfg.product_map = load_product_map(cfg.product_map_path)

    with open_text(fasta_path) as fh:
        seqs = {rec.id: rec.seq for rec in SeqIO.parse(fh, "fasta")}
    with open_text(gff_path) as fh:
        doc = parse(fh.read())

    diagnostics: list = []
    per_entry = build_entry_features(doc, seqs, cfg, diagnostics)

    common = load_common_json(common_path)
    common_dict = common.model_dump(exclude_none=True)
    if submission_category:
        common_dict["_submission_category"] = submission_category
    roles = load_sequence_roles(sequence_roles_path) if sequence_roles_path else {}

    base_source = dict(common.SOURCE or {})
    organism = base_source.get("organism", "")
    src_id_key = common.SOURCE_IDENTIFIER
    infra = base_source.get(src_id_key, "") if src_id_key else ""
    mol_type = base_source.get("mol_type", "")

    all_ids = list(per_entry.keys())
    is_wgs = all((roles.get(e) is None or roles.get(e).type == "unplaced") for e in all_ids)

    # gap annotators from common.ASSEMBLY_GAP
    gap_annotators: list[GapAnnotator] = []
    gap_cfg = common.ASSEMBLY_GAP
    if gap_cfg:
        cfgs = gap_cfg if isinstance(gap_cfg, list) else [gap_cfg]
        gap_annotators = [GapAnnotator(linkage_evidence=c.linkage_evidence,
                                       min_gap_length=c.min_gap_length,
                                       max_gap_length=c.max_gap_length,
                                       gap_type=c.gap_type,
                                       estimated_length=c.estimated_length)
                          for c in cfgs if c.enabled]

    rows: list[list[str]] = list(create_common(common_dict))
    for entry_id in all_ids:
        seq = seqs[entry_id]
        length = len(seq)
        role = roles.get(entry_id)
        is_circular = role.is_circular if role is not None else False
        if is_circular:
            rows.append([entry_id, "TOPOLOGY", "", "circular", ""])
        src = dict(base_source)
        src.update(source_qualifier(role, entry_id, is_wgs))
        src["ff_definition"] = ff_definition(role, entry_id, organism, infra, mol_type, is_wgs)
        items = list(src.items())
        first_col = "" if is_circular else entry_id
        rows.append([first_col, "source", f"1..{length}", items[0][0], items[0][1]])
        for k, v in items[1:]:
            rows.append(["", "", "", k, str(v)])
        for feat in per_entry[entry_id]:
            rows.extend(feature_rows(feat))
        if gap_annotators:
            rows.extend(annotate_gaps(gap_annotators, str(seq)))

    ann_text = "\n".join("\t".join(r) for r in rows) + "\n"
    return ann_text, seqs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--mss-config", required=True)
    ap.add_argument("--common", required=True)
    ap.add_argument("--sequence-roles")
    ap.add_argument("--submission-category", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    ann_text, seqs = build_ann_text(args.gff, args.fasta, args.mss_config,
                                    args.common, args.sequence_roles, args.submission_category)
    with open(f"{args.out}.ann", "w", encoding="utf-8") as fh:
        fh.write(ann_text)
    with open(f"{args.out}.fasta", "w", encoding="utf-8") as fh:
        fh.write(emit_fasta(seqs))
    print(f"[make_ann] -> {args.out}.ann / {args.out}.fasta")


if __name__ == "__main__":
    main()
```

（`build_ann_text` はテストで `text` を返すため、戻り値の 1 つ目 `ann_text` を検証する。テストは `build_ann_text(...)[0]` を使うよう Step1 の `text = build_ann_text(...)` を `text, _ = build_ann_text(...)` に修正すること。）

- [ ] **Step 4: Step1 のテスト呼び出しを戻り値タプルに合わせて修正し、テストを通す**

Step1 の `text = build_ann_text(...)` → `text, _ = build_ann_text(...)`。

Run: `cd dev/heterosigma/scripts && uv run pytest tests/test_make_ann.py -v`
Expected: PASS（ddbj_mss_tools 有り環境）

- [ ] **Step 5: Commit**

```bash
git add dev/heterosigma/scripts/make_ann.py dev/heterosigma/scripts/tests/test_make_ann.py
git commit -m "feat(heterosigma): make_ann.py adapter (mss features + ddbj_mss_tools common COMMON/source)"
```

---

## Task 11: 設定ファイル（mss toml / common json / sequence_roles）

**Files:**
- Create: `dev/heterosigma/nuclear.mss.toml`, `dev/heterosigma/organelle.mss.toml`, `dev/heterosigma/common_nuclear.json`, `dev/heterosigma/common_organelle.json`, `dev/heterosigma/sequence_roles.tsv`

**Interfaces:**
- Consumes: Task 10 の `make_ann.py`、Task 1 の `MssConfig`（`[product] map`）

- [ ] **Step 1: `nuclear.mss.toml` を作成**

```toml
[source]
organism = "Heterosigma akashiwo"
mol_type = "genomic DNA"

[locus_tag]
prefix = "HAKA"   # 仮。BioSample 登録後の正式 prefix に差し替え
width = 6
start = 10
step = 10

[cds]
transl_table = 1

[product]
default = "hypothetical protein"
map = "product_map.tsv"
```

- [ ] **Step 2: `organelle.mss.toml` を作成**

```toml
[source]
organism = "Heterosigma akashiwo"
mol_type = "genomic DNA"

[locus_tag]
prefix = "HAKA"
width = 6
start = 10
step = 10

[cds]
transl_table = 1   # MT=1。CP は CDS の transl_table=11 属性が優先される

[product]
default = "hypothetical protein"
```

- [ ] **Step 3: `common_nuclear.json` を作成（プレースホルダ）**

```json
{
  "DATATYPE": {"type": "WGS"},
  "KEYWORD": {"keyword": ["WGS", "STANDARD_DRAFT"]},
  "DBLINK": {"project": "PRJDBxxxxxx", "biosample": "SAMDxxxxxxxx"},
  "SUBMITTER": {"ab_name": ["Author,A."], "contact": "", "email": "", "institute": ""},
  "REFERENCE": [{"title": "Draft genome of Heterosigma akashiwo", "ab_name": ["Author,A."], "status": "Unpublished"}],
  "DATE": {"hold_date": ""},
  "SOURCE": {"organism": "Heterosigma akashiwo", "strain": "", "mol_type": "genomic DNA"},
  "SOURCE_IDENTIFIER": "strain",
  "ASSEMBLY_GAP": {"linkage_evidence": "paired-ends", "min_gap_length": 10, "gap_type": "within scaffold", "estimated_length": "known"}
}
```

- [ ] **Step 4: `common_organelle.json` を作成**

```json
{
  "DBLINK": {"project": "PRJDBxxxxxx", "biosample": "SAMDxxxxxxxx"},
  "SUBMITTER": {"ab_name": ["Author,A."], "contact": "", "email": "", "institute": ""},
  "REFERENCE": [{"title": "Organelle genomes of Heterosigma akashiwo", "ab_name": ["Author,A."], "status": "Unpublished"}],
  "DATE": {"hold_date": ""},
  "SOURCE": {"organism": "Heterosigma akashiwo", "strain": "", "mol_type": "genomic DNA"},
  "SOURCE_IDENTIFIER": "strain"
}
```

- [ ] **Step 5: `sequence_roles.tsv` を作成（organelle のみ記載）**

```
#seq_id	type	seq_name	status	topology
MT	organelle	mitochondrion	complete	circular
CP	organelle	plastid:chloroplast	complete	circular
```

- [ ] **Step 6: Commit（設定は dev。コミット可否はユーザ方針に従う。ここではローカル作成のみでコミットは任意）**

```bash
# 設定ファイルは dev/heterosigma。リポジトリ追加不要方針のためコミットは任意。
ls dev/heterosigma/*.toml dev/heterosigma/common_*.json dev/heterosigma/sequence_roles.tsv
```

---

## Task 12: `verify_agat.py` + Step1 AGAT 実行

**Files:**
- Create: `dev/heterosigma/scripts/verify_agat.py`
- Test: `dev/heterosigma/scripts/tests/test_verify_agat.py`

**Interfaces:**
- CLI: `python verify_agat.py --gff standardized.agat.gff3 [--organelle-seqids MT,CP]` → 検査結果を stdout、問題があれば exit code 1
- Produces: `check(lines) -> dict`（`orphans`, `dup_ids`, `empty_ids`, `organelle_gene_without_mrna`, `parentless_rna` の件数）

- [ ] **Step 1: 失敗テストを書く**

```python
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from verify_agat import check


def test_check_flags_orphan_and_parentless_rna():
    lines = [
        "##gff-version 3\n",
        "c1\tS\tgene\t1\t9\t.\t+\t.\tID=g1\n",
        "c1\tS\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1\n",
        "c1\tS\texon\t1\t9\t.\t+\t.\tID=e1;Parent=missing\n",   # orphan
        "c1\tInfernal\tncRNA\t20\t30\t.\t+\t.\tID=n1\n",         # parentless RNA
    ]
    r = check(lines, {"MT", "CP"})
    assert r["orphans"] >= 1
    assert r["parentless_rna"] >= 1
```

- [ ] **Step 2: 失敗確認**

Run: `cd dev/heterosigma/scripts && uv run pytest tests/test_verify_agat.py -v`
Expected: FAIL

- [ ] **Step 3: 実装**

```python
#!/usr/bin/env python3
"""Verify AGAT-standardized GFF3: hierarchy integrity + organelle/parentless-RNA coverage."""
from __future__ import annotations

import argparse
import re
import sys


def _attr(col9: str, key: str) -> str | None:
    m = re.search(rf"(?:^|;){key}=([^;]+)", col9)
    return m.group(1) if m else None


def check(lines: list[str], organelle_seqids: set[str]) -> dict:
    ids: set[str] = set()
    dup = 0
    empty = 0
    records = []  # (seqid, type, id, parents)
    for line in lines:
        if line.startswith("#") or "\t" not in line:
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        fid = _attr(cols[8], "ID")
        parents = _attr(cols[8], "Parent")
        parents = parents.split(",") if parents else []
        records.append((cols[0], cols[2], fid, parents))
        if cols[2] in ("gene", "mRNA", "tRNA", "rRNA", "ncRNA", "CDS", "exon"):
            if fid == "" or fid is None:
                if cols[2] in ("gene", "mRNA"):
                    empty += 1
            elif fid in ids:
                dup += 1
            else:
                ids.add(fid)

    orphans = sum(1 for _, _, _, ps in records if ps and not all(p in ids for p in ps))
    parentless_rna = sum(1 for _, t, _, ps in records if t in ("ncRNA", "tRNA", "rRNA") and not ps)
    gene_ids_with_mrna = {p for _, t, _, ps in records if t == "mRNA" for p in ps}
    organelle_gene_without_mrna = sum(
        1 for sid, t, fid, _ in records
        if t == "gene" and sid in organelle_seqids and fid not in gene_ids_with_mrna
    )
    return {
        "dup_ids": dup, "empty_ids": empty, "orphans": orphans,
        "parentless_rna": parentless_rna,
        "organelle_gene_without_mrna": organelle_gene_without_mrna,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--organelle-seqids", default="MT,CP")
    args = ap.parse_args()
    org = {s.strip() for s in args.organelle_seqids.split(",") if s.strip()}
    with open(args.gff, encoding="utf-8", errors="replace") as fh:
        r = check(fh.readlines(), org)
    for k, v in r.items():
        print(f"{k}: {v}")
    problems = r["dup_ids"] or r["empty_ids"] or r["orphans"]
    print("PASS" if not problems else "FAIL: hierarchy problems (consider fallback strategy 2/1)")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: テストが通ることを確認**

Run: `cd dev/heterosigma/scripts && uv run pytest tests/test_verify_agat.py -v`
Expected: PASS

- [ ] **Step 5: Step1 AGAT を実行（amd64 コンテナ）+ 検証ゲート**

Run（コンテナ内、要 AGAT）:
```bash
cd dev/heterosigma
agat_convert_sp_gxf2gxf.pl -g braker_with_ncRNA_mtcp.gff3 -o standardized.agat.gff3
uv run python scripts/verify_agat.py --gff standardized.agat.gff3
```
Expected: `PASS`。`transcript`→`mRNA` 改名の有無・organelle_gene_without_mrna=0・parentless_rna の減少を確認。
- FAIL の場合: **戦略2**（`agat_sp_...` を organelle 抽出分のみに適用し核はそのまま結合）→ なお不足なら **戦略1**（Step3 の `pass_coerce_transcript_to_mrna` + Task 4 の parentless-RNA 経路で構造を吸収）。判断メモを `dev/heterosigma/agat_notes.md` に残す。

- [ ] **Step 6: Commit（スクリプトとテストのみ）**

```bash
git add dev/heterosigma/scripts/verify_agat.py dev/heterosigma/scripts/tests/test_verify_agat.py
git commit -m "feat(heterosigma): verify_agat.py hierarchy/coverage gate for Step1"
```

---

## Task 13: エンドツーエンド統合（小サブセット → 全ゲノム）

**Files:**
- Create: `dev/heterosigma/scripts/run_pipeline.sh`（オーケストレーション）
- Test: `tests/test_heterosigma_e2e.py`（`@pytest.mark.slow`、小サブセット）

**Interfaces:**
- Consumes: Task 8–12 の全成果物

- [ ] **Step 1: 小サブセットの slow 統合テストを書く**

`tests/test_heterosigma_e2e.py`（数 scaffold + MT + CP を含む小 fixture を `tests/fixtures/heterosigma_mini.gff3` / `.fa` として用意し、Step2→3→4 を通す）:

```python
import os, subprocess, sys
import pytest

pytestmark = pytest.mark.slow
SCRIPTS = os.path.abspath("dev/heterosigma/scripts")
_MSS_TOOLS = os.environ.get("DDBJ_MSS_TOOLS_SRC",
    os.path.abspath("../ddbj_mss_tools/src"))


@pytest.mark.skipif(not os.path.isdir(os.path.join(_MSS_TOOLS, "common")),
                    reason="ddbj_mss_tools not available")
def test_mini_pipeline_produces_ann(tmp_path):
    gff = "tests/fixtures/heterosigma_mini.gff3"
    fa = "tests/fixtures/heterosigma_mini.fa"
    nuc = tmp_path / "nuclear.gff3"
    org = tmp_path / "organelle.gff3"
    subprocess.run([sys.executable, f"{SCRIPTS}/split_by_compartment.py",
                    "--gff", gff, "--nuclear", str(nuc), "--organelle", str(org),
                    "--drop-pseudogene"], check=True)
    out = tmp_path / "organelle"
    subprocess.run([sys.executable, f"{SCRIPTS}/make_ann.py",
                    "--gff", str(org), "--fasta", fa,
                    "--mss-config", "dev/heterosigma/organelle.mss.toml",
                    "--common", "dev/heterosigma/common_organelle.json",
                    "--sequence-roles", "dev/heterosigma/sequence_roles.tsv",
                    "--submission-category", "GNM", "--out", str(out)], check=True)
    text = (tmp_path / "organelle.ann").read_text()
    assert "COMMON" in text and "source" in text
```

（fixture は本物の mtcp GFF/FASTA から MT・CP 全体 + scaffold_1 の先頭数遺伝子を抜き出して作成する。作成手順を Step2 に記述。）

- [ ] **Step 2: 小 fixture を作成**

Run:
```bash
cd dev/heterosigma
# MT/CP 全行 + scaffold_1 の先頭 200 行を抽出
{ grep -P '^(MT|CP)\t' braker_with_ncRNA_mtcp.gff3; grep -P '^scaffold_1\t' braker_with_ncRNA_mtcp.gff3 | head -200; } > /tmp/mini.body
{ echo '##gff-version 3'; cat /tmp/mini.body; } > ../../tests/fixtures/heterosigma_mini.gff3
# 対応 FASTA（MT, CP, scaffold_1）を抽出
uv run python -c "from Bio import SeqIO; import gzip; recs={r.id:r for r in SeqIO.parse(gzip.open('Haka_JPv1.fa.gz','rt'),'fasta')}; SeqIO.write([recs[i] for i in ('MT','CP','scaffold_1') if i in recs], '../../tests/fixtures/heterosigma_mini.fa','fasta')"
```
Expected: `tests/fixtures/heterosigma_mini.{gff3,fa}` 生成

- [ ] **Step 3: 統合テストを実行**

Run: `uv run pytest tests/test_heterosigma_e2e.py -v -m slow`
Expected: PASS

- [ ] **Step 4: `run_pipeline.sh` を作成（全ゲノム用オーケストレーション）**

```bash
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."   # dev/heterosigma
SRC=../../src
# Step1 は別途 AGAT 実行済み前提: standardized.agat.gff3
python scripts/verify_agat.py --gff standardized.agat.gff3

# Step2
python scripts/build_product_map.py --kaas annotation_kaas.tsv --out product_map.tsv
python scripts/split_by_compartment.py --gff standardized.agat.gff3 \
    --nuclear nuclear.gff3 --organelle organelle.gff3 --drop-pseudogene

# Step3 (normalize) — gzip FASTA
PYTHONPATH=$SRC python -m ddbj_gff.normalize --gff nuclear.gff3 --fasta Haka_JPv1.fa.gz \
    --taxid 2829 --out nuclear.normalized.gff3 --report nuclear.normalize.txt
PYTHONPATH=$SRC python -m ddbj_gff.normalize --gff organelle.gff3 --fasta Haka_JPv1.fa.gz \
    --taxid 2829 --out organelle.normalized.gff3 --report organelle.normalize.txt

# Step4 (make_ann)
mkdir -p submission
python scripts/make_ann.py --gff nuclear.normalized.gff3 --fasta Haka_JPv1.fa.gz \
    --mss-config nuclear.mss.toml --common common_nuclear.json \
    --sequence-roles sequence_roles.tsv --submission-category WGS --out submission/nuclear
python scripts/make_ann.py --gff organelle.normalized.gff3 --fasta Haka_JPv1.fa.gz \
    --mss-config organelle.mss.toml --common common_organelle.json \
    --sequence-roles sequence_roles.tsv --submission-category GNM --out submission/organelle
echo "Done: submission/nuclear.ann, submission/organelle.ann"
```

- [ ] **Step 5: 全ゲノム実走（amd64 コンテナ）+ 目視/validate 確認**

Run（コンテナ）:
```bash
cd dev/heterosigma && bash scripts/run_pipeline.sh
PYTHONPATH=../../src python -m ddbj_gff.validate --gff nuclear.normalized.gff3 || true
head -50 submission/organelle.ann
grep -c $'\tCDS\t' submission/nuclear.ann
grep -c $'\tmisc_feature\t' submission/organelle.ann   # 内部stop 遺伝子
```
Expected: 核 CDS 多数、organelle に MT/CP の CDS/tRNA、TOPOLOGY circular、hypothetical protein 多数、EC 除去済み product。診断（no-cds/internal-stop 等）をレビュー。

- [ ] **Step 6: Commit（スクリプトと fixture、テスト）**

```bash
git add dev/heterosigma/scripts/run_pipeline.sh tests/test_heterosigma_e2e.py tests/fixtures/heterosigma_mini.gff3 tests/fixtures/heterosigma_mini.fa
git commit -m "test(heterosigma): end-to-end mini pipeline + full-genome orchestration"
```

---

## Task 14: 最終検証・登録前チェックリスト

**Files:**
- Create: `dev/heterosigma/SUBMISSION_CHECKLIST.md`

- [ ] **Step 1: 全 mss/normalize テストを実行**

Run: `uv run pytest -q`
Expected: PASS（slow 除く）。加えて `uv run pytest -q -m slow`（ddbj_mss_tools 有り環境）。

- [ ] **Step 2: 登録前チェックリストを作成**

`dev/heterosigma/SUBMISSION_CHECKLIST.md` に以下を列挙:
- [ ] locus_tag prefix `HAKA` を正式 prefix に差し替え（`*.mss.toml`）
- [ ] `common_*.json` の BioProject/BioSample/submitter/reference/hold_date を実値化
- [ ] MT の table 1 翻訳で内部 stop が出た遺伝子（misc_feature 化）の妥当性レビュー
- [ ] AGAT 検証ゲートの結果と採用戦略（3/2/1）を記録
- [ ] `submission/*.ann` を DDBJ の MSS チェック（trans-check 等）にかける
- [ ] FASTA seqid が .ann の entry と一致

- [ ] **Step 3: Commit**

```bash
git add dev/heterosigma/SUBMISSION_CHECKLIST.md
git commit -m "docs(heterosigma): submission pre-flight checklist"
```

---

## Self-Review

**Spec coverage:**
- §Step1 AGAT + 検証ゲート → Task 12 ✅
- §Step2 build_product_map / split(+pseudogene 除外) → Task 8, 9 ✅
- §Step3 transcript→mRNA + gzip + directive → Task 7（directive は既存 normalize で付与）✅
- §Step4 feature=mss / COMMON・source=common → Task 3(feature-only), 10(adapter) ✅
- §5 product ルール → Task 1 ✅ / feature-only API → Task 3 ✅ / 親なし RNA → Task 4 ✅ / pseudogene → Task 5 ✅ / misc_feature → Task 6 ✅ / gzip・transl_table 属性優先(既存挙動) → Task 7 ✅
- §6 COMMON/source(common.json + sequence_roles) → Task 10, 11 ✅
- §7 設定 → Task 11 ✅
- §9 テスト → 各 Task に内包 + Task 13 ✅
- §2 run 構成（核 WGS / organelle GNM） → Task 11, 13 ✅

**Placeholder scan:** COMMON json の `PRJDBxxxxxx` 等は「登録前に実値化」する意図的プレースホルダ（Task 14 で明示）。コード steps に TBD なし。

**Type consistency:** `_product(mrna, gene, cfg)`、`build_entry_features(doc, seqs, cfg, diagnostics) -> dict`、`build_rna_feature(rna, locus_tag, seqlen, gene_id, tx_id, gene_name=None)`、`feature_rows(feat) -> list[list[str]]`、`load_product_map(path) -> dict`、`build_ann_text(...) -> (str, dict)` — Task 間で一致。`MssConfig.product_map`/`product_map_path` は Task 1 で定義し Task 3/10 で使用。

**修正済み:** Task 13 Step2 の FASTA 抽出 one-liner から不要な `bgzf` import を削除済み。
