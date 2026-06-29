# GFF3 → DDBJ MSS — Phase 2 残り機能（複数 transcript / 非コードRNA）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 既存 `ddbj_gff.mss`（minimal MVP, merged）を拡張し、複数 transcript（`minimal`/`nonredundant`/`full` 切替）と非コードRNA（miRNA/pre_miRNA ほか）に対応する。

**Architecture:** 案A。`convert` の per-gene 処理を新 `build_gene_features(gene, mode, assigner, genome_seq, cfg, diagnostics)` に置換。protein-coding はモード別（minimal=代表のみ / full=全 transcript / nonredundant=全 mRNA＋CDSを遺伝子座内 location 重複排除）、non-coding は型→MSS feature マッピング。既存ビルダ（`build_mrna_feature`/`build_cds_feature`/`build_insdc_location`/`collect_spans`/`_submitter_note`/`_representative_mrna`）を再利用。

**Tech Stack:** Python 3.11+ / 既存 `ddbj_gff`＋`ddbj_gff.mss` / Biopython / stdlib tomllib / pytest / dev コンテナ。

## Global Constraints

- Python `>=3.11`、`from __future__ import annotations`。実行依存 biopython のみ（+stdlib tomllib）。
- 変更は `src/ddbj_gff/mss/` 配下のみ。診断は `ddbj_gff.errors`（`Severity`/`Diagnostic`/`GffParseError`）を再利用。
- `transcript_mode` ∈ {`minimal`,`nonredundant`,`full`}、**既定 `nonredundant`**。config `[transcript].mode` ＋ CLI `--mode`（CLI 優先）。不正値 → `GffParseError`（code `invalid-mode`）。
- UTR feature は出力しない。オルガネラ特殊ケース（trans_splicing/transl_except/環状）は対象外（Phase 3）。
- CDS 重複排除（nonredundant）は**遺伝子座内・CDS の location 文字列キー**。共有時は CDS の `note` を全 transcript id 列挙に書き換え。
- locus_tag は1遺伝子座につき1つ（ゲノム全体連番アサイナ、既存）。feature を1つも出さない遺伝子は locus_tag を消費しない（lazy 割当）。
- **テストは dev コンテナ内で実行**: 各 `pytest`/`python` を `docker exec ddbj-gff-dev uv run …` で。`git` は host。
- 各タスクは「失敗テスト→失敗確認→最小実装→成功確認→コミット」。

## 既存コードの前提（変更対象の現状）

- `src/ddbj_gff/mss/convert.py`: `collect_spans` / `_ordered` / `build_insdc_location` / `extract_seq` / `build_source_feature` / `_submitter_note(gene, feat)` / `mrna_partial_flags` / `build_mrna_feature(mrna, gene, locus_tag, seqlen)` / `_product` / `build_cds_feature(mrna, gene, locus_tag, genome_seq, cfg, diagnostics)->MssFeature|None` / `_representative_mrna(gene, diagnostics)->Feature|None`（`.1`優先・複数で `multi-transcript` WARNING）/ `_span_start(feature)` / `convert(doc, seqs, cfg, common_rows, *, strict=False)`。`assigner = LocusTagAssigner.from_config(cfg)` は seqid ループ**前**にある（ゲノム全体連番）。
- `src/ddbj_gff/mss/config.py`: `@dataclass MssConfig`（mutable）、`load_config(path)->(MssConfig, list[Diagnostic])`、`_KNOWN_SECTIONS`。
- `src/ddbj_gff/mss/cli.py`: `main(argv)->int`。
- Phase 1 `Feature`: `.id/.type/.children/.spans/.attributes/.gene/.product/.note/.locus_tag`。`Feature.note`→`attributes.get("Note",[])`。

---

## File Structure

| ファイル | 変更 |
|---|---|
| `src/ddbj_gff/mss/config.py` | `transcript_mode` フィールド＋`[transcript].mode` 読込＋検証 |
| `src/ddbj_gff/mss/convert.py` | `_RNA_MAP`/`build_noncoding_features`/`_set_submitter_transcripts`/`build_gene_features` 追加、`convert` 内ループ置換 |
| `src/ddbj_gff/mss/cli.py` | `--mode` 引数＋cfg 上書き |
| `tests/test_mss_config.py`, `tests/test_mss_convert.py`, `tests/test_mss_cli.py`, `tests/test_mss_snapshot.py`, `tests/test_mss_integration.py`, `tests/mss_fixtures/` | テスト・フィクスチャ拡張 |

---

## Task 1: config に transcript_mode

**Files:** Modify `src/ddbj_gff/mss/config.py`; Test `tests/test_mss_config.py`

**Interfaces:**
- Produces: `MssConfig.transcript_mode: str = "nonredundant"`; `load_config` reads `[transcript].mode`, validates ∈ {minimal,nonredundant,full} (else `GffParseError` code `invalid-mode`), defaults `nonredundant`; `"transcript"` added to `_KNOWN_SECTIONS`.

- [ ] **Step 1: 失敗するテストを追加** (`tests/test_mss_config.py` に追記)

```python
def test_transcript_mode_default_and_explicit(tmp_path):
    p = tmp_path / "c.toml"
    p.write_text('[source]\norganism="O"\nmol_type="genomic DNA"\n[locus_tag]\nprefix="P"\n')
    cfg, _ = load_config(str(p))
    assert cfg.transcript_mode == "nonredundant"
    p.write_text('[source]\norganism="O"\nmol_type="genomic DNA"\n[locus_tag]\nprefix="P"\n[transcript]\nmode="full"\n')
    cfg, _ = load_config(str(p))
    assert cfg.transcript_mode == "full"


def test_invalid_transcript_mode_raises(tmp_path):
    p = tmp_path / "c.toml"
    p.write_text('[source]\norganism="O"\n[locus_tag]\nprefix="P"\n[transcript]\nmode="bogus"\n')
    import pytest
    from ddbj_gff.errors import GffParseError
    with pytest.raises(GffParseError):
        load_config(str(p))
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_mss_config.py -k transcript -v` → FAIL（`transcript_mode` 属性なし）

- [ ] **Step 3: 実装**

`config.py` の `MssConfig` に末尾フィールドを追加:
```python
    transcript_mode: str = "nonredundant"
```
`_KNOWN_SECTIONS` に `"transcript"` を追加:
```python
_KNOWN_SECTIONS = {"source", "locus_tag", "cds", "assembly_gap", "product", "transcript"}
```
`load_config` の `cfg = MssConfig(...)` 構築の直前に追加:
```python
    transcript = data.get("transcript", {})
    mode = transcript.get("mode", "nonredundant")
    if mode not in ("minimal", "nonredundant", "full"):
        raise GffParseError(Diagnostic(Severity.ERROR, None, "invalid-mode",
                                       f"transcript mode {mode!r} must be one of minimal/nonredundant/full"))
```
そして `MssConfig(...)` 呼び出しに `transcript_mode=mode,` を追加。

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_mss_config.py -v` → all pass

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/mss/config.py tests/test_mss_config.py
git commit -m "feat(mss): transcript_mode config setting with validation"
```

---

## Task 2: 非コードRNA ビルダ

**Files:** Modify `src/ddbj_gff/mss/convert.py`; Test `tests/test_mss_noncoding.py`

**Interfaces:**
- Consumes: `collect_spans`/`build_insdc_location`/`_submitter_note`（同モジュール）, `MssFeature`/`MssQualifier`。
- Produces（`convert.py` に追記）:
  - module-level `_RNA_MAP: dict[str,str]`.
  - `build_noncoding_features(gene, locus_tag: str, seqlen: int, cfg) -> list[MssFeature]`：`gene.children` のうち `_RNA_MAP` にある型を MSS feature 化（pre_miRNA→precursor_RNA、miRNA→ncRNA[+ncRNA_class]、tRNA/rRNA/ncRNA…）。location は exon 連結 or feature 自身の span。partial 判定なし。recognized RNA 子が無ければ `[]`。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_mss_noncoding.py`:
```python
from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.convert import build_noncoding_features


def cfg():
    return MssConfig(source={})


def test_mirna_gene_maps_to_precursor_and_ncrna():
    gene = Feature("Mp1g00675", "S", "gene", [Span("chr1", 603255, 603422, "-")],
                   {"gene_biotype": ["miRNA"]}, [])
    pre = Feature("Mp1g00675.pre", "S", "pre_miRNA", [Span("chr1", 603255, 603422, "-")],
                  {"Note": ["Mpo-pre-miR11669"]}, [])
    m1 = Feature("Mp1g00675.1", "S", "miRNA", [Span("chr1", 603382, 603402, "-")],
                 {"Note": ["Mpo-miR11669.1"]}, [])
    m2 = Feature("Mp1g00675.2", "S", "miRNA", [Span("chr1", 603384, 603405, "-")],
                 {"Note": ["Mpo-miR11669.2"]}, [])
    gene.children = [pre, m1, m2]
    feats = build_noncoding_features(gene, "PFX_000010", 700000, cfg())
    assert [f.key for f in feats] == ["precursor_RNA", "ncRNA", "ncRNA"]
    assert feats[0].location == "complement(603255..603422)"
    nc = {q.key: q.value for q in feats[1].qualifiers}
    assert nc["locus_tag"] == "PFX_000010"
    assert nc["ncRNA_class"] == "miRNA"
    assert any(q.key == "note" and q.value == "Mpo-miR11669.1" for q in feats[1].qualifiers)
    assert any(q.key == "note" and "submitter_gene_id: Mp1g00675" in q.value for q in feats[1].qualifiers)


def test_no_recognized_rna_children_returns_empty():
    gene = Feature("g", "S", "gene", [Span("chr1", 1, 9, "+")], {}, [])
    gene.children = []
    assert build_noncoding_features(gene, "PFX_000010", 1000, cfg()) == []
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_mss_noncoding.py -v` → FAIL（`ImportError`）

- [ ] **Step 3: 実装** — `convert.py` 末尾に追記:
```python
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


def build_noncoding_features(gene, locus_tag: str, seqlen: int, cfg) -> list:
    features = []
    for rna in gene.children:
        feat_key = _RNA_MAP.get(rna.type)
        if feat_key is None:
            continue
        spans = collect_spans(rna, "exon") or rna.spans
        location = build_insdc_location(spans, seqlen)
        quals = [MssQualifier("locus_tag", locus_tag)]
        if feat_key == "ncRNA":
            quals.append(MssQualifier("ncRNA_class", rna.type if rna.type != "ncRNA" else "other"))
        if rna.product:
            quals.append(MssQualifier("product", rna.product))
        if gene.gene or rna.gene:
            quals.append(MssQualifier("gene", gene.gene or rna.gene))
        for note_val in rna.note:
            quals.append(MssQualifier("note", note_val))
        quals.append(_submitter_note(gene, rna))
        features.append(MssFeature(feat_key, location, quals))
    return features
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_mss_noncoding.py -v` → 2 passed

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/mss/convert.py tests/test_mss_noncoding.py
git commit -m "feat(mss): non-coding RNA feature builder (miRNA/pre_miRNA/etc.)"
```

---

## Task 3: build_gene_features（3モード＋CDS重複排除）

**Files:** Modify `src/ddbj_gff/mss/convert.py`; Test `tests/test_mss_gene_features.py`

**Interfaces:**
- Consumes: `build_mrna_feature`/`build_cds_feature`/`collect_spans`/`_representative_mrna`/`_RNA_MAP`/`build_noncoding_features`/`_submitter_note`, `LocusTagAssigner`, `Diagnostic`/`Severity`。
- Produces（`convert.py` に追記）:
  - `_set_submitter_transcripts(cds: MssFeature, gene, transcript_ids: list[str]) -> None`：CDS の submitter `note` 値を全 transcript id 列挙に書き換え。
  - `build_gene_features(gene, mode: str, assigner, genome_seq, cfg, diagnostics) -> list[MssFeature]`：mRNA 子があれば protein-coding（mode別）、無ければ非コード（`_RNA_MAP` 子→`build_noncoding_features`、無ければ `no-rna` WARNING）。locus_tag は feature を出す時のみ lazy 割当。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_mss_gene_features.py`:
```python
from Bio.Seq import Seq
from ddbj_gff.model import Feature, Span
from ddbj_gff.mss.config import MssConfig
from ddbj_gff.mss.locus_tag import LocusTagAssigner
from ddbj_gff.mss.convert import build_gene_features


def cfg():
    return MssConfig(source={}, locus_tag_prefix="PFX")


def assigner():
    return LocusTagAssigner.from_config(cfg())


def two_transcript_gene(cds_same=True):
    # gene with two mRNAs; both share identical CDS (1..9) but differ in exon (UTR)
    gene = Feature("g", "S", "gene", [Span("c", 1, 40, "+")], {}, [])
    t1 = Feature("g.1", "S", "mRNA", [Span("c", 1, 30, "+")], {}, [])
    t1.children = [Feature("e1", "S", "exon", [Span("c", 1, 30, "+")], {}, []),
                   Feature("c1", "S", "CDS", [Span("c", 1, 9, "+", 0)], {}, [])]
    t2 = Feature("g.2", "S", "mRNA", [Span("c", 1, 40, "+")], {}, [])
    cds2 = Span("c", 1, 9, "+", 0) if cds_same else Span("c", 1, 12, "+", 0)
    t2.children = [Feature("e2", "S", "exon", [Span("c", 1, 40, "+")], {}, []),
                   Feature("c2", "S", "CDS", [cds2], {}, [])]
    gene.children = [t1, t2]
    return gene


def keys(feats):
    return [f.key for f in feats]


def test_minimal_keeps_one_transcript():
    g = two_transcript_gene()
    feats = build_gene_features(g, "minimal", assigner(), Seq("ATGAAATAA" + "C" * 31), cfg(), [])
    assert keys(feats) == ["mRNA", "CDS"]


def test_full_keeps_all_transcripts_and_cds():
    g = two_transcript_gene()
    feats = build_gene_features(g, "full", assigner(), Seq("ATGAAATAA" + "C" * 31), cfg(), [])
    assert keys(feats) == ["mRNA", "CDS", "mRNA", "CDS"]


def test_nonredundant_dedupes_shared_cds():
    g = two_transcript_gene(cds_same=True)
    feats = build_gene_features(g, "nonredundant", assigner(), Seq("ATGAAATAA" + "C" * 31), cfg(), [])
    assert keys(feats) == ["mRNA", "mRNA", "CDS"]   # both mRNAs, one shared CDS
    cds = feats[2]
    note = [q.value for q in cds.qualifiers if q.key == "note"][0]
    assert "g.1" in note and "g.2" in note            # note lists both source transcripts


def test_nonredundant_keeps_distinct_cds():
    g = two_transcript_gene(cds_same=False)
    feats = build_gene_features(g, "nonredundant", assigner(), Seq("ATGAAATAA" + "C" * 31), cfg(), [])
    assert keys(feats) == ["mRNA", "mRNA", "CDS", "CDS"]  # different CDS -> not deduped


def test_noncoding_gene_dispatch():
    gene = Feature("m", "S", "gene", [Span("c", 1, 50, "-")], {"gene_biotype": ["miRNA"]}, [])
    gene.children = [Feature("m.pre", "S", "pre_miRNA", [Span("c", 1, 50, "-")], {}, []),
                     Feature("m.1", "S", "miRNA", [Span("c", 10, 30, "-")], {}, [])]
    feats = build_gene_features(gene, "nonredundant", assigner(), Seq("A" * 100), cfg(), [])
    assert keys(feats) == ["precursor_RNA", "ncRNA"]


def test_gene_with_no_rna_or_mrna_warns():
    gene = Feature("x", "S", "gene", [Span("c", 1, 9, "+")], {}, [])
    gene.children = []
    diags = []
    feats = build_gene_features(gene, "nonredundant", assigner(), Seq("A" * 9), cfg(), diags)
    assert feats == []
    assert any(d.code == "no-rna" for d in diags)
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_mss_gene_features.py -v` → FAIL（`ImportError`）

- [ ] **Step 3: 実装** — `convert.py` 末尾に追記:
```python
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
        if not any(c.type in _RNA_MAP for c in gene.children):
            diagnostics.append(Diagnostic(Severity.WARNING, None, "no-rna",
                                          f"gene {gene.id!r} has no mRNA or recognized RNA child; skipped"))
            return []
        return build_noncoding_features(gene, assigner.assign(gene), len(genome_seq), cfg)

    transcripts = sorted(transcripts, key=lambda m: m.id or "")
    if mode == "minimal":
        rep = _representative_mrna(gene, diagnostics)
        transcripts = [rep] if rep is not None else []

    features = []
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
        features.append(build_mrna_feature(mrna, gene, locus_tag, len(genome_seq)))
        cds = build_cds_feature(mrna, gene, locus_tag, genome_seq, cfg, diagnostics)
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
    return features
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_mss_gene_features.py -v` → 6 passed

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/mss/convert.py tests/test_mss_gene_features.py
git commit -m "feat(mss): mode-aware per-gene builder with CDS dedup and non-coding dispatch"
```

---

## Task 4: convert 配線 ＋ CLI --mode

**Files:** Modify `src/ddbj_gff/mss/convert.py`, `src/ddbj_gff/mss/cli.py`, `tests/test_mss_convert.py`, `tests/test_mss_cli.py`

**Interfaces:**
- Consumes: `build_gene_features`（Task 3）, `cfg.transcript_mode`（Task 1）。
- Produces: `convert` の per-gene 処理を `build_gene_features(gene, cfg.transcript_mode, assigner, genome_seq, cfg, diagnostics)` に置換。CLI に `--mode {minimal,nonredundant,full}`（指定時 `cfg.transcript_mode` を上書き）。

- [ ] **Step 1: テストを更新/追加**

`tests/test_mss_convert.py` の既存 `test_multi_transcript_warns_and_keeps_first` を **minimal モード明示**に変更（既定が nonredundant になったため）。その関数の `convert(doc, seqs, cfg(), ["COMMON"])` を次に変更:
```python
    c = cfg(); c.transcript_mode = "minimal"
    mss, diags = convert(doc, seqs, c, ["COMMON"])
```
（残りのアサーション `multi-transcript` WARNING と mRNA/CDS 各1 はそのまま。）

`tests/test_mss_convert.py` に追記:
```python
def test_convert_default_mode_is_nonredundant_keeps_all_transcripts():
    gff = (
        "##gff-version 3\n"
        "chr1\tS\tgene\t1\t40\t.\t+\t.\tID=g\n"
        "chr1\tS\tmRNA\t1\t30\t.\t+\t.\tID=g.1;Parent=g\n"
        "chr1\tS\texon\t1\t30\t.\t+\t.\tID=e1;Parent=g.1\n"
        "chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=g.1\n"
        "chr1\tS\tmRNA\t1\t40\t.\t+\t.\tID=g.2;Parent=g\n"
        "chr1\tS\texon\t1\t40\t.\t+\t.\tID=e2;Parent=g.2\n"
        "chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=c2;Parent=g.2\n"
    )
    doc = parse(gff)
    mss, diags = convert(doc, {"chr1": Seq("ATGAAATAA" + "C" * 31)}, cfg(), ["COMMON"])
    keys = [f.key for f in mss.entries[0].features]
    assert keys == ["source", "mRNA", "mRNA", "CDS"]  # nonredundant: 2 mRNA, shared CDS once
    assert not any(d.code == "multi-transcript" for d in diags)  # no warning in nonredundant
```

`tests/test_mss_cli.py` に追記（既存 GFF/CONFIG/COMMON 定数を再利用）:
```python
def test_cli_mode_flag_overrides(tmp_path):
    (tmp_path / "g.gff").write_text(GFF)
    (tmp_path / "g.fa").write_text(FASTA)
    (tmp_path / "c.toml").write_text(CONFIG)
    (tmp_path / "common.tsv").write_text(COMMON)
    out = tmp_path / "r"
    rc = main(["--gff", str(tmp_path/"g.gff"), "--fasta", str(tmp_path/"g.fa"),
               "--config", str(tmp_path/"c.toml"), "--common", str(tmp_path/"common.tsv"),
               "--out", str(out), "--mode", "minimal"])
    assert rc == 0  # GFF has a single transcript so output is identical, but --mode is accepted
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_mss_convert.py tests/test_mss_cli.py -v` → 新規/更新テストが FAIL（`--mode` 未知、default 挙動が minimal のまま等）

- [ ] **Step 3: 実装**

`convert.py`：`convert` 内の `for gene in genes:` ブロック（`_representative_mrna`～`build_cds_feature` の一連）を次に置換:
```python
        for gene in genes:
            features.extend(build_gene_features(gene, cfg.transcript_mode, assigner,
                                                genome_seq, cfg, diagnostics))
```
（`assigner = LocusTagAssigner.from_config(cfg)` は seqid ループ前のまま。`_representative_mrna`/`_span_start` は残す。）

`cli.py`：`ap.add_argument("--out", required=True)` の後に追加:
```python
    ap.add_argument("--mode", choices=["minimal", "nonredundant", "full"], default=None)
```
`cfg, cfg_diags = load_config(args.config)` の直後に追加:
```python
    if args.mode:
        cfg.transcript_mode = args.mode
```

- [ ] **Step 4: 成功確認**
`docker exec ddbj-gff-dev uv run pytest tests/test_mss_convert.py tests/test_mss_cli.py -v` → all pass。
全体回帰: `docker exec ddbj-gff-dev uv run pytest -q` → all pass。

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/mss/convert.py src/ddbj_gff/mss/cli.py tests/test_mss_convert.py tests/test_mss_cli.py
git commit -m "feat(mss): wire transcript_mode into convert and CLI --mode"
```

---

## Task 5: ゴールドスナップショット拡張

**Files:** Modify `tests/mss_fixtures/mini.gff3`, `tests/mss_fixtures/mini.fa`; regenerate `tests/mss_fixtures/expected.ann`/`expected.fasta`; Test `tests/test_mss_snapshot.py`

**Interfaces:**
- Consumes: `cli.main`。Produces: 既定 `nonredundant` で複数 transcript＋miRNA を含む E2E ゴールド。

- [ ] **Step 1: フィクスチャ拡張（TAB 厳守）**

`tests/mss_fixtures/mini.gff3` に、既存 g1/g2 の後（FASTA ディレクティブは無い）へ次を追記（列はタブ。座標は既存 60bp 配列に収める）— **複数 transcript 遺伝子 g3（chr1 上、+鎖、CDS同一/UTR差の2 transcript）** と **miRNA 遺伝子 g4**。`##sequence-region` の終端と `mini.fa` 長を新座標に合わせて拡張する必要があるため、`mini.fa` を 120bp に拡張する:

`tests/mss_fixtures/mini.gff3`（全文を次で置き換え。タブ区切り）:
```text
##gff-version 3
##sequence-region chr1 1 120
chr1	S	gene	1	30	.	+	.	ID=g1;gene=MpX
chr1	S	mRNA	1	30	.	+	.	ID=g1.1;Parent=g1;product=widget protein
chr1	S	exon	1	30	.	+	.	ID=e1;Parent=g1.1
chr1	S	CDS	1	9	.	+	0	ID=c1;Parent=g1.1
chr1	S	gene	40	60	.	-	.	ID=g2;locus_tag=PRE_999
chr1	S	mRNA	40	60	.	-	.	ID=g2.1;Parent=g2
chr1	S	exon	40	60	.	-	.	ID=e2;Parent=g2.1
chr1	S	CDS	49	60	.	-	0	ID=c3;Parent=g2.1
chr1	S	gene	61	90	.	+	.	ID=g3
chr1	S	mRNA	61	84	.	+	.	ID=g3.1;Parent=g3
chr1	S	exon	61	84	.	+	.	ID=e31;Parent=g3.1
chr1	S	CDS	61	69	.	+	0	ID=c31;Parent=g3.1
chr1	S	mRNA	61	90	.	+	.	ID=g3.2;Parent=g3
chr1	S	exon	61	90	.	+	.	ID=e32;Parent=g3.2
chr1	S	CDS	61	69	.	+	0	ID=c32;Parent=g3.2
chr1	S	gene	100	120	.	+	.	ID=g4;gene_biotype=miRNA
chr1	S	pre_miRNA	100	120	.	+	.	ID=g4.pre;Note=demo-pre;Parent=g4
chr1	S	miRNA	105	115	.	+	.	ID=g4.1;Note=demo-miR;Parent=g4
```

- [ ] **Step 2: mini.fa を 120bp に整合（CDS が翻訳できるよう塩基配置）**

`tests/mss_fixtures/mini.fa` を次で置き換え（chr1, 120bp）。1..9=`ATGAAATAA`(g1 完全CDS)、49..60 の −鎖 revcomp が `ATGAAAGCATAA`(g2)、61..69=`ATGAAATAA`(g3 共有CDS 完全)、残りは A 埋め:
```text
>chr1
ATGAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTATGCTTTCATATGAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```
検証（コドンが意図通りか）:
```
docker exec ddbj-gff-dev uv run python - <<'PY'
from ddbj_gff import parse
from Bio import SeqIO
from ddbj_gff.mss.convert import extract_seq, collect_spans
doc = parse(open("tests/mss_fixtures/mini.gff3").read())
seq = next(SeqIO.parse("tests/mss_fixtures/mini.fa","fasta")).seq
assert len(seq) == 120, len(seq)
for g in doc.roots:
    for m in [c for c in g.children if c.type=="mRNA"]:
        cds = collect_spans(m,"CDS")
        if cds:
            s = extract_seq(cds, seq); print(g.id, m.id, str(s), str(s.translate()))
PY
```
期待: g1.1=`ATGAAATAA`(MK*)、g2.1=`ATGAAAGCATAA`(MKA*)、g3.1/g3.2=`ATGAAATAA`(MK*)。一致しなければ **mini.fa の該当塩基のみ調整**（CDS座標は変えない）。

- [ ] **Step 3: 既定モードで期待出力を再生成しレビュー**
```
docker exec ddbj-gff-dev uv run python -m ddbj_gff.mss \
  --gff tests/mss_fixtures/mini.gff3 --fasta tests/mss_fixtures/mini.fa \
  --config tests/mss_fixtures/config.toml --common tests/mss_fixtures/common.metadata.tsv \
  --out tests/mss_fixtures/expected
```
生成された `expected.ann` を Read で確認: g1（locus_tag MPTK1_000010, CDS 1..9）、g2（PRE_999, complement）、**g3（locus_tag MPTK1_000020, mRNA 2本＋共有 CDS 1本で note に g3.1, g3.2）**、**g4（MPTK1_000030, precursor_RNA＋ncRNA[ncRNA_class=miRNA]）** が妥当であること。妥当ならゴールドとする。

- [ ] **Step 4: スナップショットテスト確認**（`tests/test_mss_snapshot.py` は既存のまま — 同じ fixtures→expected を比較）
`docker exec ddbj-gff-dev uv run pytest tests/test_mss_snapshot.py -v` → 1 passed（got == expected）

- [ ] **Step 5: Commit**
```bash
git add tests/mss_fixtures/mini.gff3 tests/mss_fixtures/mini.fa tests/mss_fixtures/expected.ann tests/mss_fixtures/expected.fasta
git commit -m "test(mss): extend golden snapshot with multi-transcript and miRNA genes"
```

---

## Task 6: 統合テスト更新（実 marchantia, 既定 nonredundant）

**Files:** Modify `tests/test_mss_integration.py`

**Interfaces:** Consumes: `convert`（既定 nonredundant）, 実 marchantia ファイル。

- [ ] **Step 1: テストを更新**

`tests/test_mss_integration.py` の `test_marchantia_converts_without_errors` の本体に、変換後の追加アサーションを足す（`MssConfig(...)` は `transcript_mode` 未指定で既定 `nonredundant` になる点を活かす）。`assert total_cds > 1000` の後に追記:
```python
    # nonredundant default: all transcripts emit an mRNA; CDS is deduped per locus
    total_mrna = sum(1 for e in mss.entries for f in e.features if f.key == "mRNA")
    total_cds = sum(1 for e in mss.entries for f in e.features if f.key == "CDS")
    assert total_mrna >= total_cds  # dedup never produces more CDS than mRNA
    assert total_mrna > 20000       # all ~22k transcripts emit an mRNA
    # miRNA genes produce ncRNA / precursor_RNA features
    rna_keys = {f.key for e in mss.entries for f in e.features}
    assert "ncRNA" in rna_keys and "precursor_RNA" in rna_keys
```

- [ ] **Step 2: slow 実行** — `docker exec ddbj-gff-dev uv run pytest tests/test_mss_integration.py -m slow -v`（実ゲノム、数分可）→ 1 passed。もし `total_mrna > 20000` 等が想定とずれたら、実際の値をログして閾値を実測に合わせて調整（ただし「全 transcript の mRNA」「ncRNA/precursor_RNA 存在」「ERROR 0」の主旨は維持）。

- [ ] **Step 3: 全体確認** — `docker exec ddbj-gff-dev uv run pytest -q`（slow除外）→ all pass。

- [ ] **Step 4: Commit**
```bash
git add tests/test_mss_integration.py
git commit -m "test(mss): integration asserts multi-transcript + ncRNA under nonredundant default"
```

---

## Self-Review

**1. Spec coverage**（spec §→タスク）:
- §3 モード config/CLI → Task 1, 4
- §4.1 protein-coding 3モード＋CDS dedup → Task 3、配線 Task 4
- §4.2 非コードRNA マッピング → Task 2、dispatch Task 3
- §5 診断（multi-transcript=minimalのみ／no-rna／invalid-mode）→ Task 1, 3, 4
- §6 テスト（config/convert/cli/snapshot/integration）→ Task 1-6
- §7 スコープ: オルガネラ等は対象外（実装なし）

ギャップ: なし。

**2. Placeholder scan**: 各ステップに実コード。Task 5 の expected.* は「生成→レビュー→固定」手順（プレースホルダでない）。"TBD"等なし。

**3. Type consistency**: `transcript_mode`（config/cli/convert）、`build_gene_features(gene, mode, assigner, genome_seq, cfg, diagnostics)`、`build_noncoding_features(gene, locus_tag, seqlen, cfg)`、`_set_submitter_transcripts(cds, gene, transcript_ids)`、`_RNA_MAP`、既存 `build_mrna_feature`/`build_cds_feature`/`collect_spans`/`_representative_mrna`/`_submitter_note`/`LocusTagAssigner.from_config` — タスク間で一致。CDS dedup キー = `MssFeature.location`（文字列）。`MssConfig`/`MssQualifier` は mutable dataclass（属性代入可）。
