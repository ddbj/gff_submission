# 正規化器 feature-type 修正（Finding A/B）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `pass_so_terms` の SO-term collapse をコア型 whitelist `{gene, mRNA, CDS, exon, intron}` に絞り、INSDC feature-key への書換をやめて Finding A/B（ゼニゴケ実走で露呈した feature-type-not-insdc 大量発生・pre_miRNA 劣化）を是正する。

**Architecture:** `src/ddbj_gff/normalize/passes.py` の `pass_so_terms` に1つの判定（書換先が whitelist に含まれる時のみ collapse）を追加するだけ。SO→INSDC-feature 変換は Phase 2 生成（既存 `_RNA_MAP` 等）の責務に一本化。`vocab.py`/`pass_directives`/`normalize()`/CLI/`report.py` は不変。

**Tech Stack:** Python 3.11+ / 既存 `ddbj_gff`（Phase1・3-A・3-B） / pytest / dev コンテナ。

## Global Constraints

- 変更は `src/ddbj_gff/normalize/passes.py` の `pass_so_terms` ＋ 関連テストのみ。`vocab.py`（`feature_qualifiers`/dedup）・`pass_directives`・`NormalizeContext`・`normalize()`・CLI・`report.py`・`src/ddbj_gff/validate/data/` は**変更しない**。
- collapse 対象集合は厳密に `_COLLAPSE_TARGETS = {"gene", "mRNA", "CDS", "exon", "intron"}`。
- collapse 条件: `target = vocab.insdc_map.get(f.type)`。`target is None` → `unmapped-type` report（unresolved）して不変。`target == f.type` または `target not in _COLLAPSE_TARGETS` → **放置（書換えない・Change を出さない）**。それ以外（`target in _COLLAPSE_TARGETS and target != f.type`）→ 書換＋qualifier 付与。
- qualifier 付与ロジック（具体値→属性、プレースホルダ→`needs-manual`、既存キー非上書き、フラグ→`["true"]`）は**そのまま維持**（防御的に残す。whitelisted target には実データ上プレースホルダ qualifier が無いため `needs-manual` は実質不発火）。
- 放置型（miRNA/pre_miRNA/UTR/mobile_element/binding_site 等、target が whitelist 外）は report に出さない。`unmapped-type`（既知 SO term でない）のみ report。
- **テストは dev コンテナ内**: `docker exec ddbj-gff-dev uv run pytest …`。`git` は host。
- TDD: 失敗テスト→失敗確認→最小実装→成功確認→コミット。

## File Structure

| ファイル | 責務 | 変更 |
|---|---|---|
| `src/ddbj_gff/normalize/passes.py` | `_COLLAPSE_TARGETS` 追加＋`pass_so_terms` の collapse 判定変更 | Modify |
| `tests/test_normalize_pass_so_terms.py` | ユニットテスト全面差し替え（collapse/放置/未マップ） | Replace |
| `tests/test_normalize_integration.py` | 受け入れテスト1件追加（Finding A 回帰防止） | Append |

## 既存インターフェース（参照）

- `pass_so_terms(doc, ctx) -> list[Change]`（`ctx.vocab` は `Vocab` で `.insdc_map: dict[str,str]` / `.feature_qualifiers: dict[str,tuple[str,...]]`）。ヘルパ `_is_placeholder(qual)`・`_qualifier_to_attr(qual)` は既存・不変。
- `report.Change(action, target, message)`。`Feature.type`（可変）/`.attributes: dict[str,list[str]]`/`.id`。
- `ddbj_gff.parse`、`ddbj_gff.validate.validate(doc) -> list[Diagnostic]`（`d.code`）、`ddbj_gff.normalize.normalize(doc, *, seq_lengths=None, config=None) -> (doc, report)`、`NormalizeConfig(taxid=…)`。
- 確認済みマッピング（`feature-mapping.tsv`）: `coding_exon→exon`, `pseudogenic_CDS→CDS`(+`/pseudo`), `processed_pseudogene→gene`(+`/pseudogene="processed"`), `spliceosomal_intron→intron`, `miRNA→ncRNA`, `pre_miRNA→ncRNA`, `five_prime_UTR→5'UTR`, `binding_site→misc_binding`, `mobile_genetic_element→mobile_element`。全 SO term は validator の `feature_types`（列2）に含まれる。

---

### Task 1: whitelist ゲート ＋ ユニットテスト差し替え

**Files:**
- Modify: `src/ddbj_gff/normalize/passes.py`（`pass_so_terms` 末尾関数）
- Test: `tests/test_normalize_pass_so_terms.py`（全面差し替え）

**Interfaces:**
- Produces: `_COLLAPSE_TARGETS` 定数。`pass_so_terms` の挙動: whitelist 内 target のみ collapse、それ以外放置、未マップのみ report。

- [ ] **Step 1: テストを差し替え（失敗するテスト）**

`tests/test_normalize_pass_so_terms.py` を以下で**全面置換**:
```python
from ddbj_gff.model import Feature, Span, GffDocument
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.normalize.passes import NormalizeContext, pass_so_terms
from ddbj_gff.validate.vocab import load_vocab


def _ctx():
    return NormalizeContext(vocab=load_vocab(), seq_lengths=None, config=NormalizeConfig())


def _doc(*feats):
    return GffDocument(features=list(feats))


# --- collapse to a core whitelist type (gene/mRNA/CDS/exon/intron) ---

def test_pseudogenic_cds_collapses_to_cds_with_pseudo_flag():
    f = Feature("c", "S", "pseudogenic_CDS", [Span("chr1", 1, 9, "+", 0)], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "CDS"
    assert f.attributes.get("pseudo") == ["true"]
    assert any(c.action == "rename-type" for c in changes)
    assert any(c.action == "add-qualifier" for c in changes)


def test_processed_pseudogene_collapses_to_gene_with_qualifier():
    f = Feature("g", "S", "processed_pseudogene", [Span("chr1", 1, 9, "+")], {}, [])
    pass_so_terms(_doc(f), _ctx())
    assert f.type == "gene"
    assert f.attributes.get("pseudogene") == ["processed"]


def test_coding_exon_collapses_to_exon():
    f = Feature("e", "S", "coding_exon", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "exon"
    assert any(c.action == "rename-type" for c in changes)


def test_spliceosomal_intron_collapses_to_intron():
    f = Feature("i", "S", "spliceosomal_intron", [Span("chr1", 1, 9, "+")], {}, [])
    pass_so_terms(_doc(f), _ctx())
    assert f.type == "intron"


def test_same_name_is_noop():
    f = Feature("c", "S", "CDS", [Span("chr1", 1, 9, "+", 0)], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "CDS"
    assert changes == []


def test_existing_attribute_not_clobbered():
    f = Feature("c", "S", "pseudogenic_CDS", [Span("chr1", 1, 9, "+", 0)], {"pseudo": ["existing"]}, [])
    pass_so_terms(_doc(f), _ctx())
    assert f.attributes["pseudo"] == ["existing"]


# --- non-core targets are LEFT for Phase 2 (Finding A/B fix) ---

def test_mirna_left_for_phase2():
    f = Feature("r", "S", "miRNA", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "miRNA"        # unchanged: Phase 2 maps it to ncRNA[ncRNA_class=miRNA]
    assert changes == []


def test_pre_mirna_left_for_phase2():
    f = Feature("r", "S", "pre_miRNA", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "pre_miRNA"    # unchanged: Phase 2 maps it to precursor_RNA
    assert changes == []


def test_five_prime_utr_left_alone():
    f = Feature("u", "S", "five_prime_UTR", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "five_prime_UTR"
    assert changes == []


def test_non_core_insdc_target_left_alone():
    # binding_site -> misc_binding, mobile_genetic_element -> mobile_element:
    # targets not in the core whitelist -> left unchanged, no fabricated qualifier, no Change
    for t in ("binding_site", "mobile_genetic_element"):
        f = Feature("x", "S", t, [Span("chr1", 1, 9, "+")], {}, [])
        changes = pass_so_terms(_doc(f), _ctx())
        assert f.type == t
        assert f.attributes == {}
        assert changes == []


def test_unmapped_type_reported_unchanged():
    f = Feature("x", "S", "totally_made_up_type", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "totally_made_up_type"
    assert any(c.action == "unmapped-type" for c in changes)
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_pass_so_terms.py -v` → `test_mirna_left_for_phase2` 等が FAIL（現状は miRNA→ncRNA に書換えるため）。

- [ ] **Step 3: 実装** — `src/ddbj_gff/normalize/passes.py` の `pass_so_terms` を以下に置換（ヘルパ `_is_placeholder`/`_qualifier_to_attr` は変更しない。`_COLLAPSE_TARGETS` を `pass_so_terms` の直前に追加）:
```python
_COLLAPSE_TARGETS = {"gene", "mRNA", "CDS", "exon", "intron"}


def pass_so_terms(doc, ctx) -> list:
    vocab = ctx.vocab
    changes: list = []
    for f in doc.features:
        target = vocab.insdc_map.get(f.type)
        if target is None:
            changes.append(Change("unmapped-type", f.id or "?",
                                  f"feature type {f.type!r} is not a known SO term; left unchanged"))
            continue
        if target == f.type or target not in _COLLAPSE_TARGETS:
            # already a core type, or maps to a non-core INSDC feature that Phase 2
            # handles during generation (ncRNA/precursor_RNA/5'UTR/...) -> leave as-is
            continue
        old = f.type
        f.type = target
        changes.append(Change("rename-type", f.id or "?", f"{old} -> {target}"))
        for qual in vocab.feature_qualifiers.get(old, ()):
            if _is_placeholder(qual):
                changes.append(Change("needs-manual", f.id or "?",
                                      f"qualifier {qual} for {old} needs a manual value (not added)"))
                continue
            key, val = _qualifier_to_attr(qual)
            if key in f.attributes:
                continue
            f.attributes[key] = [val if val is not None else "true"]
            changes.append(Change("add-qualifier", f.id or "?",
                                  f"added {key}={f.attributes[key][0]}"))
    return changes
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_pass_so_terms.py -v` → 10 passed。全体 `docker exec ddbj-gff-dev uv run pytest -q`（回帰なし）。

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/normalize/passes.py tests/test_normalize_pass_so_terms.py
git commit -m "fix(normalize): collapse SO terms only to core whitelist types (Finding A/B)"
```

---

### Task 2: 受け入れ統合テスト（Finding A 回帰防止）

**Files:**
- Test: `tests/test_normalize_integration.py`（テスト関数を1件追記）

**Interfaces:**
- Consumes: `parse`/`validate`/`normalize`/`NormalizeConfig`（同ファイルで既に import 済み）。

- [ ] **Step 1: 受け入れテストを追記（失敗するテスト）**

`tests/test_normalize_integration.py` の末尾に追記:
```python
def test_finding_a_no_feature_type_warning_for_left_types():
    # coding_exon/pseudogenic_CDS collapse to core types; miRNA/pre_miRNA/UTR are left for
    # Phase 2. After normalize, the validator must NOT emit feature-type-not-insdc (Finding A).
    gff = (
        "##gff-version 3\n"
        "##sequence-region chr1 1 10000\n"
        "chr1\tS\tgene\t100\t900\t.\t+\t.\tID=g1;locus_tag=ABC_1\n"
        "chr1\tS\tmRNA\t100\t900\t.\t+\t.\tID=m1;Parent=g1\n"
        "chr1\tS\tfive_prime_UTR\t100\t129\t.\t+\t.\tID=u1;Parent=m1\n"
        "chr1\tS\tcoding_exon\t130\t900\t.\t+\t.\tID=e1;Parent=m1\n"
        "chr1\tS\tpseudogenic_CDS\t130\t870\t.\t+\t0\tID=c1;Parent=m1\n"
        "chr1\tS\tmiRNA\t2000\t2100\t.\t+\t.\tID=r1\n"
        "chr1\tS\tpre_miRNA\t3000\t3200\t.\t+\t.\tID=r2\n"
    )
    norm, _ = normalize(parse(gff), seq_lengths={"chr1": 10000}, config=NormalizeConfig(taxid=3702))
    types = {f.type for f in norm.features}
    assert "exon" in types and "CDS" in types                 # coding_exon/pseudogenic_CDS collapsed
    assert "coding_exon" not in types and "pseudogenic_CDS" not in types
    assert "miRNA" in types and "pre_miRNA" in types           # RNA leaves left for Phase 2
    assert "five_prime_UTR" in types                           # UTR left alone
    codes = {d.code for d in validate(norm)}
    assert "feature-type-not-insdc" not in codes               # Finding A regression guard
```

- [ ] **Step 2: 失敗確認 → 成功確認**
Task 1 実装後はこのテストは通る想定。実行: `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_integration.py -v`（非slow が pass）。もし Task 1 前に書いた場合は FAIL（feature-type-not-insdc が出る）。

- [ ] **Step 3: 全体確認** — `docker exec ddbj-gff-dev uv run pytest -q`（slow 除外）→ 全 pass。

- [ ] **Step 4: Commit**
```bash
git add tests/test_normalize_integration.py
git commit -m "test(normalize): acceptance - no feature-type-not-insdc after whitelist collapse"
```

---

## Self-Review

**1. Spec coverage**（spec §→タスク）:
- §2 F-D2/F-D3（whitelist collapse・SO→INSDC は Phase2） → Task 1（`_COLLAPSE_TARGETS` ＋ 判定）。
- §3 型ごと挙動表 → Task 1 ユニットテスト（collapse: pseudogenic_CDS/processed_pseudogene/coding_exon/spliceosomal_intron; 放置: miRNA/pre_miRNA/five_prime_UTR/binding_site/mobile_genetic_element; 未マップ; no-op; clobber）。
- §4 受け入れ（validate に feature-type-not-insdc が出ない） → Task 2。
- §2 F-D4（qualifier ロジック維持） → Task 1 で `_is_placeholder`/`_qualifier_to_attr` 不変・clobber テスト維持。
- §4 受け入れ基準2（ゼニゴケ実走再検証） → マージ後にコントローラが手動実行（プラン外の検証ステップ）。

ギャップ/意図的事項:
- `needs-manual` 経路は whitelisted target に実データ上プレースホルダ qualifier が無いため**実質不発火**になり、専用テストを置かない（防御コードは残す）。レビューには「dead-in-practice な防御ガード」と明示。

**2. Placeholder scan**: 各ステップに実コード（全置換テスト・全置換関数・フィクスチャ文字列）。"TBD" 等なし。GFF フィクスチャは TAB 区切り（`\t`）。

**3. Type consistency**: `_COLLAPSE_TARGETS`（Task1）/`pass_so_terms(doc, ctx)->list[Change]`（不変シグネチャ）/`Change(action,target,message)`/`normalize(doc,*,seq_lengths,config)->(doc,report)`/`validate(doc)->list[Diagnostic]`(`.code`) — 全タスク・既存コードと一致。`Feature(id,source,type,spans,attributes,parent_ids)` の位置引数順も既存テストに一致。
