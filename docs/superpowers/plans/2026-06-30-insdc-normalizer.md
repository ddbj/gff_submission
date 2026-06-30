# INSDC 正規化器（フェーズ3-B・共通ケース MVP）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 任意の GFF3 を INSDC GFF3 形式へ近づける GFF→GFF 正規化器サブパッケージ `ddbj_gff.normalize` を、ディレクティブ補完 ＋ SO-term 正規化（共通ケース MVP）として実装する。

**Architecture:** 案A: パス・レジストリ。`normalize(doc, *, seq_lengths=None, config=None) -> (GffDocument, NormalizationReport)` が、入力 doc の deepcopy 作業コピーに独立した正規化パス関数群（`ALL_PASSES`）を順次適用し、変更を `Change` として集約。3-A `ddbj_gff.validate` の `ALL_RULES` と対称。Phase1 のパーサ/ライター/モデル、3-A の vocab を再利用。

**Tech Stack:** Python 3.11+ / 既存 `ddbj_gff`（Phase1）・`ddbj_gff.validate`（Phase3-A） / pytest / dev コンテナ。`normalize` 本体は stdlib のみ（FASTA 読込のみ CLI 側で Biopython）。

## Global Constraints

- Python `>=3.11`、各モジュール冒頭 `from __future__ import annotations`。実行依存は biopython のみ（FASTA 読込に使用、CLI 限定）。`normalize`/`passes`/`report`/`config` 本体は stdlib のみ。
- 変更は `src/ddbj_gff/normalize/` 配下に新規作成。例外: Task 1 のみ既存 `src/ddbj_gff/validate/vocab.py` を拡張（後方互換を保つ）。
- **detect-only ではなく transform**: `normalize()` は入力 `doc` を**変更しない**（`copy.deepcopy` した作業コピーを変換）。冪等（正規化済みを再度 normalize しても重複追加しない）。
- ディレクティブ追加時は `Directive(raw, kind, value)` の **`raw` と `value` を両方**正しく設定する（writer は `raw` を出力、`GffDocument` のプロパティと後段 `validate()` は `value` を読む）。
- 値の形（Phase1 準拠）: `gff-version`/`insdc-gff-version` → `str`、`species` → `int`(taxid)、`sequence-region` → `(seqid, start, end)` タプル、`transl_table` → `dict {"primary": N}`。
- qualifier 自動付与は具体値のみ。プレースホルダ（`<` `>` を含む）・末尾/途中 `*` を含む値は付与せず report に記録。
- **テストは dev コンテナ内**: `docker exec ddbj-gff-dev uv run pytest …`。`git` は host。slow テストは `-m slow`。
- 各タスクは「失敗テスト→失敗確認→最小実装→成功確認→コミット」。

## Phase1/3-A 既存インターフェース（参照）

- `ddbj_gff.parse(text) -> GffDocument`。`ddbj_gff.writer.write(doc, *, canonical_sort=False) -> str`（ディレクティブは `d.raw` を、その後 features を出力。`doc.fasta` があれば `##FASTA` も）。
- `model.Directive(raw: str, kind: str, value=None)`（dataclass）。`model.Feature(id, source, type, spans, attributes, parent_ids, children, parents)`、`attributes: dict[str, list[str]]`、プロパティ `.transl_table -> int|None`（`int(attributes["transl_table"][0])`）。`model.GffDocument(directives, features, feature_index, roots, fasta, sequences, diagnostics)`。
- `GffDocument` プロパティ（全て directives を読む）: `gff_version -> str|None`、`insdc_gff_version -> str|None`、`species -> int|None`（taxid。URL 中 `id=(\d+)` を int 化）、`sequence_regions -> dict[str,(int,int)]`、`transl_table_map -> dict|None`。
- `ddbj_gff.validate.validate(doc, *, severity_overrides=None) -> list[Diagnostic]`、`ddbj_gff.errors.Severity`（`.ERROR` 等、`.value` は文字列）。
- `ddbj_gff.validate.vocab.load_vocab() -> Vocab`（lru_cache）。現状 `Vocab(feature_types: frozenset, insdc_map: dict[str,str], dbxref_dbtags: frozenset)`（frozen dataclass）。`insdc_map` は `feature-mapping.tsv` 列2(SO term)→列4(INSDC Feature)。

## File Structure

| ファイル | 責務 |
|---|---|
| `src/ddbj_gff/validate/vocab.py`（変更） | `Vocab` に `feature_qualifiers` 追加、列5–6 ロード＋重複 SO-term dedup |
| `src/ddbj_gff/normalize/__init__.py` | `normalize` / `NormalizationReport` / `NormalizeConfig` 再エクスポート |
| `src/ddbj_gff/normalize/report.py` | `Change` / `NormalizationReport`（dataclass・render） |
| `src/ddbj_gff/normalize/config.py` | `NormalizeConfig` ＋ `load_normalize_config`（TOML） |
| `src/ddbj_gff/normalize/passes.py` | `NormalizeContext` / `pass_directives` / `pass_so_terms` |
| `src/ddbj_gff/normalize/normalize.py` | `normalize(...)` ＋ `ALL_PASSES` |
| `src/ddbj_gff/normalize/cli.py`, `__main__.py` | CLI |
| `tests/test_normalize_*.py`, `tests/normalize_fixtures/` | テスト・フィクスチャ |

---

### Task 1: vocab 拡張（feature_qualifiers ＋ 重複 dedup）

**Files:**
- Modify: `src/ddbj_gff/validate/vocab.py`
- Test: `tests/test_normalize_vocab.py`

**Interfaces:**
- Consumes: 同梱 `data/feature-mapping.tsv`（列2 SO term / 列4 INSDC Feature / 列5–6 Qualifier）。
- Produces: `Vocab.feature_qualifiers: dict[str, tuple[str, ...]]`（SO term → 生 qualifier 文字列群）。既存 `feature_types`/`insdc_map`/`dbxref_dbtags` は不変。`load_vocab()` シグネチャ不変。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_normalize_vocab.py`:
```python
from ddbj_gff.validate.vocab import load_vocab


def test_feature_qualifiers_loaded():
    v = load_vocab()
    # pseudogenic_CDS -> CDS, qualifier /pseudo
    assert v.insdc_map["pseudogenic_CDS"] == "CDS"
    assert any("pseudo" in q for q in v.feature_qualifiers["pseudogenic_CDS"])


def test_duplicate_so_term_prefers_concrete():
    v = load_vocab()
    # LINE_element appears twice: /mobile_element_type="LINE*" and ="LINE"; concrete (no '*') wins
    quals = v.feature_qualifiers["LINE_element"]
    assert quals
    assert all("*" not in q for q in quals)


def test_insdc_map_unchanged_for_3a():
    v = load_vocab()
    assert "CDS" in v.feature_types
    assert v.insdc_map.get("ncRNA_gene") == "ncRNA"
    # term with no qualifier maps to empty tuple
    assert v.feature_qualifiers.get("CDS", ()) == ()
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_vocab.py -v` → FAIL（`AttributeError: feature_qualifiers`）

- [ ] **Step 3: 実装** — `src/ddbj_gff/validate/vocab.py` を更新:
```python
from __future__ import annotations

import functools
from dataclasses import dataclass, field
from pathlib import Path

_DATA = Path(__file__).parent / "data"


@dataclass(frozen=True)
class Vocab:
    feature_types: frozenset[str]
    insdc_map: dict[str, str]
    dbxref_dbtags: frozenset[str]
    feature_qualifiers: dict[str, tuple[str, ...]] = field(default_factory=dict)


def _is_concrete(quals: tuple[str, ...]) -> bool:
    return all("<" not in q and ">" not in q and "*" not in q for q in quals)


def _read_feature_mapping() -> tuple[frozenset[str], dict[str, str], dict[str, tuple[str, ...]]]:
    terms: set[str] = set()
    mapping: dict[str, str] = {}
    qualifiers: dict[str, tuple[str, ...]] = {}
    with open(_DATA / "feature-mapping.tsv", encoding="utf-8") as fh:
        next(fh, None)  # header
        for line in fh:
            if not line.strip():
                continue
            cols = [c.strip() for c in line.rstrip("\n").split("\t")]
            if len(cols) < 2 or not cols[1]:
                continue
            so_term = cols[1]
            terms.add(so_term)
            insdc = cols[3] if len(cols) > 3 and cols[3] else ""
            quals = tuple(c for c in cols[4:6] if c)  # 列5,6 (Qualifier 1,2)
            if not insdc:
                continue
            if so_term not in mapping:
                mapping[so_term] = insdc
                qualifiers[so_term] = quals
            elif not _is_concrete(qualifiers.get(so_term, ())) and _is_concrete(quals):
                qualifiers[so_term] = quals  # 重複: 具体値の行を優先
    return frozenset(terms), mapping, qualifiers


def _read_dbxref() -> frozenset[str]:
    tags: set[str] = set()
    with open(_DATA / "dbxref.tsv", encoding="utf-8") as fh:
        for line in fh:
            t = line.strip()
            if t and not t.startswith("#"):
                tags.add(t)
    return frozenset(tags)


@functools.lru_cache(maxsize=1)
def load_vocab() -> Vocab:
    terms, mapping, quals = _read_feature_mapping()
    return Vocab(
        feature_types=terms,
        insdc_map=mapping,
        dbxref_dbtags=_read_dbxref(),
        feature_qualifiers=quals,
    )
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_vocab.py -v` → 3 passed。回帰確認 `docker exec ddbj-gff-dev uv run pytest tests/test_validate_vocab.py -q` → pass。

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/validate/vocab.py tests/test_normalize_vocab.py
git commit -m "feat(vocab): load feature qualifiers + dedup duplicate SO terms"
```

---

### Task 2: report.py（Change / NormalizationReport）＋ パッケージ init

**Files:**
- Create: `src/ddbj_gff/normalize/__init__.py`, `src/ddbj_gff/normalize/report.py`
- Test: `tests/test_normalize_report.py`

**Interfaces:**
- Produces: `Change(action: str, target: str, message: str)`（frozen dataclass）。`NormalizationReport(applied: list[Change], unresolved: list[Change])` ＋ `.render() -> str`。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_normalize_report.py`:
```python
from ddbj_gff.normalize.report import Change, NormalizationReport


def test_report_render_and_counts():
    r = NormalizationReport(
        applied=[Change("rename-type", "c1", "pseudogenic_CDS -> CDS")],
        unresolved=[Change("no-taxid", "species", "no taxid provided")],
    )
    text = r.render()
    assert "1 applied, 1 need attention" in text
    assert "rename-type" in text
    assert "no-taxid" in text


def test_report_defaults_empty():
    r = NormalizationReport()
    assert r.applied == [] and r.unresolved == []
    assert "0 applied, 0 need attention" in r.render()
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_report.py -v` → FAIL（ModuleNotFound）

- [ ] **Step 3: 実装**

`src/ddbj_gff/normalize/__init__.py`:
```python
"""ddbj_gff.normalize: INSDC GFF3 normalizer (phase 3B, common-case MVP)."""

__all__ = []
```

`src/ddbj_gff/normalize/report.py`:
```python
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Change:
    action: str    # add-directive | rename-type | add-qualifier | approx-region
                   # | unmapped-type | needs-manual | no-taxid
    target: str    # feature id / seqid / directive 名
    message: str


@dataclass
class NormalizationReport:
    applied: list = field(default_factory=list)
    unresolved: list = field(default_factory=list)

    def render(self) -> str:
        lines = [f"normalization: {len(self.applied)} applied, {len(self.unresolved)} need attention"]
        for c in self.applied:
            lines.append(f"  [applied]   {c.action} {c.target}: {c.message}")
        for c in self.unresolved:
            lines.append(f"  [attention] {c.action} {c.target}: {c.message}")
        return "\n".join(lines) + "\n"
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_report.py -v` → 2 passed

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/normalize/__init__.py src/ddbj_gff/normalize/report.py tests/test_normalize_report.py
git commit -m "feat(normalize): Change/NormalizationReport model"
```

---

### Task 3: config.py（NormalizeConfig ＋ TOML ロード）

**Files:**
- Create: `src/ddbj_gff/normalize/config.py`
- Test: `tests/test_normalize_config.py`

**Interfaces:**
- Produces: `NormalizeConfig(taxid: int|None=None, transl_table: int=1, insdc_gff_version: str="1.0.0")`（非frozen dataclass）。`load_normalize_config(path: str) -> NormalizeConfig`（TOML `[normalize]` を読む）。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_normalize_config.py`:
```python
from ddbj_gff.normalize.config import NormalizeConfig, load_normalize_config


def test_defaults():
    c = NormalizeConfig()
    assert c.taxid is None
    assert c.transl_table == 1
    assert c.insdc_gff_version == "1.0.0"


def test_load_from_toml(tmp_path):
    p = tmp_path / "n.toml"
    p.write_text('[normalize]\ntaxid = 3702\ntransl_table = 11\ninsdc_gff_version = "1.0.0"\n')
    c = load_normalize_config(str(p))
    assert c.taxid == 3702
    assert c.transl_table == 11


def test_load_missing_section_uses_defaults(tmp_path):
    p = tmp_path / "n.toml"
    p.write_text('[other]\nx = 1\n')
    c = load_normalize_config(str(p))
    assert c.taxid is None and c.transl_table == 1
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_config.py -v` → FAIL

- [ ] **Step 3: 実装**

`src/ddbj_gff/normalize/config.py`:
```python
from __future__ import annotations

import tomllib
from dataclasses import dataclass


@dataclass
class NormalizeConfig:
    taxid: int | None = None
    transl_table: int = 1
    insdc_gff_version: str = "1.0.0"


def load_normalize_config(path: str) -> NormalizeConfig:
    with open(path, "rb") as fh:
        data = tomllib.load(fh)
    n = data.get("normalize", {})
    return NormalizeConfig(
        taxid=n.get("taxid"),
        transl_table=n.get("transl_table", 1),
        insdc_gff_version=n.get("insdc_gff_version", "1.0.0"),
    )
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_config.py -v` → 3 passed

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/normalize/config.py tests/test_normalize_config.py
git commit -m "feat(normalize): NormalizeConfig + TOML loader"
```

---

### Task 4: pass_directives（ディレクティブ補完）＋ NormalizeContext

**Files:**
- Create: `src/ddbj_gff/normalize/passes.py`
- Test: `tests/test_normalize_pass_directives.py`

**Interfaces:**
- Consumes: `report.Change`、`config.NormalizeConfig`、`model.Directive`、`GffDocument` プロパティ。
- Produces: `NormalizeContext(vocab, seq_lengths, config)`（dataclass）。`pass_directives(doc, ctx) -> list[Change]`（doc を in-place 変換）。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_normalize_pass_directives.py`:
```python
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
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_pass_directives.py -v` → FAIL

- [ ] **Step 3: 実装**

`src/ddbj_gff/normalize/passes.py`:
```python
from __future__ import annotations

from dataclasses import dataclass

from ..model import Directive
from .report import Change

_SPECIES_URL = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}"


@dataclass
class NormalizeContext:
    vocab: object             # validate.vocab.Vocab
    seq_lengths: dict | None
    config: object            # NormalizeConfig


def pass_directives(doc, ctx) -> list:
    cfg = ctx.config
    changes: list = []

    if doc.gff_version is None:
        doc.directives.insert(0, Directive("##gff-version 3", "gff-version", "3"))
        changes.append(Change("add-directive", "gff-version", "added ##gff-version 3"))

    if doc.insdc_gff_version is None:
        v = cfg.insdc_gff_version
        doc.directives.append(Directive(f"#!insdc-gff-version {v}", "insdc-gff-version", v))
        changes.append(Change("add-directive", "insdc-gff-version", f"added #!insdc-gff-version {v}"))

    if not isinstance(doc.species, int):
        if cfg.taxid is not None:
            url = _SPECIES_URL.format(taxid=cfg.taxid)
            doc.directives.append(Directive(f"##species {url}", "species", cfg.taxid))
            changes.append(Change("add-directive", "species", f"added ##species (taxid {cfg.taxid})"))
        else:
            changes.append(Change("no-taxid", "species",
                                  "##species not added: no taxid (set [normalize].taxid or --taxid)"))

    have = set(doc.sequence_regions)
    seqids: list = []
    for f in doc.features:
        for s in f.spans:
            if s.seqid not in have and s.seqid not in seqids:
                seqids.append(s.seqid)
    for seqid in seqids:
        length = ctx.seq_lengths.get(seqid) if ctx.seq_lengths else None
        approx = length is None
        if approx:
            length = max((s.end for f in doc.features for s in f.spans if s.seqid == seqid), default=1)
        doc.directives.append(
            Directive(f"##sequence-region {seqid} 1 {length}", "sequence-region", (seqid, 1, length)))
        if approx:
            changes.append(Change("approx-region", seqid,
                                  f"added ##sequence-region {seqid} 1 {length} "
                                  f"(length approximated from max feature end; provide --fasta for true length)"))
        else:
            changes.append(Change("add-directive", seqid, f"added ##sequence-region {seqid} 1 {length}"))

    if doc.transl_table_map is None:
        if any(f.type == "CDS" and f.transl_table is None for f in doc.features):
            vals = {f.transl_table for f in doc.features if f.type == "CDS" and f.transl_table is not None}
            n = vals.pop() if len(vals) == 1 else cfg.transl_table
            doc.directives.append(Directive(f"#!transl_table primary:{n}", "transl_table", {"primary": n}))
            changes.append(Change("add-directive", "transl_table", f"added #!transl_table primary:{n}"))

    return changes
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_pass_directives.py -v` → 7 passed

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/normalize/passes.py tests/test_normalize_pass_directives.py
git commit -m "feat(normalize): pass_directives + NormalizeContext"
```

---

### Task 5: pass_so_terms（SO-term 正規化）

**Files:**
- Modify: `src/ddbj_gff/normalize/passes.py`
- Test: `tests/test_normalize_pass_so_terms.py`

**Interfaces:**
- Consumes: `Vocab.insdc_map` / `Vocab.feature_qualifiers`（Task 1）、`report.Change`、`NormalizeContext`。
- Produces: `pass_so_terms(doc, ctx) -> list[Change]`（feature.type を INSDC 名に書換＋具体 qualifier を属性付与）。ヘルパ `_qualifier_to_attr`、`_is_placeholder`。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_normalize_pass_so_terms.py`:
```python
from ddbj_gff.model import Feature, Span, GffDocument
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.normalize.passes import NormalizeContext, pass_so_terms
from ddbj_gff.validate.vocab import load_vocab


def _ctx():
    return NormalizeContext(vocab=load_vocab(), seq_lengths=None, config=NormalizeConfig())


def _doc(*feats):
    return GffDocument(features=list(feats))


def test_rename_pseudogenic_cds_adds_pseudo_flag():
    f = Feature("c", "S", "pseudogenic_CDS", [Span("chr1", 1, 9, "+", 0)], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "CDS"
    assert f.attributes.get("pseudo") == ["true"]
    assert any(c.action == "rename-type" for c in changes)
    assert any(c.action == "add-qualifier" for c in changes)


def test_rename_processed_pseudogene_adds_keyed_qualifier():
    f = Feature("g", "S", "processed_pseudogene", [Span("chr1", 1, 9, "+")], {}, [])
    pass_so_terms(_doc(f), _ctx())
    assert f.type == "gene"
    assert f.attributes.get("pseudogene") == ["processed"]


def test_same_name_is_noop():
    f = Feature("c", "S", "CDS", [Span("chr1", 1, 9, "+", 0)], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "CDS"
    assert changes == []


def test_unmapped_type_reported_unchanged():
    f = Feature("x", "S", "totally_made_up_type", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "totally_made_up_type"
    assert any(c.action == "unmapped-type" for c in changes)


def test_placeholder_qualifier_not_fabricated():
    # mobile_genetic_element -> mobile_element with /mobile_element_type="other:<NAME>" (placeholder)
    f = Feature("m", "S", "mobile_genetic_element", [Span("chr1", 1, 9, "+")], {}, [])
    changes = pass_so_terms(_doc(f), _ctx())
    assert f.type == "mobile_element"
    assert "mobile_element_type" not in f.attributes  # placeholder value NOT added
    assert any(c.action == "needs-manual" for c in changes)


def test_existing_attribute_not_clobbered():
    f = Feature("c", "S", "pseudogenic_CDS", [Span("chr1", 1, 9, "+", 0)], {"pseudo": ["existing"]}, [])
    pass_so_terms(_doc(f), _ctx())
    assert f.attributes["pseudo"] == ["existing"]
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_pass_so_terms.py -v` → FAIL

- [ ] **Step 3: 実装** — `src/ddbj_gff/normalize/passes.py` 末尾に追記:
```python
def _is_placeholder(qual: str) -> bool:
    return "<" in qual or ">" in qual or "*" in qual


def _qualifier_to_attr(qual: str) -> tuple[str, str | None]:
    body = qual.lstrip("/")
    if "=" in body:
        key, val = body.split("=", 1)
        return key.strip(), val.strip().strip('"')
    return body.strip(), None  # valueless flag (e.g. /pseudo)


def pass_so_terms(doc, ctx) -> list:
    vocab = ctx.vocab
    changes: list = []
    for f in doc.features:
        target = vocab.insdc_map.get(f.type)
        if target is None:
            changes.append(Change("unmapped-type", f.id or "?",
                                  f"feature type {f.type!r} is not a known SO term; left unchanged"))
            continue
        if target == f.type:
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
                continue  # don't clobber existing
            f.attributes[key] = [val if val is not None else "true"]
            changes.append(Change("add-qualifier", f.id or "?",
                                  f"added {key}={f.attributes[key][0]}"))
    return changes
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_pass_so_terms.py -v` → 6 passed

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/normalize/passes.py tests/test_normalize_pass_so_terms.py
git commit -m "feat(normalize): pass_so_terms (SO term -> INSDC + qualifiers)"
```

---

### Task 6: normalize() 集約 ＋ ALL_PASSES ＋ 再エクスポート

**Files:**
- Create: `src/ddbj_gff/normalize/normalize.py`
- Modify: `src/ddbj_gff/normalize/__init__.py`
- Test: `tests/test_normalize.py`

**Interfaces:**
- Consumes: `passes.ALL_PASSES` 相当（`pass_directives`/`pass_so_terms`）、`NormalizeContext`、`report.NormalizationReport`、`config.NormalizeConfig`、`ddbj_gff.validate.vocab.load_vocab`。
- Produces: `normalize(doc, *, seq_lengths=None, config=None) -> tuple[GffDocument, NormalizationReport]`。`ALL_PASSES` リスト。`__init__` から `normalize`/`NormalizationReport`/`NormalizeConfig` 再エクスポート。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_normalize.py`:
```python
from ddbj_gff import parse
from ddbj_gff.errors import Severity
from ddbj_gff.validate import validate
from ddbj_gff.normalize import normalize, NormalizeConfig

GFF_MESSY = (
    "##gff-version 3\n"
    "chr1\tS\tgene\t100\t900\t.\t+\t.\tID=g;locus_tag=ABC_1\n"
    "chr1\tS\tCDS\t130\t870\t.\t+\t0\tID=c;Parent=g\n"
)  # missing insdc-gff-version, species, sequence-region, transl_table


def test_normalize_clears_targeted_errors():
    norm, report = normalize(parse(GFF_MESSY),
                             seq_lengths={"chr1": 10000}, config=NormalizeConfig(taxid=3702))
    codes = {d.code for d in validate(norm) if d.severity == Severity.ERROR}
    assert "missing-insdc-gff-version" not in codes
    assert "missing-species-taxid" not in codes
    assert "missing-sequence-region" not in codes
    assert "cds-missing-transl-table" not in codes


def test_normalize_does_not_mutate_input():
    doc = parse(GFF_MESSY)
    before = len(doc.directives)
    normalize(doc, config=NormalizeConfig(taxid=3702))
    assert len(doc.directives) == before  # input untouched (works on a copy)


def test_report_separates_applied_and_unresolved():
    # no taxid -> species unresolved; an unmapped feature -> unresolved
    g = "##gff-version 3\nchr1\tS\tmade_up\t1\t9\t.\t+\t.\tID=x\n"
    _, report = normalize(parse(g))
    assert any(c.action == "no-taxid" for c in report.unresolved)
    assert any(c.action == "unmapped-type" for c in report.unresolved)
    assert all(c.action in ("add-directive",) for c in report.applied)


def test_idempotent():
    norm1, _ = normalize(parse(GFF_MESSY), seq_lengths={"chr1": 10000}, config=NormalizeConfig(taxid=3702))
    norm2, report2 = normalize(norm1, seq_lengths={"chr1": 10000}, config=NormalizeConfig(taxid=3702))
    assert len(norm2.directives) == len(norm1.directives)
    assert report2.applied == []  # nothing left to change
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize.py -v` → FAIL（ImportError: normalize）

- [ ] **Step 3: 実装**

`src/ddbj_gff/normalize/normalize.py`:
```python
from __future__ import annotations

import copy

from ..validate.vocab import load_vocab
from .config import NormalizeConfig
from .passes import NormalizeContext, pass_directives, pass_so_terms
from .report import NormalizationReport

ALL_PASSES = [pass_directives, pass_so_terms]

# actions that represent a clean applied change; everything else needs human attention
_APPLIED = {"add-directive", "rename-type", "add-qualifier"}


def normalize(doc, *, seq_lengths=None, config=None) -> tuple:
    config = config or NormalizeConfig()
    work = copy.deepcopy(doc)
    ctx = NormalizeContext(vocab=load_vocab(), seq_lengths=seq_lengths, config=config)
    applied: list = []
    unresolved: list = []
    for run_pass in ALL_PASSES:
        for change in run_pass(work, ctx):
            (applied if change.action in _APPLIED else unresolved).append(change)
    return work, NormalizationReport(applied=applied, unresolved=unresolved)
```

`src/ddbj_gff/normalize/__init__.py` を更新:
```python
"""ddbj_gff.normalize: INSDC GFF3 normalizer (phase 3B, common-case MVP)."""

from .config import NormalizeConfig
from .normalize import normalize
from .report import NormalizationReport

__all__ = ["normalize", "NormalizationReport", "NormalizeConfig"]
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize.py -v` → 4 passed。全体 `docker exec ddbj-gff-dev uv run pytest -q`。

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/normalize/normalize.py src/ddbj_gff/normalize/__init__.py tests/test_normalize.py
git commit -m "feat(normalize): normalize() aggregation + round-trip oracle"
```

---

### Task 7: CLI

**Files:**
- Create: `src/ddbj_gff/normalize/cli.py`, `src/ddbj_gff/normalize/__main__.py`
- Test: `tests/test_normalize_cli.py`

**Interfaces:**
- Consumes: `parse`（Phase1）、`writer.write`、`normalize`（Task 6）、`NormalizeConfig`/`load_normalize_config`（Task 3）。
- Produces: `cli.main(argv=None) -> int`（`--gff` 必須、`--fasta`/`--config`/`--taxid`/`--transl-table`/`--insdc-gff-version`/`--out`/`--report` 任意）。`python -m ddbj_gff.normalize` 実行可。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_normalize_cli.py`:
```python
import pytest
from ddbj_gff.normalize.cli import main

GFF = "##gff-version 3\nchr1\tS\tgene\t1\t9\t.\t+\t.\tID=g;locus_tag=X_1\n"


def test_cli_normalizes_to_stdout(tmp_path, capsys):
    p = tmp_path / "in.gff"
    p.write_text(GFF)
    rc = main(["--gff", str(p), "--taxid", "3702"])
    assert rc == 0
    out = capsys.readouterr()
    assert "#!insdc-gff-version 1.0.0" in out.out
    assert "##species" in out.out and "id=3702" in out.out
    assert "##sequence-region chr1 1" in out.out
    assert "[applied]" in out.err  # report to stderr


def test_cli_writes_out_file(tmp_path):
    p = tmp_path / "in.gff"
    p.write_text(GFF)
    outp = tmp_path / "out.gff"
    rc = main(["--gff", str(p), "--taxid", "3702", "--out", str(outp)])
    assert rc == 0
    assert "#!insdc-gff-version" in outp.read_text()


def test_cli_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        main(["--gff", str(tmp_path / "nope.gff")])
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_cli.py -v` → FAIL

- [ ] **Step 3: 実装**

`src/ddbj_gff/normalize/cli.py`:
```python
from __future__ import annotations

import argparse
import sys

from .. import parse
from ..writer import write
from .config import NormalizeConfig, load_normalize_config
from .normalize import normalize


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="ddbj_gff.normalize",
                                 description="Normalize a GFF3 toward the INSDC profile")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta")
    ap.add_argument("--config")
    ap.add_argument("--taxid", type=int)
    ap.add_argument("--transl-table", type=int, dest="transl_table")
    ap.add_argument("--insdc-gff-version", dest="insdc_gff_version")
    ap.add_argument("--out")
    ap.add_argument("--report")
    args = ap.parse_args(argv)

    cfg = load_normalize_config(args.config) if args.config else NormalizeConfig()
    if args.taxid is not None:
        cfg.taxid = args.taxid
    if args.transl_table is not None:
        cfg.transl_table = args.transl_table
    if args.insdc_gff_version is not None:
        cfg.insdc_gff_version = args.insdc_gff_version

    seq_lengths = None
    if args.fasta:
        from Bio import SeqIO
        seq_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(args.fasta, "fasta")}

    with open(args.gff, encoding="ascii", errors="replace") as fh:
        doc = parse(fh.read())

    norm, report = normalize(doc, seq_lengths=seq_lengths, config=cfg)

    out_text = write(norm)
    if args.out:
        with open(args.out, "w", encoding="ascii") as fh:
            fh.write(out_text)
    else:
        sys.stdout.write(out_text)

    report_text = report.render()
    if args.report:
        with open(args.report, "w", encoding="ascii") as fh:
            fh.write(report_text)
    else:
        sys.stderr.write(report_text)

    return 0
```

`src/ddbj_gff/normalize/__main__.py`:
```python
import sys

from .cli import main

sys.exit(main())
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_normalize_cli.py -v` → 3 passed。全体 `docker exec ddbj-gff-dev uv run pytest -q`。

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/normalize/cli.py src/ddbj_gff/normalize/__main__.py tests/test_normalize_cli.py
git commit -m "feat(normalize): CLI entry point"
```

---

### Task 8: フィクスチャ ＋ 実ファイル統合テスト（slow）

**Files:**
- Create: `tests/normalize_fixtures/messy_input.gff3`
- Test: `tests/test_normalize_integration.py`

**Interfaces:** Consumes `parse`/`normalize`/`validate`/`writer.write`、実 example（存在時のみ）。

- [ ] **Step 1: フィクスチャ＋テストを書く**

`tests/normalize_fixtures/messy_input.gff3`（TAB 区切り。必須ディレクティブ欠落＋要マッピング SO-term）:
```text
##gff-version 3
chr1	S	gene	100	900	.	+	.	ID=g1;locus_tag=ABC_000010
chr1	S	coding_exon	100	900	.	+	.	ID=e1;Parent=g1
chr1	S	pseudogenic_CDS	130	870	.	+	0	ID=c1;Parent=g1
```

`tests/test_normalize_integration.py`:
```python
from pathlib import Path
import pytest
from ddbj_gff import parse
from ddbj_gff.errors import Severity
from ddbj_gff.validate import validate
from ddbj_gff.writer import write
from ddbj_gff.normalize import normalize, NormalizeConfig

FIX = Path(__file__).parent / "normalize_fixtures"
ROOT = Path(__file__).resolve().parents[1]


def test_messy_fixture_normalizes_and_validates():
    doc = parse((FIX / "messy_input.gff3").read_text())
    norm, report = normalize(doc, seq_lengths={"chr1": 5000}, config=NormalizeConfig(taxid=3702))
    # SO terms normalized to INSDC names
    types = {f.type for f in norm.features}
    assert "exon" in types and "CDS" in types          # coding_exon->exon, pseudogenic_CDS->CDS
    assert "coding_exon" not in types and "pseudogenic_CDS" not in types
    # targeted ERRORs cleared by directive completion
    errors = {d.code for d in validate(norm) if d.severity == Severity.ERROR}
    assert "missing-insdc-gff-version" not in errors
    assert "missing-species-taxid" not in errors
    assert "missing-sequence-region" not in errors
    # output is writable GFF text and round-trips through parse
    text = write(norm)
    assert "#!insdc-gff-version 1.0.0" in text
    reparsed = parse(text)
    assert reparsed.insdc_gff_version == "1.0.0"


@pytest.mark.slow
def test_rice_cp_normalize_clears_version_keeps_special_case():
    p = ROOT / "examples" / "rice_cp" / "rice_cp.gff3"
    if not p.exists():
        pytest.skip(f"missing {p}")
    norm, _ = normalize(parse(p.read_text(errors="replace")), config=NormalizeConfig(taxid=39947))
    diags = validate(norm)
    errors = {d.code for d in diags if d.severity == Severity.ERROR}
    codes = {d.code for d in diags}
    assert "missing-insdc-gff-version" not in errors          # 3B added the directive
    assert "noncanonical-special-case" in codes               # trans-splicing still flagged (deferred to 3B-full)
```

- [ ] **Step 2: 実行**
`docker exec ddbj-gff-dev uv run pytest tests/test_normalize_integration.py -v`（非slow の messy フィクスチャは必ず実行・pass）。
`docker exec ddbj-gff-dev uv run pytest tests/test_normalize_integration.py -m slow -v`（rice_cp 実ファイル）。
もし messy フィクスチャで想定外挙動が出たら、フィクスチャ/期待値を調べる（パス実装は変えない）。

- [ ] **Step 3: 全体確認** — `docker exec ddbj-gff-dev uv run pytest -q`（slow 除外）→ 全 pass。

- [ ] **Step 4: Commit**
```bash
git add tests/normalize_fixtures/messy_input.gff3 tests/test_normalize_integration.py
git commit -m "test(normalize): fixture + slow integration (normalize->validate oracle)"
```

---

## Self-Review

**1. Spec coverage**（spec §→タスク）:
- §3 構成・vocab 拡張（feature_qualifiers・dedup） → Task 1。`ALL_PASSES` → Task 6。
- §4 pass_directives（gff-version/insdc-gff-version/species/sequence-region/transl_table） → Task 4。
- §5 pass_so_terms（rename＋具体 qualifier＋placeholder/unmapped 報告） → Task 5。
- §6 Change/Report → Task 2、API normalize() → Task 6、CLI → Task 7、NormalizeConfig/TOML → Task 3。
- §7 テスト（vocab/passes/normalize＋round-trip オラクル/cli/integration）・フィクスチャ → Task 1–8。

ギャップ/意図的逸脱:
- spec §6 の `Change.action` に挙げた **`dup-mapping` は per-run の Change として出さない**。重複 SO-term の解決は vocab ロード時に「具体値優先」で確定的・静的に行うため（Task 1）、利用者入力に依存しない。MVP では report ノイズを避け silent とする。
- spec §7 の「小さい golden 正規化出力」は**採用せず**、構造アサート＋round-trip（normalize→write→parse）＋validate オラクルで検証（Task 8）。golden ファイルは脆く、オラクルの方が強い回帰検出になるため。

**2. Placeholder scan**: 各ステップに実コード。"TBD"等なし。フィクスチャは実 SO-term（coding_exon/pseudogenic_CDS）。

**3. Type consistency**: `Vocab(feature_types, insdc_map, dbxref_dbtags, feature_qualifiers=…)`（Task 1、4つ目はデフォルト付きで 3-A の3引数構築と互換）、`Change(action,target,message)`（Task 2）、`NormalizeConfig(taxid,transl_table,insdc_gff_version)`（Task 3）、`NormalizeContext(vocab,seq_lengths,config)`＋`pass_directives(doc,ctx)`/`pass_so_terms(doc,ctx)->list[Change]`（Task 4,5）、`normalize(doc,*,seq_lengths,config)->(GffDocument,NormalizationReport)`＋`ALL_PASSES`（Task 6）、`cli.main(argv)->int`（Task 7）— タスク間一致。`Directive`/`Feature`/`GffDocument`/`write`/`parse`/`validate`/`Severity` は Phase1・3-A から再利用。
