# INSDC プロファイル検証器（フェーズ3-A）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** GFF3（Phase1 オブジェクトモデル）を INSDC GFF3 プロファイルに対して検証し診断を返す **detect-only** な検証器サブパッケージ `ddbj_gff.validate` を、公式 `feature-mapping.tsv`（gff3tools, Apache-2.0）同梱で実装する。

**Architecture:** 案A: ルールレジストリ。`validate(doc, *, severity_overrides=None) -> list[Diagnostic]` が独立ルール関数群を `GffDocument` に実行。各ルールは安定 code＋既定 severity で `Diagnostic`（Phase1 `errors` 再利用）を返す。統制語彙は同梱 TSV を `vocab.load_vocab()` で読む。自動修正はしない（修正は将来の 3B）。

**Tech Stack:** Python 3.11+ / 既存 `ddbj_gff`（Phase1）/ pytest / dev コンテナ。Biopython 不要。

## Global Constraints

- Python `>=3.11`、`from __future__ import annotations`。実行依存 biopython のみ（検証器は stdlib のみ使用）。
- 変更は `src/ddbj_gff/validate/` 配下。診断は Phase1 `ddbj_gff.errors`（`Severity`/`Diagnostic`/`GffParseError`）を再利用。
- **detect-only**: `validate()` は入力 `doc` を変更しない。自動修正なし。
- rule code → 既定 severity 表を持ち、`severity_overrides`（code→`error`/`warning`/`info`/`off`）で上書き。`off` は当該ルール無効化。
- 同梱データは `src/ddbj_gff/validate/data/`：`feature-mapping.tsv`（gff3tools 由来）、`dbxref.tsv`（INSDC dbxref 語彙の curated スナップショット）、`PROVENANCE.md`、`NOTICE`（Apache-2.0 帰属）。
- 特殊ケースの canonical 変換は対象外（3B）。3A は `noncanonical-special-case`(INFO) 検出のみ。
- ID必須/一意・同ID異type は Phase1 のパース時に enforce 済みで `doc.diagnostics` に出る → CLI で統合提示。検証器の独立ルールは観測可能な `multiple-parents`/`dangling-parent` を担当。
- **テストは dev コンテナ内**: 各 `pytest` を `docker exec ddbj-gff-dev uv run pytest …`。`git` は host。データ取得の `curl` は host で可（network 必要）。
- 各タスクは「失敗テスト→失敗確認→最小実装→成功確認→コミット」。

## Phase1 モデル前提（参照）
`Feature`: `.id|.type|.source|.spans|.attributes|.parent_ids|.children|.gene|.product|.note(list)|.locus_tag|.transl_table(int|None)|.is_circular(bool)|.is_trans_spliced(bool)|.dbxref(list)`。`Span`: `.seqid|.start|.end|.strand|.phase`。`GffDocument`: `.features|.roots|.feature_index|.directives|.diagnostics` ＋ properties `gff_version(str|None)`/`insdc_gff_version(str|None)`/`species(int|None)`/`sequence_regions(dict[str,(int,int)])`/`transl_table_map(dict|None)`。

---

## File Structure
| ファイル | 責務 |
|---|---|
| `src/ddbj_gff/validate/__init__.py` | `validate` 再エクスポート |
| `src/ddbj_gff/validate/vocab.py` | `Vocab` / `load_vocab()`（同梱TSV読込） |
| `src/ddbj_gff/validate/severities.py` | `DEFAULT_SEVERITIES` / `make_diagnostic` / レベル解決 |
| `src/ddbj_gff/validate/rules.py` | 個別ルール関数群＋`ALL_RULES` |
| `src/ddbj_gff/validate/validate.py` | `validate(doc, *, severity_overrides)` |
| `src/ddbj_gff/validate/cli.py`, `__main__.py` | CLI |
| `src/ddbj_gff/validate/data/` | feature-mapping.tsv / dbxref.tsv / PROVENANCE.md / NOTICE |
| `tests/test_validate_*.py`, `tests/validate_fixtures/` | テスト・フィクスチャ |

---

## Task 1: 同梱データ取得 ＋ vocab ローダ

**Files:** Create `src/ddbj_gff/validate/__init__.py`, `src/ddbj_gff/validate/vocab.py`, `src/ddbj_gff/validate/data/{feature-mapping.tsv,dbxref.tsv,PROVENANCE.md,NOTICE}`; Test `tests/test_validate_vocab.py`

**Interfaces:**
- Produces: `@dataclass(frozen=True) class Vocab(feature_types: frozenset[str], insdc_map: dict[str,str], dbxref_dbtags: frozenset[str])`; `load_vocab() -> Vocab`（同梱 data/ から読み、lru_cache）。

- [ ] **Step 1: 公式データを取得・配置（host、network 必要）**

```bash
mkdir -p src/ddbj_gff/validate/data
curl -fsSL https://raw.githubusercontent.com/enasequence/gff3tools/main/src/main/resources/feature-mapping.tsv \
  -o src/ddbj_gff/validate/data/feature-mapping.tsv
head -1 src/ddbj_gff/validate/data/feature-mapping.tsv   # 確認: 先頭が "SOID<TAB>SO term<TAB>Definition<TAB>Feature..."
wc -l src/ddbj_gff/validate/data/feature-mapping.tsv     # 確認: 100行以上
```
取得できない場合は BLOCKED 報告（controller が手当て）。

`src/ddbj_gff/validate/data/dbxref.tsv`（INSDC dbxref 語彙の curated スナップショット。1行1 DBTAG、`#` コメント可）:
```text
# Curated snapshot of common INSDC db_xref DBTAGs.
# Source: https://www.insdc.org/submitting-standards/dbxref-qualifier-vocabulary/
# This is a documented SUBSET; the dbxref-unknown-dbtag rule is WARN, so omissions are tolerable.
GenBank
RefSeq
EMBL
DDBJ
UniProtKB/Swiss-Prot
UniProtKB/TrEMBL
GeneID
GO
InterPro
PFAM
Pfam
RFAM
Rfam
EC
ECOCYC
ASAP
taxon
HGNC
MGI
FLYBASE
SGD
TAIR
Araport
Ensembl
EnsEMBL
miRBase
PDB
PubMed
dbSNP
ATCC
CDD
COG
PMID
```

`src/ddbj_gff/validate/data/PROVENANCE.md`:
```markdown
# Bundled controlled-vocabulary data — provenance

- `feature-mapping.tsv`: snapshot of `src/main/resources/feature-mapping.tsv` from
  https://github.com/enasequence/gff3tools (branch `main`). Columns:
  SOID, SO term, Definition, Feature (INSDC), Qualifier 1, Qualifier 2 (tab-separated, has header).
  License: Apache-2.0 (see NOTICE). Refresh: re-run the curl in the implementation plan Task 1.
- `dbxref.tsv`: curated subset of the INSDC db_xref vocabulary
  (https://www.insdc.org/submitting-standards/dbxref-qualifier-vocabulary/). One DBTAG per line.
```

`src/ddbj_gff/validate/data/NOTICE`:
```text
This product bundles feature-mapping.tsv from enasequence/gff3tools,
which is licensed under the Apache License, Version 2.0.
https://github.com/enasequence/gff3tools  (LICENSE: Apache-2.0)
```

`src/ddbj_gff/validate/__init__.py`:
```python
"""ddbj_gff.validate: INSDC GFF3 profile validator (phase 3A, detect-only)."""

from .validate import validate

__all__ = ["validate"]
```
（注: `validate.py` は Task 5 で作成。Task 1 ではこの `__init__` を置くと import エラーになるため、Task 1 時点では `__init__.py` を `"""..."""` のみ（`__all__ = []`）にしておき、Task 5 で再エクスポートを追加する。）

Task 1 用 `src/ddbj_gff/validate/__init__.py`（最小）:
```python
"""ddbj_gff.validate: INSDC GFF3 profile validator (phase 3A, detect-only)."""

__all__ = []
```

- [ ] **Step 2: 失敗するテストを書く**

`tests/test_validate_vocab.py`:
```python
from ddbj_gff.validate.vocab import Vocab, load_vocab


def test_load_vocab_feature_types():
    v = load_vocab()
    assert isinstance(v, Vocab)
    for term in ("CDS", "mRNA", "gene", "exon"):
        assert term in v.feature_types
    assert "totally_not_a_real_so_term" not in v.feature_types


def test_load_vocab_insdc_map_and_dbxref():
    v = load_vocab()
    assert v.insdc_map.get("CDS") == "CDS"   # SO term CDS maps to INSDC feature CDS
    assert "GenBank" in v.dbxref_dbtags
    assert "GeneID" in v.dbxref_dbtags
```

- [ ] **Step 3: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate_vocab.py -v` → FAIL（`ModuleNotFoundError: ddbj_gff.validate`）

- [ ] **Step 4: 実装**

`src/ddbj_gff/validate/vocab.py`:
```python
from __future__ import annotations

import functools
from dataclasses import dataclass
from pathlib import Path

_DATA = Path(__file__).parent / "data"


@dataclass(frozen=True)
class Vocab:
    feature_types: frozenset[str]
    insdc_map: dict[str, str]
    dbxref_dbtags: frozenset[str]


def _read_feature_mapping() -> tuple[frozenset[str], dict[str, str]]:
    terms: set[str] = set()
    mapping: dict[str, str] = {}
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
            if insdc:
                mapping[so_term] = insdc
    return frozenset(terms), mapping


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
    terms, mapping = _read_feature_mapping()
    return Vocab(feature_types=terms, insdc_map=mapping, dbxref_dbtags=_read_dbxref())
```

- [ ] **Step 5: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate_vocab.py -v` → 2 passed

- [ ] **Step 6: Commit**
```bash
git add src/ddbj_gff/validate/__init__.py src/ddbj_gff/validate/vocab.py src/ddbj_gff/validate/data tests/test_validate_vocab.py
git commit -m "feat(validate): bundle SO-INSDC feature mapping + vocab loader"
```

---

## Task 2: severities（既定重大度表 ＋ 診断ファクトリ）

**Files:** Create `src/ddbj_gff/validate/severities.py`; Test `tests/test_validate_severities.py`

**Interfaces:**
- Consumes: Phase1 `Severity`/`Diagnostic`。
- Produces: `DEFAULT_SEVERITIES: dict[str, Severity]`（全 rule code）; `make_diagnostic(code: str, message: str, line_no: int|None=None) -> Diagnostic`（既定 severity を引いて生成）; `resolve_level(name: str) -> Severity|None`（"error"/"warning"/"info"→Severity、"off"→None、不正→ValueError）。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_validate_severities.py`:
```python
import pytest
from ddbj_gff.errors import Severity
from ddbj_gff.validate.severities import DEFAULT_SEVERITIES, make_diagnostic, resolve_level


def test_default_severities_cover_known_codes():
    for code in ("missing-insdc-gff-version", "feature-type-not-insdc",
                 "noncanonical-special-case", "gene-missing-locus-tag"):
        assert code in DEFAULT_SEVERITIES


def test_make_diagnostic_uses_default_severity():
    d = make_diagnostic("missing-insdc-gff-version", "no insdc version")
    assert d.severity == Severity.ERROR
    assert d.code == "missing-insdc-gff-version"
    d2 = make_diagnostic("feature-type-not-insdc", "x")
    assert d2.severity == Severity.WARNING
    d3 = make_diagnostic("noncanonical-special-case", "x")
    assert d3.severity == Severity.INFO


def test_resolve_level():
    assert resolve_level("error") == Severity.ERROR
    assert resolve_level("warning") == Severity.WARNING
    assert resolve_level("info") == Severity.INFO
    assert resolve_level("off") is None
    with pytest.raises(ValueError):
        resolve_level("bogus")
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate_severities.py -v` → FAIL

- [ ] **Step 3: 実装**

`src/ddbj_gff/validate/severities.py`:
```python
from __future__ import annotations

from ..errors import Diagnostic, Severity

DEFAULT_SEVERITIES: dict[str, Severity] = {
    "missing-gff-version": Severity.ERROR,
    "missing-insdc-gff-version": Severity.ERROR,
    "missing-species-taxid": Severity.ERROR,
    "missing-sequence-region": Severity.ERROR,
    "duplicate-sequence-region": Severity.ERROR,
    "non-ascii": Severity.ERROR,
    "undefined-seqid": Severity.ERROR,
    "feature-outside-region": Severity.ERROR,
    "start-gt-end": Severity.ERROR,
    "feature-type-not-insdc": Severity.WARNING,
    "multiple-parents": Severity.ERROR,
    "dangling-parent": Severity.ERROR,
    "cds-missing-transl-table": Severity.ERROR,
    "cds-invalid-phase": Severity.ERROR,
    "gene-missing-locus-tag": Severity.WARNING,
    "dbxref-unknown-dbtag": Severity.WARNING,
    "noncanonical-special-case": Severity.INFO,
}

_LEVELS = {"error": Severity.ERROR, "warning": Severity.WARNING, "info": Severity.INFO}


def make_diagnostic(code: str, message: str, line_no: int | None = None) -> Diagnostic:
    return Diagnostic(DEFAULT_SEVERITIES.get(code, Severity.WARNING), line_no, code, message)


def resolve_level(name: str) -> Severity | None:
    key = name.lower()
    if key == "off":
        return None
    if key in _LEVELS:
        return _LEVELS[key]
    raise ValueError(f"invalid severity level: {name!r} (use error/warning/info/off)")
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate_severities.py -v` → 3 passed

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/validate/severities.py tests/test_validate_severities.py
git commit -m "feat(validate): default rule severities and diagnostic factory"
```

---

## Task 3: ルール群1（ヘッダ・エンコード・座標）

**Files:** Create `src/ddbj_gff/validate/rules.py`; Test `tests/test_validate_rules_header.py`

**Interfaces:**
- Consumes: `make_diagnostic`（Task 2）, `Vocab`（Task 1）, Phase1 `GffDocument`/`Feature`/`Span`。
- Produces（`rules.py`）: `rule_directives(doc, vocab)`, `rule_ascii(doc, vocab)`, `rule_seqid_bounds(doc, vocab)`, `rule_start_gt_end(doc, vocab)` — 各 `-> list[Diagnostic]`。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_validate_rules_header.py`:
```python
from ddbj_gff import parse
from ddbj_gff.model import Feature, Span, Directive, GffDocument
from ddbj_gff.validate.vocab import Vocab
from ddbj_gff.validate import rules

V = Vocab(frozenset({"gene", "mRNA", "CDS", "exon"}), {}, frozenset({"GenBank"}))


def codes(diags):
    return {d.code for d in diags}


def test_directives_all_missing():
    doc = GffDocument()  # no directives
    c = codes(rules.rule_directives(doc, V))
    assert {"missing-gff-version", "missing-insdc-gff-version",
            "missing-species-taxid", "missing-sequence-region"} <= c


def test_directives_present_ok_and_duplicate_region():
    doc = GffDocument(directives=[
        Directive("##gff-version 3", "gff-version", "3"),
        Directive("#!insdc-gff-version 1.0.0", "insdc-gff-version", "1.0.0"),
        Directive("##species ...", "species", 4530),
        Directive("##sequence-region c 1 100", "sequence-region", ("c", 1, 100)),
        Directive("##sequence-region c 1 200", "sequence-region", ("c", 1, 200)),  # dup seqid
    ])
    c = codes(rules.rule_directives(doc, V))
    assert "missing-gff-version" not in c and "missing-insdc-gff-version" not in c
    assert "missing-species-taxid" not in c and "missing-sequence-region" not in c
    assert "duplicate-sequence-region" in c


def test_ascii_flags_non_ascii_attribute():
    f = Feature("g", "S", "gene", [Span("c", 1, 9, "+")], {"product": ["protéin"]}, [])
    doc = GffDocument(features=[f])
    assert "non-ascii" in codes(rules.rule_ascii(doc, V))


def test_seqid_bounds_undefined_and_outside():
    seq_dir = Directive("x", "sequence-region", ("c", 1, 100))
    f_out = Feature("a", "S", "gene", [Span("c", 50, 150, "+")], {}, [])      # 150 > 100
    f_undef = Feature("b", "S", "gene", [Span("z", 1, 9, "+")], {}, [])       # seqid z undefined
    doc = GffDocument(directives=[seq_dir], features=[f_out, f_undef])
    c = codes(rules.rule_seqid_bounds(doc, V))
    assert "feature-outside-region" in c
    assert "undefined-seqid" in c


def test_seqid_bounds_circular_origin_spanning_allowed():
    seq_dir = Directive("x", "sequence-region", ("c", 1, 100))
    f = Feature("a", "S", "CDS", [Span("c", 90, 130, "+")], {"Is_circular": ["true"]}, [])  # end>len ok if circular
    doc = GffDocument(directives=[seq_dir], features=[f])
    assert "feature-outside-region" not in codes(rules.rule_seqid_bounds(doc, V))


def test_start_gt_end():
    f = Feature("a", "S", "gene", [Span("c", 50, 10, "+")], {}, [])
    doc = GffDocument(features=[f])
    assert "start-gt-end" in codes(rules.rule_start_gt_end(doc, V))
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate_rules_header.py -v` → FAIL（`ImportError`）

- [ ] **Step 3: 実装**

`src/ddbj_gff/validate/rules.py`:
```python
from __future__ import annotations

from .severities import make_diagnostic


def _all_strings(feature) -> list[str]:
    out = [feature.type, feature.source]
    out += [s.seqid for s in feature.spans]
    for key, values in feature.attributes.items():
        out.append(key)
        out.extend(values)
    return out


def rule_directives(doc, vocab) -> list:
    diags = []
    if doc.gff_version is None:
        diags.append(make_diagnostic("missing-gff-version", "##gff-version directive is missing"))
    if doc.insdc_gff_version is None:
        diags.append(make_diagnostic("missing-insdc-gff-version",
                                     "#!insdc-gff-version directive is missing"))
    if not isinstance(doc.species, int):
        diags.append(make_diagnostic("missing-species-taxid",
                                     "##species directive with NCBI taxid is missing"))
    seqids = [d.value[0] for d in doc.directives
              if d.kind == "sequence-region" and d.value]
    if not seqids:
        diags.append(make_diagnostic("missing-sequence-region",
                                     "##sequence-region directive is missing"))
    seen = set()
    for sid in seqids:
        if sid in seen:
            diags.append(make_diagnostic("duplicate-sequence-region",
                                         f"duplicate ##sequence-region for seqid {sid!r}"))
        seen.add(sid)
    return diags


def rule_ascii(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        if any(not s.isascii() for s in _all_strings(f)):
            diags.append(make_diagnostic("non-ascii",
                                         f"feature {f.id!r} contains non-ASCII characters"))
    return diags


def rule_seqid_bounds(doc, vocab) -> list:
    diags = []
    regions = doc.sequence_regions
    for f in doc.features:
        circular = f.is_circular
        for s in f.spans:
            if s.seqid not in regions:
                diags.append(make_diagnostic("undefined-seqid",
                                             f"feature {f.id!r} references seqid {s.seqid!r} "
                                             f"with no ##sequence-region"))
                continue
            lo, hi = regions[s.seqid]
            if (s.start < lo or s.end > hi) and not circular:
                diags.append(make_diagnostic("feature-outside-region",
                                             f"feature {f.id!r} span {s.start}..{s.end} is outside "
                                             f"sequence-region {s.seqid}:{lo}..{hi}"))
    return diags


def rule_start_gt_end(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        for s in f.spans:
            if s.start > s.end:
                diags.append(make_diagnostic("start-gt-end",
                                             f"feature {f.id!r} has start>end ({s.start}>{s.end})"))
    return diags
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate_rules_header.py -v` → 6 passed

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/validate/rules.py tests/test_validate_rules_header.py
git commit -m "feat(validate): header/encoding/coordinate rules"
```

---

## Task 4: ルール群2（SO-term・Parent・CDS・gene・Dbxref・特殊ケース）

**Files:** Modify `src/ddbj_gff/validate/rules.py`; Test `tests/test_validate_rules_body.py`

**Interfaces:**
- Produces（`rules.py` に追記）: `rule_feature_type`, `rule_parents`, `rule_cds`, `rule_gene_locus_tag`, `rule_dbxref`, `rule_special_case` — 各 `(doc, vocab) -> list[Diagnostic]`。末尾に `ALL_RULES = [...]`（全10ルール）。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_validate_rules_body.py`:
```python
from ddbj_gff.model import Feature, Span, GffDocument
from ddbj_gff.validate.vocab import Vocab
from ddbj_gff.validate import rules

V = Vocab(frozenset({"gene", "mRNA", "CDS", "exon"}), {"CDS": "CDS"}, frozenset({"GenBank", "GeneID"}))


def codes(diags):
    return {d.code for d in diags}


def test_feature_type_not_insdc():
    f = Feature("a", "S", "weird_type", [Span("c", 1, 9, "+")], {}, [])
    assert "feature-type-not-insdc" in codes(rules.rule_feature_type(GffDocument(features=[f]), V))
    f2 = Feature("b", "S", "gene", [Span("c", 1, 9, "+")], {}, [])
    assert "feature-type-not-insdc" not in codes(rules.rule_feature_type(GffDocument(features=[f2]), V))


def test_multiple_parents_and_dangling():
    parent = Feature("g", "S", "gene", [Span("c", 1, 99, "+")], {}, [])
    multi = Feature("m", "S", "mRNA", [Span("c", 1, 99, "+")], {}, ["g", "h"])   # 2 parents
    dangling = Feature("m2", "S", "mRNA", [Span("c", 1, 99, "+")], {}, ["ghost"])
    doc = GffDocument(features=[parent, multi, dangling],
                      feature_index={"g": parent, "m": multi, "m2": dangling})
    c = codes(rules.rule_parents(doc, V))
    assert "multiple-parents" in c
    assert "dangling-parent" in c


def test_cds_missing_transl_table():
    cds = Feature("c1", "S", "CDS", [Span("c", 1, 9, "+", 0)], {}, [])
    doc = GffDocument(features=[cds])  # no transl_table attr, no file-level #!transl_table
    assert "cds-missing-transl-table" in codes(rules.rule_cds(doc, V))


def test_cds_transl_table_satisfied_by_file_pragma():
    from ddbj_gff.model import Directive
    cds = Feature("c1", "S", "CDS", [Span("c", 1, 9, "+", 0)], {}, [])
    doc = GffDocument(directives=[Directive("#!transl_table primary:1", "transl_table", {"primary": 1})],
                      features=[cds])
    assert "cds-missing-transl-table" not in codes(rules.rule_cds(doc, V))


def test_cds_invalid_phase():
    cds = Feature("c1", "S", "CDS", [Span("c", 1, 9, "+", None)], {"transl_table": ["11"]}, [])  # phase None
    assert "cds-invalid-phase" in codes(rules.rule_cds(doc=GffDocument(features=[cds]), vocab=V))


def test_gene_missing_locus_tag():
    g = Feature("g", "S", "gene", [Span("c", 1, 9, "+")], {}, [])
    assert "gene-missing-locus-tag" in codes(rules.rule_gene_locus_tag(GffDocument(features=[g]), V))
    g2 = Feature("g2", "S", "gene", [Span("c", 1, 9, "+")], {"locus_tag": ["X_1"]}, [])
    assert "gene-missing-locus-tag" not in codes(rules.rule_gene_locus_tag(GffDocument(features=[g2]), V))


def test_dbxref_unknown_dbtag():
    f = Feature("a", "S", "gene", [Span("c", 1, 9, "+")], {"Dbxref": ["GenBank:X1", "WeirdDB:9"]}, [])
    c = codes(rules.rule_dbxref(GffDocument(features=[f]), V))
    assert "dbxref-unknown-dbtag" in c   # WeirdDB unknown


def test_special_case_detection():
    ts = Feature("g", "S", "gene", [Span("c", 1, 9, "-")], {"exception": ["trans-splicing"]}, [])
    te = Feature("c1", "S", "CDS", [Span("c", 1, 9, "+", 0)],
                 {"transl_except": ["(pos:1..3,aa:Sec)"], "transl_table": ["11"]}, [])
    c = codes(rules.rule_special_case(GffDocument(features=[ts, te]), V))
    assert "noncanonical-special-case" in c


def test_all_rules_list_complete():
    assert len(rules.ALL_RULES) == 10
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate_rules_body.py -v` → FAIL

- [ ] **Step 3: 実装** — `rules.py` 末尾に追記:
```python
def rule_feature_type(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        if f.type not in vocab.feature_types:
            diags.append(make_diagnostic("feature-type-not-insdc",
                                         f"feature {f.id!r} type {f.type!r} is not an INSDC-supported "
                                         f"SO term"))
    return diags


def rule_parents(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        if len(f.parent_ids) > 1:
            diags.append(make_diagnostic("multiple-parents",
                                         f"feature {f.id!r} has {len(f.parent_ids)} parents "
                                         f"(INSDC allows a single parent per row)"))
        for pid in f.parent_ids:
            if pid not in doc.feature_index:
                diags.append(make_diagnostic("dangling-parent",
                                             f"feature {f.id!r} references missing Parent {pid!r}"))
    return diags


def rule_cds(doc, vocab) -> list:
    diags = []
    has_file_table = doc.transl_table_map is not None
    for f in doc.features:
        if f.type != "CDS":
            continue
        if f.transl_table is None and not has_file_table:
            diags.append(make_diagnostic("cds-missing-transl-table",
                                         f"CDS {f.id!r} lacks transl_table and no file-level "
                                         f"#!transl_table is present"))
        for s in f.spans:
            if s.phase not in (0, 1, 2):
                diags.append(make_diagnostic("cds-invalid-phase",
                                             f"CDS {f.id!r} has invalid phase {s.phase!r}"))
    return diags


def rule_gene_locus_tag(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        if f.type == "gene" and not f.locus_tag:
            diags.append(make_diagnostic("gene-missing-locus-tag",
                                         f"gene {f.id!r} has no locus_tag"))
    return diags


def rule_dbxref(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        for xref in f.dbxref:
            dbtag = xref.split(":", 1)[0]
            if dbtag and dbtag not in vocab.dbxref_dbtags:
                diags.append(make_diagnostic("dbxref-unknown-dbtag",
                                             f"feature {f.id!r} Dbxref DBTAG {dbtag!r} is not in the "
                                             f"INSDC vocabulary"))
    return diags


def rule_special_case(doc, vocab) -> list:
    diags = []
    for f in doc.features:
        if f.is_trans_spliced and "location" not in f.attributes:
            diags.append(make_diagnostic("noncanonical-special-case",
                                         f"feature {f.id!r} uses non-canonical trans-splicing "
                                         f"representation (no location= attribute)"))
        if "transl_except" in f.attributes:
            diags.append(make_diagnostic("noncanonical-special-case",
                                         f"feature {f.id!r} uses transl_except attribute "
                                         f"(canonical form is a recoded_codon child feature)"))
        if "anticodon" in f.attributes:
            diags.append(make_diagnostic("noncanonical-special-case",
                                         f"feature {f.id!r} uses anticodon attribute "
                                         f"(canonical form is an anticodon child feature)"))
    return diags


ALL_RULES = [
    rule_directives,
    rule_ascii,
    rule_seqid_bounds,
    rule_start_gt_end,
    rule_feature_type,
    rule_parents,
    rule_cds,
    rule_gene_locus_tag,
    rule_dbxref,
    rule_special_case,
]
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate_rules_body.py -v` → 9 passed

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/validate/rules.py tests/test_validate_rules_body.py
git commit -m "feat(validate): SO-term/parent/CDS/gene/dbxref/special-case rules + ALL_RULES"
```

---

## Task 5: validate() 集約 ＋ severity 上書き

**Files:** Create `src/ddbj_gff/validate/validate.py`; Modify `src/ddbj_gff/validate/__init__.py`; Test `tests/test_validate.py`

**Interfaces:**
- Consumes: `ALL_RULES`（Task 4）, `load_vocab`（Task 1）, `resolve_level`（Task 2）。
- Produces: `validate(doc, *, severity_overrides: dict[str,str]|None=None) -> list[Diagnostic]`（全ルール実行、override 適用[off で除外/他は severity 差替]、(line_no, code) 順整列、入力 doc 不変）。`__init__` から `validate` 再エクスポート。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_validate.py`:
```python
from ddbj_gff import parse
from ddbj_gff.errors import Severity
from ddbj_gff.validate import validate

GFF_BAD = (
    "##gff-version 3\n"
    "chr1\tS\tgene\t1\t9\t.\t+\t.\tID=g\n"   # no #!insdc-gff-version, no ##species, no ##sequence-region, gene no locus_tag
)


def test_validate_runs_all_and_finds_errors():
    diags = validate(parse(GFF_BAD))
    codes = {d.code for d in diags}
    assert "missing-insdc-gff-version" in codes
    assert "missing-sequence-region" in codes
    assert "gene-missing-locus-tag" in codes


def test_severity_override_off_and_promote():
    base = {d.code for d in validate(parse(GFF_BAD))}
    assert "gene-missing-locus-tag" in base
    # off removes it
    off = validate(parse(GFF_BAD), severity_overrides={"gene-missing-locus-tag": "off"})
    assert "gene-missing-locus-tag" not in {d.code for d in off}
    # promote to error
    promoted = validate(parse(GFF_BAD), severity_overrides={"gene-missing-locus-tag": "error"})
    g = [d for d in promoted if d.code == "gene-missing-locus-tag"][0]
    assert g.severity == Severity.ERROR


def test_validate_sorted_and_does_not_mutate_doc():
    doc = parse(GFF_BAD)
    before = len(doc.diagnostics)
    diags = validate(doc)
    keys = [(d.line_no if d.line_no is not None else -1, d.code) for d in diags]
    assert keys == sorted(keys)
    assert len(doc.diagnostics) == before  # detect-only: validate must not append to doc.diagnostics
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate.py -v` → FAIL（`ImportError: cannot import name 'validate'`）

- [ ] **Step 3: 実装**

`src/ddbj_gff/validate/validate.py`:
```python
from __future__ import annotations

import dataclasses

from .rules import ALL_RULES
from .severities import resolve_level
from .vocab import load_vocab


def validate(doc, *, severity_overrides: dict[str, str] | None = None) -> list:
    overrides = severity_overrides or {}
    vocab = load_vocab()
    diags: list = []
    for rule in ALL_RULES:
        diags.extend(rule(doc, vocab))

    out: list = []
    for d in diags:
        if d.code in overrides:
            level = resolve_level(overrides[d.code])
            if level is None:        # "off"
                continue
            d = dataclasses.replace(d, severity=level)
        out.append(d)

    out.sort(key=lambda d: (d.line_no if d.line_no is not None else -1, d.code))
    return out
```

`src/ddbj_gff/validate/__init__.py` を更新:
```python
"""ddbj_gff.validate: INSDC GFF3 profile validator (phase 3A, detect-only)."""

from .validate import validate

__all__ = ["validate"]
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate.py -v` → 3 passed。全体 `docker exec ddbj-gff-dev uv run pytest -q`。

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/validate/validate.py src/ddbj_gff/validate/__init__.py tests/test_validate.py
git commit -m "feat(validate): validate() aggregation with severity overrides"
```

---

## Task 6: CLI

**Files:** Create `src/ddbj_gff/validate/cli.py`, `src/ddbj_gff/validate/__main__.py`; Test `tests/test_validate_cli.py`

**Interfaces:**
- Consumes: `parse`（Phase1）, `validate`（Task 5）, `Severity`。
- Produces: `cli.main(argv=None) -> int`（`--gff` 必須、`--severity CODE=LEVEL`（複数可）。parse の `doc.diagnostics` ＋ `validate(doc, overrides)` を統合し重大度別サマリを stderr に。ERROR があれば 1）。`python -m ddbj_gff.validate` 実行可。

- [ ] **Step 1: 失敗するテストを書く**

`tests/test_validate_cli.py`:
```python
from ddbj_gff.validate.cli import main

VALID = (
    "##gff-version 3\n"
    "#!insdc-gff-version 1.0.0\n"
    "##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=3702\n"
    "##sequence-region chr1 1 1000\n"
    "chr1\tS\tgene\t1\t99\t.\t+\t.\tID=g;locus_tag=ABC_1\n"
    "chr1\tS\tmRNA\t1\t99\t.\t+\t.\tID=m;Parent=g\n"
    "chr1\tS\texon\t1\t99\t.\t+\t.\tID=e;Parent=m\n"
    "chr1\tS\tCDS\t1\t9\t.\t+\t0\tID=c;Parent=m;transl_table=1\n"
)
BAD = "##gff-version 3\nchr1\tS\tgene\t1\t9\t.\t+\t.\tID=g\n"


def test_cli_valid_returns_zero(tmp_path):
    p = tmp_path / "v.gff"; p.write_text(VALID)
    assert main(["--gff", str(p)]) == 0


def test_cli_bad_returns_one_and_reports(tmp_path, capsys):
    p = tmp_path / "b.gff"; p.write_text(BAD)
    rc = main(["--gff", str(p)])
    assert rc == 1
    err = capsys.readouterr().err
    assert "ERROR" in err
    assert "missing-insdc-gff-version" in err


def test_cli_severity_override(tmp_path):
    # turning the only ERRORs off should make rc 0 — promote nothing, just off the blockers
    p = tmp_path / "b.gff"; p.write_text(BAD)
    rc = main(["--gff", str(p),
               "--severity", "missing-insdc-gff-version=off",
               "--severity", "missing-species-taxid=off",
               "--severity", "missing-sequence-region=off",
               "--severity", "undefined-seqid=off"])
    assert rc == 0
```

- [ ] **Step 2: 失敗確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate_cli.py -v` → FAIL

- [ ] **Step 3: 実装**

`src/ddbj_gff/validate/cli.py`:
```python
from __future__ import annotations

import argparse
import sys

from .. import parse
from ..errors import Severity
from .validate import validate


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="ddbj_gff.validate",
                                 description="Validate a GFF3 file against the INSDC profile")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--severity", action="append", default=[],
                    metavar="CODE=LEVEL", help="override a rule severity (error/warning/info/off)")
    args = ap.parse_args(argv)

    overrides: dict[str, str] = {}
    for item in args.severity:
        if "=" not in item:
            ap.error(f"--severity expects CODE=LEVEL, got {item!r}")
        code, level = item.split("=", 1)
        overrides[code.strip()] = level.strip()

    with open(args.gff, encoding="ascii", errors="replace") as fh:
        doc = parse(fh.read())

    diags = list(doc.diagnostics) + validate(doc, severity_overrides=overrides)

    counts: dict[str, int] = {}
    for d in diags:
        counts[d.severity.value] = counts.get(d.severity.value, 0) + 1
        sys.stderr.write(f"{d.severity.value}\t{d.code}\t{d.message}\n")
    sys.stderr.write("summary: "
                     + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())) + "\n")
    return 1 if counts.get(Severity.ERROR.value) else 0
```

`src/ddbj_gff/validate/__main__.py`:
```python
import sys

from .cli import main

sys.exit(main())
```

- [ ] **Step 4: 成功確認** — `docker exec ddbj-gff-dev uv run pytest tests/test_validate_cli.py -v` → 3 passed。全体 `docker exec ddbj-gff-dev uv run pytest -q`。

- [ ] **Step 5: Commit**
```bash
git add src/ddbj_gff/validate/cli.py src/ddbj_gff/validate/__main__.py tests/test_validate_cli.py
git commit -m "feat(validate): CLI entry point"
```

---

## Task 7: フィクスチャ ＋ 実ファイル統合テスト（slow）

**Files:** Create `tests/validate_fixtures/valid_insdc.gff3`; Test `tests/test_validate_integration.py`

**Interfaces:** Consumes `parse`/`validate`、実 example（存在時のみ）。

- [ ] **Step 1: 有効フィクスチャ＋テストを書く**

`tests/validate_fixtures/valid_insdc.gff3`（TAB区切り、全必須ディレクティブ＋許可SO term＋locus_tag/transl_table）:
```text
##gff-version 3
#!insdc-gff-version 1.0.0
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=3702
##sequence-region chr1 1 10000
chr1	S	gene	100	900	.	+	.	ID=g1;locus_tag=ABC_000010
chr1	S	mRNA	100	900	.	+	.	ID=m1;Parent=g1
chr1	S	exon	100	900	.	+	.	ID=e1;Parent=m1
chr1	S	CDS	130	870	.	+	0	ID=c1;Parent=m1;transl_table=1;Dbxref=GeneID:123
```

`tests/test_validate_integration.py`:
```python
import gzip
from pathlib import Path
import pytest
from ddbj_gff import parse
from ddbj_gff.errors import Severity
from ddbj_gff.validate import validate

FIX = Path(__file__).parent / "validate_fixtures"
ROOT = Path(__file__).resolve().parents[1]


def test_valid_insdc_fixture_has_no_errors():
    diags = validate(parse((FIX / "valid_insdc.gff3").read_text()))
    errors = [d for d in diags if d.severity == Severity.ERROR]
    assert errors == [], f"unexpected errors: {[(d.code, d.message) for d in errors]}"


@pytest.mark.slow
def test_rice_cp_flags_insdc_violations():
    p = ROOT / "examples" / "rice_cp" / "rice_cp.gff3"
    if not p.exists():
        pytest.skip(f"missing {p}")
    codes = {d.code for d in validate(parse(p.read_text(errors="replace")))}
    # NCBI-style file lacks the INSDC version directive
    assert "missing-insdc-gff-version" in codes
    # rps12 uses exception=trans-splicing + part= (non-canonical)
    assert "noncanonical-special-case" in codes


@pytest.mark.slow
def test_ecoli_flags_transl_except_noncanonical():
    p = ROOT / "examples" / "ecoli" / "GCF_000005845.2_ASM584v2_genomic.gff.gz"
    if not p.exists():
        pytest.skip(f"missing {p}")
    text = gzip.decompress(p.read_bytes()).decode("ascii", errors="replace")
    codes = {d.code for d in validate(parse(text))}
    assert "missing-insdc-gff-version" in codes
    assert "noncanonical-special-case" in codes   # transl_except attribute present
```

- [ ] **Step 2: 実行**
`docker exec ddbj-gff-dev uv run pytest tests/test_validate_integration.py -v` （非slowの valid フィクスチャは必ず実行・pass）。
`docker exec ddbj-gff-dev uv run pytest tests/test_validate_integration.py -m slow -v` （rice_cp/ecoli 実ファイルで期待 code を検証）。
もし valid フィクスチャに想定外 ERROR が出たら、その code を調べ、フィクスチャを INSDC 準拠に直す（規則実装は変えない）。期待: `missing-insdc-gff-version` は rice_cp/ecoli で出る、`noncanonical-special-case` も出る。

- [ ] **Step 3: 全体確認** — `docker exec ddbj-gff-dev uv run pytest -q`（slow除外）→ 全 pass。

- [ ] **Step 4: Commit**
```bash
git add tests/validate_fixtures/valid_insdc.gff3 tests/test_validate_integration.py
git commit -m "test(validate): valid fixture + slow integration on real examples"
```

---

## Self-Review

**1. Spec coverage**（spec §→タスク）:
- §3 構成・同梱データ・vocab → Task 1
- §4 規則セット・重大度 → Task 2(severities)/Task 3(header/encoding/seqid)/Task 4(so-term/parent/cds/gene/dbxref/special)
- §5 validate API・severity 上書き・CLI・統合提示・exit → Task 5/Task 6
- §6 テスト（vocab/rules/validate/cli/integration）・valid フィクスチャ → Task 1-7

ギャップ: spec §4 の `missing-id-with-children`/`duplicate-id-different-type` は Phase1 のパース時 enforce・診断であり CLI で `doc.diagnostics` として統合提示する設計（Global Constraints に明記）。検証器の独立ルールは観測可能な `multiple-parents`/`dangling-parent` で代替。spec の検査意図は満たす。

**2. Placeholder scan**: 各ステップに実コード。"TBD"等なし。Task 1 のデータは curl 取得＋検証手順、Task 7 の valid フィクスチャは「想定外 ERROR が出たら準拠に直す」手順を明示（プレースホルダでない）。

**3. Type consistency**: `Vocab(feature_types, insdc_map, dbxref_dbtags)`、`load_vocab()`、`make_diagnostic(code,message,line_no)`、`resolve_level(name)->Severity|None`、各 `rule_*(doc, vocab)->list[Diagnostic]`、`ALL_RULES`（10件）、`validate(doc,*,severity_overrides)`、`cli.main(argv)->int` — タスク間一致。`Diagnostic`/`Severity` は Phase1 再利用、override は `dataclasses.replace` で severity 差替（Diagnostic は frozen）。
