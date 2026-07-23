# gff2mss Migration Implementation Plan (責務分割: ddbj-gff 正準化 / ddbj_mss_tools 変換)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `src/ddbj_gff/mss/`(MSS feature 生成)と heterosigma `make_ann.py` を `ddbj_mss_tools` の新サブツール `gff2mss` へ移し、`gff_submission` を GFF 正準化ライブラリ `ddbj-gff` に絞る。境界=正準 INSDC GFF、依存は一方向 `gff2mss → {ddbj-gff, common}`。

**Architecture:** `mss/` を丸ごと `ddbj_mss_tools/src/gff2mss/` へコピーし、親パッケージ参照(`..errors`/`.. import aa_names`/`.. parse`/`..io`)だけを `ddbj_gff.*` に張り替える。`make_ann.build_ann_text` を `gff2mss/assemble.py` + `gff2mss/cli.py` に昇格。既存 mss テストを移設して緑を維持。gff_submission からは `mss/` を撤去、`validate` は残す。

**Tech Stack:** Python 3.11+/3.12, BioPython, pydantic>=2(common), hatchling。テストは amd64 コンテナ `ddbj-gff-dev`。DDBJ 検証は `ghcr.io/ddbj/ddbj-validator:0.1.4-beta`。

## Global Constraints
- 依存は **一方向**: `gff2mss → {ddbj-gff, common}`。`src/ddbj_gff/` 配下に `import common` / mss_tools 参照を作らない(循環禁止)。
- 移動する `mss/` の**挙動は不変**(既存 mss テストがそのまま緑であること)。
- 親パッケージ参照の張り替えは次の6行のみ(他は intra-package 相対のまま):
  `config.py:6 from ..errors→ddbj_gff.errors` / `convert.py:13 from ..errors→ddbj_gff.errors` / `convert.py:14 from .. import aa_names→from ddbj_gff import aa_names` / `cli.py:8 from .. import parse→from ddbj_gff import parse` / `cli.py:9 from ..io import open_text→ddbj_gff.io` / `cli.py:14 from ..errors import Severity→ddbj_gff.errors`。
- 正準 GFF の契約: 常に gene→mRNA→exon/CDS 3階層、SO→INSDC 型正規化済み、directives 付与済み。2階層化・transl_table・locus_tag 等は出力側(gff2mss の MssConfig / sequence_roles / common.json)。
- テスト環境(コンテナ `ddbj-gff-dev`): `ddbj_gff` は `/workspace/src`(gff_submission bind mount)、`gff2mss`+`common` は `/opt/mss_src`(`docker cp` で同期)、`pydantic>=2` 導入済み。
  テストコマンド雛形:
  `docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src ddbj-gff-dev bash -lc 'cd /workspace && uv run pytest <args>'`
- heterosigma を移行のテストベッドにする(nuclear/organelle が ddbj-validator で書式・構造エラー 0 を維持)。

## File Structure
- `ddbj_mss_tools/src/gff2mss/` (新規パッケージ):
  - `convert.py, emit.py, config.py, gaps.py, locus_tag.py, translate.py, product_map.py, model.py` … 現 `mss/` からコピー(親参照のみ張替)
  - `assemble.py` (新) … 旧 `make_ann.build_ann_text`(common で COMMON/source/gap)
  - `cli.py` (新) … `gff2mss` の引数解析 → assemble
  - `__init__.py, __main__.py`
  - 旧 `mss/cli.py`(自前 COMMON の standalone .ann) → `gff2mss/mss_ann.py` にリネーム保持(既存 test_mss_cli/snapshot 用)
- `ddbj_mss_tools/pyproject.toml` (修正) … `ddbj-gff` 依存, `[project.scripts] gff2mss`, wheel packages
- `ddbj_mss_tools/tests/` … 現 `gff_submission/tests/test_mss_*.py` を移設(import 張替)
- `gff_submission/src/ddbj_gff/mss/` (削除), `gff_submission/pyproject.toml`+README (スコープ更新)
- `gff_submission/dev/heterosigma/scripts/make_ann.py` (撤去; `gff2mss` CLI を使用)

---

## Task 1: ddbj-gff の一方向性を確認し wheel 化可能にする
**Files:** Modify: `gff_submission/pyproject.toml`(description/keywords をライブラリ向けに)。Test: 既存スイート。

**Interfaces:** Produces: `ddbj-gff` wheel(`ddbj_gff` パッケージ: parse/GffDocument/Feature/Span、normalize、validate、io、aa_names、errors)。

- [ ] **Step 1: 一方向依存の確認(mss を除いた ddbj_gff が common/mss_tools に非依存)**
```bash
grep -rnE "import common|ddbj_mss_tools|from \.mss|import mss" src/ddbj_gff | grep -v "/mss/"
```
Expected: 出力なし(mss/ 以外は common/mss_tools に非依存)。
- [ ] **Step 2: wheel ビルド確認(コンテナ)**
```bash
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && uv build --wheel 2>&1 | tail -3 && ls dist/*.whl'
```
Expected: `ddbj_gff-*.whl` が生成。
- [ ] **Step 3: Commit(あれば)**
```bash
git add pyproject.toml && git commit -m "chore(ddbj-gff): position as standalone canonicalization library" || echo "no change"
```

## Task 2: gff2mss パッケージ作成(mss/ コピー + 親参照張替)
**Files:** Create: `ddbj_mss_tools/src/gff2mss/{convert,emit,config,gaps,locus_tag,translate,product_map,model,__init__}.py`, `gff2mss/mss_ann.py`(旧 cli), `gff2mss/__main__.py`。

**Interfaces:** Produces: `gff2mss.convert.build_entry_features(doc, seqs, cfg, diagnostics)`, `gff2mss.emit.feature_rows/emit_fasta`, `gff2mss.config.load_config/MssConfig`, `gff2mss.product_map.load_product_map`, `gff2mss.model.MssFeature/MssQualifier/MssDocument/MssEntry`。すべて `ddbj_gff.{errors,aa_names,parse,io}` に依存。

- [ ] **Step 1: mss/ をコピー**
```bash
mkdir -p ../ddbj_mss_tools/src/gff2mss
cp src/ddbj_gff/mss/{convert,emit,config,gaps,locus_tag,translate,product_map,model,__init__}.py ../ddbj_mss_tools/src/gff2mss/
cp src/ddbj_gff/mss/cli.py ../ddbj_mss_tools/src/gff2mss/mss_ann.py
```
- [ ] **Step 2: 親参照 6 行を `ddbj_gff.*` へ張替**（config.py/convert.py/mss_ann.py の該当行。Global Constraints 参照）。`mss_ann.py` は `from .convert import convert` / `from .emit import ...` の intra 参照は据え置き。
- [ ] **Step 3: `__main__.py` を gff2mss.cli に向ける**
```python
from .cli import main
if __name__ == "__main__":
    raise SystemExit(main())
```
- [ ] **Step 4: import 健全性チェック(コンテナ; gff2mss を /opt/mss_src へ同期)**
```bash
docker cp ../ddbj_mss_tools/src/gff2mss ddbj-gff-dev:/opt/mss_src/gff2mss
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src ddbj-gff-dev bash -lc 'cd /workspace && uv run python -c "import sys; sys.path.insert(0,\"/workspace/src\"); from gff2mss.convert import build_entry_features; from gff2mss.emit import feature_rows; from gff2mss.config import load_config; print(\"gff2mss import OK\")"'
```
Expected: `gff2mss import OK`。

## Task 3: gff2mss/assemble.py + cli.py（make_ann の昇格）
**Files:** Create: `ddbj_mss_tools/src/gff2mss/assemble.py`, `ddbj_mss_tools/src/gff2mss/cli.py`。

**Interfaces:** Produces: `gff2mss.assemble.build_ann_text(gff, fasta, mss_config, common, sequence_roles, submission_category, locus_tag_start=None) -> (ann_text:str, out_seqs:dict)`; `gff2mss.cli.main(argv=None)`。Consumes: `gff2mss.{convert,emit,config,product_map}`, `ddbj_gff.{parse,io}`, `common.{models,common_builder,source_builder,gap_annotator,submission_category}`。

- [ ] **Step 1: `assemble.py` を作成** = 現 `dev/heterosigma/scripts/make_ann.py` の `build_ann_text` をそのまま移植し、import を `ddbj_gff.mss.*` → `gff2mss.*` に変更(`from gff2mss.config import load_config` 等)。`sys.path` 挿入は削除(パッケージ依存で解決)。全 FASTA 配列を出力する現行実装(source-only エントリ含む)を維持。
- [ ] **Step 2: `cli.py` を作成**（argparse: `--gff --fasta --mss-config --common --sequence-roles --submission-category --locus-tag-start --out` → `assemble.build_ann_text` → `.ann`/`.fasta` 書き出し。現 make_ann の `main()` 相当)。
- [ ] **Step 3: 動作確認(小 fixture、コンテナ)** — organelle circular + 核 WGS + source-only の3点を assemble で検証(現 `dev/heterosigma/scripts/tests/test_make_ann.py` を `gff2mss` import に張替えて実行):
```bash
docker cp ../ddbj_mss_tools/src/gff2mss ddbj-gff-dev:/opt/mss_src/gff2mss
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev bash -lc 'cd /workspace && uv run pytest ../ddbj_mss_tools/tests/test_gff2mss_assemble.py -v -p no:cacheprovider' 2>&1 | tail -8
```
Expected: PASS(TOPOLOGY circular / submitter_seqid / product / source-only を確認)。

## Task 4: pyproject 配線（依存・スクリプト・パッケージ）
**Files:** Modify: `ddbj_mss_tools/pyproject.toml`（`[project.dependencies]` に `ddbj-gff`、`[project.scripts]` に `gff2mss = "gff2mss.cli:main"`、`[tool.hatch.build.targets.wheel] packages` に `src/gff2mss`）。`ddbj_mss_tools/Dockerfile`/requirements に ddbj-gff の入手手段(当面はローカル wheel を pip install、将来 PyPI)。

- [ ] **Step 1: pyproject に追記**（dependencies に `"ddbj-gff"`、scripts に gff2mss、wheel packages に `src/gff2mss`）。
- [ ] **Step 2: ローカル依存の解決方法を明記**（`[tool.uv.sources] ddbj-gff = { path = "../gff_submission", editable = true }` を追加、または CI/コンテナでは Task1 の wheel を `pip install`)。
- [ ] **Step 3: 検証** — mss-tools 環境(または ddbj-gff-dev + PYTHONPATH)で `python -c "import gff2mss.cli"` と `gff2mss --help` 相当が通ること。

## Task 5: テスト移設
**Files:** Move: `gff_submission/tests/test_mss_*.py` → `ddbj_mss_tools/tests/`（import `ddbj_gff.mss.*` → `gff2mss.*`；`test_mss_cli.py`/`test_mss_snapshot.py` は `gff2mss.mss_ann` を参照）。Create: `ddbj_mss_tools/tests/test_gff2mss_assemble.py`（make_ann テスト移植）。

- [ ] **Step 1: 21 個の test_mss_*.py をコピーし import を張替**（`from ddbj_gff.mss.X` → `from gff2mss.X`）。
- [ ] **Step 2: 全 mss テストを新配置で実行(コンテナ)**
```bash
docker cp ../ddbj_mss_tools/src/gff2mss ddbj-gff-dev:/opt/mss_src/gff2mss
docker cp ../ddbj_mss_tools/tests ddbj-gff-dev:/opt/mss_tests
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv -e PYTHONPATH=/opt/mss_src:/workspace/src ddbj-gff-dev bash -lc 'uv run pytest /opt/mss_tests -q'
```
Expected: 旧 mss テスト数（約 100+）が緑（slow 含む marchantia は fixture パスに注意、必要なら skip）。

## Task 6: gff_submission から mss/ 撤去 + スコープ更新
**Files:** Delete: `gff_submission/src/ddbj_gff/mss/`, `gff_submission/tests/test_mss_*.py`。Modify: `gff_submission/pyproject.toml`(description=canonicalizer)、`README`/`docs/project_goal` にスコープ「正準化まで」。

- [ ] **Step 1: mss/ と test_mss_* を削除**
```bash
git rm -r src/ddbj_gff/mss tests/test_mss_*.py
```
- [ ] **Step 2: ddbj-gff 残存スイートが緑（parser/model/normalize/validate/io）**
```bash
docker exec -e UV_PROJECT_ENVIRONMENT=/opt/ddbj-venv ddbj-gff-dev bash -lc 'cd /workspace && uv run pytest -q'
```
Expected: mss 以外の全テスト緑。
- [ ] **Step 3: Commit**
```bash
git add -A && git commit -m "refactor(ddbj-gff): drop mss/ (moved to ddbj_mss_tools/gff2mss); scope = GFF canonicalization"
```

## Task 7: heterosigma を新構成へ + end-to-end 検証
**Files:** Modify: `gff_submission/dev/heterosigma/scripts/run_pipeline.sh`(make_ann → `gff2mss` 実行)、`make_ann.py` 撤去。

- [ ] **Step 1: run_pipeline を `gff2mss` 実行に変更**（`python scripts/make_ann.py …` → `gff2mss --gff … --fasta … --mss-config … --common … [--sequence-roles …] --submission-category … [--locus-tag-start …] --out …`）。正準化(normalize)は `ddbj-gff` の `python -m ddbj_gff.normalize` のまま。
- [ ] **Step 2: nuclear/organelle を再生成し ddbj-validator で確認**
```bash
# gff2mss + ddbj-gff が import 可能な環境で run_pipeline を実行後:
cd dev/heterosigma && docker run --rm -v "$PWD":/data -w /data ghcr.io/ddbj/ddbj-validator:0.1.4-beta -o submission/nuc_out submission/nuc
```
Expected: 書式・構造エラー 0（残るのは BioProject/BioSample の実DB照合のみ）。nuclear.fasta=66/organelle.fasta=2、organelle mRNA=0、locus_tag 衝突 0 を維持。
- [ ] **Step 3: Commit(dev は方針により未コミットでも可)**

---

## Self-Review
- **Spec coverage**: 目標アーキテクチャ(gff2mss 新設=Task2-4 / ddbj-gff 残置=Task6 / 一方向依存=Task1,4)/契約(3階層・出力側 emit_mrna=移動コードで保持)/移す・残す(Task2,5,6)/verification(一方向 grep=Task1、ddbj-validator=Task7)。網羅。
- **Placeholder scan**: TBD なし。各 Step にコマンド/対象を明記。
- **Type consistency**: `build_ann_text(...) -> (str, dict)`、`build_entry_features(doc, seqs, cfg, diagnostics)`、`feature_rows(feat)`、`MssConfig`/`load_config` は移動前後で同一(挙動不変制約)。親参照張替は6行に限定。
- **注意**: mss-tools の本番 CI/コンテナ(app スタック)への統合は本計画では ddbj-gff-dev コンテナ検証までとし、mss-tools 環境での CI 配線は follow-up。ローカル依存解決(Task4 Step2)は path/editable を既定、PyPI 公開は将来。
