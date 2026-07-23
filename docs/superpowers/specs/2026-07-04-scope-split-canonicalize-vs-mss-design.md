# 責務分割設計 — gff_submission=正準化 / ddbj_mss_tools=MSS・record 変換

- 日付: 2026-07-04
- 種別: アーキテクチャ決定(境界の確定)。実移行は別 spec / writing-plans で段階実装。

## Context
GFF→DDBJ 登録ツールの最終形として、当初は「gff_submission 全体を ddbj_mss_tools へ移植」を想定していた。
本決定でこれを見直し、**境界を『正準 INSDC GFF』に置いて責務を2分割**する(結合はライブラリ依存)。

狙い:
- gff_submission を **DDBJ 非依存の汎用 GFF 正準化ライブラリ**にし、他コンシューマ(NCBI/EBI 等)へも再利用可能にする
- `.ann` / DDBJ record など **出力形式の知識を ddbj_mss_tools に一元化**(既存 `common` と同居、二重管理を解消)
- 2階層/3階層(`emit_mrna`)・transl_table・locus_tag 採番などの「出力の決定」を、分子種メタと共に出力側へ集約

egapx2mss(『正規化済み ASN.1→.tbl → .ann』)と同じ構図。`gff2mss` =『正準 GFF → .ann / record』。

## 目標アーキテクチャ
- **gff_submission = `ddbj-gff`(正準化ライブラリ)**: raw GFF(各ツール方言)→ 正準 INSDC GFF。
  構成: `src/ddbj_gff/{parser,model,writer,errors,aa_names,attributes,io}.py` + `normalize/` + `validate/`。
  将来の profile 機構(braker/maker/funannotate 等の方言吸収)もここ。
- **ddbj_mss_tools = 出力変換**: 正準 GFF + FASTA + メタ → `.ann`(MSS) / DDBJ record。
  新サブツール `src/gff2mss/`(`gff2mss.cli:main`)。`common`(COMMON/source/gap)を共有。
- **結合 = 一方向のライブラリ依存**: `gff2mss` → { `ddbj-gff`(parser/model/aa_names), `common` }。逆依存・循環は禁止。

## 「正準 INSDC GFF」の契約(境界インターフェース)
- 常に **gene→mRNA→exon/CDS の3階層**(2階層化は出力側 `emit_mrna` の責務)。
- SO→INSDC 型正規化済み、孤児 RNA 解決済み、pseudogene 整理済み、属性は INSDC 修飾子語彙。
- directives: `##gff-version 3` / `##sequence-region` / `##species`(taxid) / `#!transl_table` 等。
- 分子種(prokaryote/organelle/eukaryote)・topology・genetic code など**出力に効くメタ**は GFF 本体でなく
  `sequence_roles.tsv` / `common.json` / gff2mss config で `gff2mss` に渡す(現行 heterosigma と同方式)。

## 移すもの / 残すもの
- **移す → `ddbj_mss_tools/src/gff2mss/`**: 現 `src/ddbj_gff/mss/` 一式
  (`convert`=build_entry_features / build_cds_feature / build_rna_feature / product ルール / emit_mrna / locus_tag、
  `emit`, `config`=MssConfig, `gaps`, `translate`, `product_map`, `model`=MssDocument/MssFeature)
  + heterosigma `make_ann.py` を `gff2mss/cli.py` + `assemble.py`(common 利用)へ昇格。
- **残す → `ddbj-gff`**: parser/model/writer/errors/aa_names/attributes/io + normalize + validate。
  `validate/` は正準 GFF の INSDC/SO 準拠 QA として残置(MSS レベル検証は外部 ddbj-validator)。

## 実行アウトライン(follow-up 実装で詳細化)
1. `ddbj_mss_tools/src/gff2mss/` を新設、`mss/` をコピー、`..model` 等の相対 import を `ddbj_gff.*` へ張り替え。
2. `ddbj_mss_tools/pyproject.toml`: 依存に `ddbj-gff`(当面ローカル path/editable、将来 PyPI)、`[project.scripts]` に `gff2mss`、wheel packages に `src/gff2mss`。
3. `gff2mss/cli.py`(+`assemble.py`)= 旧 make_ann。入力: `--gff(正準) --fasta --config --common --sequence-roles --submission-category --locus-tag-start`。
4. `gff_submission` から `src/ddbj_gff/mss/` 撤去(または deprecate シム)、`validate` 残置、README/スコープを「正準化まで」に更新。
5. テスト移設: `tests/test_mss_*` → `ddbj_mss_tools/tests/`(import 張り替え)。ddbj-gff 側は parser/model/normalize/validate テストのみ。
6. heterosigma を新構成へ: 正準化=`ddbj-gff`、変換=`gff2mss`。dev/heterosigma の前処理スクリプト(split_by_compartment / build_product_map / drop_redundant_ncrna / split_fasta 等)は正準化の前処理=ddbj-gff 側 or profile として整理。

## Verification
- **一方向依存**: `src/ddbj_gff` 配下に `import common` が無い(ddbj-gff は mss_tools に非依存)。
- **ddbj-gff 単体**: parse→normalize→validate が正準 GFF を出す(既存 not-slow + slow 緑)。
- **gff2mss**: 正準 GFF + FASTA + メタ → `.ann`、`ghcr.io/ddbj/ddbj-validator:0.1.4-beta` で書式・構造エラー 0(heterosigma nuclear/organelle を移行テストベッドに使用)。
- **契約固定**: 正準 GFF のスナップショット / round-trip テストで境界を固定。

## リスク・注意
- 依存方向を厳守(mss_tools → ddbj-gff の一方向)。circular import 禁止。
- `aa_names` は normalize(transl_except)と gff2mss(翻訳)双方が使う → `ddbj-gff` 側に残し gff2mss から import。
- 実移行は本決定の承認後に別途 writing-plans で段階実装。
