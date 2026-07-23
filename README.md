# ddbj-gff

INSDC / Sequence Ontology 準拠の **GFF3 正準化ライブラリ**（パーサ・オブジェクトモデル・normalize・validate）。

各種アノテーションツール（AUGUSTUS / BRAKER / GMAP / LiftOff など）が出力する方言的な GFF3 を、
**`gene → mRNA → CDS/exon` の正準 INSDC GFF3** に整えることを目的とします。DDBJ 登録形式（MSS `.ann` /
DDBJ Record）への変換は本リポジトリの範囲外で、別リポジトリ [`ddbj_mss_tools`](https://github.com/ddbj/ddbj_mss_tools)
の `gff2mss` サブツールが担います（`gff2mss` は本ライブラリを利用します）。

- **パッケージ名**: `ddbj-gff`（import 名 `ddbj_gff`） / **version**: `0.1.0`（開発中）
- **依存**: Python **3.11 以上**、`biopython>=1.83`
- CLI（console-script）は持たず、**モジュール CLI**（`python -m ddbj_gff.normalize` など）とライブラリ API を提供します。

## できること

| 機能 | モジュール | 概要 |
|---|---|---|
| パース / 書き出し | `parser` / `writer` | GFF3 テキスト ⇄ オブジェクトモデル（`GffDocument` / `Feature` / `Span` / `Directive`） |
| 正準化 | `normalize` | 方言吸収パス群（transcript→mRNA、gene 直下 CDS/exon の mRNA 化・再ペアレント、trans-splicing、circular origin、SO 型、transl_except、anticodon、重複 locus 統合 など） |
| 検証 | `validate` | INSDC / SO プロファイルに対する QA 診断（重大度カスタマイズ可） |
| flatfile → GFF | `flatfile` | DDBJ フラットファイル（GenBank 形式）→ 正準 GFF3（flatfile↔GFF ラウンドトリップの前半） |

## インストール

2 つのリポジトリ（`gff_submission` と `ddbj_mss_tools`）を **同じ親ディレクトリに隣接**して置くと、`gff2mss`
開発時に editable path 参照が効きます。詳細な手順・別マシンへの移行は [`docs/development-setup.md`](docs/development-setup.md) を参照。

```bash
git clone git@github.com:ddbj/gff_submission.git
cd gff_submission

# uv（推奨）
uv sync                       # 依存 + dev(pytest) を .venv に構築

# もしくは pip
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e .              # biopython>=1.83
```

## 使い方

### CLI（`python -m …`）

```bash
# 正準化: 方言 GFF3 → 正準 GFF3（レポートを別途出力）
python -m ddbj_gff.normalize \
    --gff input.gff3 --fasta genome.fa \
    --config normalize.toml \
    --out normalized.gff3 --report normalize.txt

# 検証: INSDC プロファイルに対する診断（stderr に出力、重大度上書き可）
python -m ddbj_gff.validate --gff normalized.gff3
python -m ddbj_gff.validate --gff normalized.gff3 --severity SO0001=warning

# DDBJ フラットファイル → 正準 GFF3
python -m ddbj_gff.flatfile.cli --in record.gbk --out record.gff3
```

`normalize` の主なオプション: `--gff`(必須) `--fasta` `--config` `--taxid` `--transl-table`
`--insdc-gff-version` `--out` `--report`。`validate`: `--gff`(必須) `--severity CODE=LEVEL`(繰り返し可)。

### ライブラリ API

```python
from ddbj_gff import parse, write
from ddbj_gff.normalize.normalize import normalize
from ddbj_gff.normalize.config import NormalizeConfig
from ddbj_gff.validate import validate

doc = parse(open("input.gff3").read())

# 正準化（seq_lengths は circular origin 判定などに使用。無くても可）
norm, report = normalize(doc, config=NormalizeConfig(taxid=4100))
print(len(report.applied), "applied /", len(report.unresolved), "need attention")

# 検証
for d in validate(norm):
    print(d.severity.value, d.code, d.message)

# 書き出し
open("normalized.gff3", "w").write(write(norm))
```

DDBJ フラットファイルからの変換は `from ddbj_gff.flatfile.convert import flatfile_to_gff`。

### normalize 設定（`--config` の `[normalize]` テーブル）

```toml
[normalize]
taxid = 4100                       # NCBI taxid（##species / directive に使用）
transl_table = 1
insdc_gff_version = "1.0.0"
coerce_transcript_to_mrna = true   # CDS を持つ transcript を mRNA に
wrap_cds_in_mrna = true            # gene 直下 CDS/exon（mRNA 無し）を mRNA で包む
reparent_gene_children = true      # gene 直下 CDS/exon（空 mRNA あり）を mRNA へ再ペアレント
merge_overlapping_loci = false     # 同一鎖で mRNA が重なる gene を 1 locus に統合（opt-in）
merge_overlap_min_fraction = 0.0   # 統合の重なり閾値 overlap/min(len)
```

## テスト

```bash
uv run pytest        # uv の場合（既定で slow マーカーを除外）
# または venv を有効化して
pytest               # addopts = -m 'not slow'
pytest -m ""         # slow を含む全テスト
```

## リポジトリ構成

```
src/ddbj_gff/        # ライブラリ本体
  parser.py model.py writer.py io.py errors.py aa_names.py attributes.py
  normalize/         # 正準化パス群 + config + CLI
  validate/          # INSDC/SO プロファイル検証（+ data/ 語彙）
  flatfile/          # DDBJ フラットファイル → GFF3
tests/               # pytest スイート（small fixture 込み）
docs/                # 設計文書・セットアップ手順
experimental/        # 旧 ddbj/gff（2024）の実験的 GFF→MSS スクリプト（参考保管）
```

> `experimental/` は本リポジトリの前身（旧 `ddbj/gff`、2024）の実験的スクリプト・サンプルデータで、
> 参考のため保管しています。現行の正準化ライブラリ（`src/ddbj_gff/`）とは独立です。

## 関連

- [`ddbj_mss_tools`](https://github.com/ddbj/ddbj_mss_tools) — `gff2mss`（正準 GFF → MSS `.ann`）ほか。連携方法は [`docs/mss-tools-integration.md`](docs/mss-tools-integration.md)。
- 別マシンでの環境構築 / 移行手順: [`docs/development-setup.md`](docs/development-setup.md)

## 状態

開発中（version `0.1.0` 固定運用、PyPI 未公開）。API・挙動は変わり得ます。
