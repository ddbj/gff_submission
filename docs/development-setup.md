# 開発環境の構築（別マシンへの移行手順）

`gff_submission`（Python パッケージ名 **`ddbj-gff`** / import 名 `ddbj_gff`）の開発環境を
新しいマシンで再構築するための手順。本リポジトリは GFF3 正規化・パース等を提供する**ライブラリ**です。
console-script（`pip` でインストールされる実行コマンド。`[project.scripts]`）は持ちませんが、
`python -m ddbj_gff.normalize` / `python -m ddbj_gff.validate` のモジュール CLI があります
（`src/ddbj_gff/flatfile/cli.py` も同様）。主に `ddbj_mss_tools` の `gff2mss` サブツールから
ライブラリとして利用されます。

## 前提

| 項目 | 要件 |
|---|---|
| Python | **3.11 以上**（`pyproject.toml` の `requires-python = ">=3.11"`） |
| uv | **推奨**。無い場合は pip でも可 |
| git | 必須 |

## 新マシンへの取得

このリポジトリは GitHub にあります（**`git@github.com:ddbj/gff_submission.git`**。旧名 `ddbj/gff` は自動リダイレクト）:

```bash
git clone git@github.com:ddbj/gff_submission.git
```

> **注記（`uv.lock`）**: `uv.lock` は `.gitignore` 対象で**追跡されていません**。clone には含まれないため、
> 新マシンの `uv sync` は依存を**その場で再解決**します（固定 lock には従いません）。依存は `biopython>=1.83`
> の1つのみのため実害は小さいです。再現性を厳密に固定したい場合は `.gitignore` から `uv.lock` を外して追跡してください。

`ddbj_mss_tools` の `gff2mss` を開発する場合、両リポジトリを**同じ親ディレクトリに隣接**して置きます
（`ddbj_mss_tools` 側が `../gff_submission` を editable 参照するため）:

```
<親ディレクトリ>/
├── gff_submission/     # このリポジトリ
└── ddbj_mss_tools/
```

## セットアップ

### uv を使う場合（推奨）

```bash
cd gff_submission
uv sync            # 依存 + dev グループ(pytest) を .venv に構築（uv.lock があればそれに従う。上記注記参照）
```

### pip を使う場合

```bash
cd gff_submission
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e .   # 依存: biopython>=1.83
pip install pytest # dev グループ相当
```

## テストの実行

```bash
uv run pytest      # uv の場合
# または venv を有効化して
pytest             # 既定で slow マーカーを除外（pyproject の addopts = -m 'not slow'）
pytest -m ""       # slow を含む全テストを実行したい場合
```

## ddbj_mss_tools との連携

- `ddbj-gff` は `ddbj_mss_tools` の `gff2mss` サブツール**専用の optional 依存**です。他ツールは依存しません。
- 開発時の反映（uv editable path / コンテナ用 wheel 生成）とバージョン運用は `docs/mss-tools-integration.md` を参照。
- **バージョンは現状 `0.1.0` 固定**。dev は editable、コンテナは毎回 wheel をビルドするため、更新のたびに
  version を上げる必要はありません（GitHub には公開済み。将来 タグ固定 `git+https://…@<tag>` や PyPI 公開など
  **版固定の消費方式**へ移行する際に version 運用を切り替え）。

## Docker（任意）

本リポジトリにも `Dockerfile` がありますが、`ddbj_mss_tools` の slim イメージへ組み込む場合は、
`ddbj_mss_tools` 側の `scripts/build-ddbj-gff-wheel.sh` がこのリポジトリから wheel を生成します
（`docs/mss-tools-integration.md` 参照）。
