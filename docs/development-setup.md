# 開発環境の構築（別マシンへの移行手順）

`gff_submission`（Python パッケージ名 **`ddbj-gff`** / import 名 `ddbj_gff`）の開発環境を
新しいマシンで再構築するための手順。本リポジトリは GFF3 正規化・パース等を提供する**ライブラリ**で、
CLI は持ちません（`ddbj_mss_tools` の `gff2mss` サブツールから利用されます）。

## 前提

| 項目 | 要件 |
|---|---|
| Python | **3.11 以上**（`pyproject.toml` の `requires-python = ">=3.11"`） |
| uv | **推奨**（`uv.lock` で依存を固定管理している）。無い場合は pip でも可 |
| git | 必須 |

## 新マシンへの取得（重要: git remote が無い）

このリポジトリには**現在 git のリモートが設定されていません**。新マシンへ移すには、いずれか:

1. **ディレクトリごとコピー**（`.git` を含めて丸ごと転送。履歴も保持される） — 最も簡単
2. **リモートを設定して push**（GitHub 等に空リポジトリを作成 → `git remote add origin <URL>` → `git push -u origin main`）してから新マシンで clone

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
uv sync            # uv.lock に従って依存 + dev グループ(pytest) を .venv に構築
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
  version を上げる必要はありません（将来 PyPI/GitHub 公開・版固定へ移行する際に version 運用を切り替え）。

## Docker（任意）

本リポジトリにも `Dockerfile` がありますが、`ddbj_mss_tools` の slim イメージへ組み込む場合は、
`ddbj_mss_tools` 側の `scripts/build-ddbj-gff-wheel.sh` がこのリポジトリから wheel を生成します
（`docs/mss-tools-integration.md` 参照）。
