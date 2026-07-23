# ddbj_mss_tools からの利用と開発時の運用

`ddbj-gff`（このリポジトリ）は、`ddbj_mss_tools` の **`gff2mss` サブツール専用の依存**です。
`gff2mss` 以外のツール（egapx2mss / mss_builder / mss2ff / batch_wgs_builder）は `ddbj-gff` を
一切 import せず、無くても動作します。`ddbj_mss_tools` 側では `ddbj-gff` は**コア依存ではなく
optional extra**（`pip install "ddbj-mss-tools[gff2mss]"`）として扱われています。

`ddbj-gff` は開発中で頻繁に更新される想定のため、更新を各所へ反映する運用を以下にまとめます。
（現状 version は `0.1.0` 固定のまま運用。GitHub remote は `git@github.com:ddbj/gff_submission.git`。
PyPI 公開・タグ固定はまだ無し — dev は editable path、コンテナは wheel 直接指定で消費する。）

## 1. ローカル開発（推奨・最速）

`ddbj_mss_tools/pyproject.toml` に uv の path source が設定済み:

```toml
[tool.uv.sources]
ddbj-gff = { path = "../gff_submission", editable = true }
```

- **editable install** なので、この `gff_submission` を編集すれば `ddbj_mss_tools` 側に**即反映**（再インストール不要）。
- 前提: 2つのリポジトリが**隣り合って**チェックアウトされていること（`.../ddbj/gff_submission` と `.../ddbj/ddbj_mss_tools`）。
- uv 環境で `gff2mss` の開発・テストをする場合はこれが最も手間がかからない。

## 2. コンテナ / CI（`ddbj_mss_tools` の Dockerfile.slim）

`../gff_submission` は Docker のビルドコンテキスト外なので、**ビルド前に wheel を生成**して持ち込む:

```bash
# ddbj_mss_tools 側で
scripts/build-ddbj-gff-wheel.sh          # ../gff_submission から wheel を生成しリポジトリ直下へ
docker build -t ddbj-mss-tools:slim -f Dockerfile.slim .
```

- 生成される wheel（`ddbj_gff-*.whl`）は **`ddbj_mss_tools` にコミットしない**（`.gitignore` 済み）。
- **この `gff_submission` を更新したら、上記スクリプトを再実行してからコンテナを再ビルド**すれば最新が入る。
- version が `0.1.0` 固定でも、wheel をファイル指定で `pip install` するため毎回中身が入れ替わる（pip のバージョンキャッシュに引っかからない）。

## 3. 更新を反映する手順（まとめ）

| 反映先 | 手順 |
|---|---|
| uv ローカル開発 | 何もしない（editable path で自動反映） |
| コンテナ/CI | `ddbj_mss_tools` で `scripts/build-ddbj-gff-wheel.sh` → `docker build` |
| 非 uv の pip 環境 | `pip install -e ../gff_submission`（editable）または wheel を再ビルドして `pip install --force-reinstall <whl>` |

> 注意: `version` が `0.1.0` のまま `pip install ddbj-gff==0.1.0` 的な**版固定インストールをすると、更新しても
> 同一版扱いで反映されない**。上記はいずれも editable か wheel ファイル直接指定でこれを回避している。

## 4. 将来（安定・公開後）

`ddbj-gff` が安定し、GitHub 等に push できるようになったら、`ddbj_mss_tools` 側の依存を以下いずれかへ移行するとクリーン:

- `ddbj-gff @ git+https://github.com/ddbj/gff_submission.git@<tag>`（タグ/コミット固定）
- PyPI / プライベートインデックス公開 + 版固定

その際は**更新ごとに `version` を上げる**運用に切り替える（版でキャッシュが正しく無効化される）。
