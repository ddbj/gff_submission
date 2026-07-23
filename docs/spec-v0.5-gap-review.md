# INSDC GFF3 v0.5 仕様 × 実装 レビュー

- **レビュー日**: 2026-07-02
- **対象実装**: `src/ddbj_gff/`(parser / model / writer / attributes / validate / normalize / mss)
- **突き合わせた資料**:
  1. `docs/INSDC GFF3 Specification - v0.5.docx`（本プロファイル仕様）
  2. SO GFF3 仕様 <https://github.com/the-sequence-ontology/specifications/blob/master/gff3.md>（base 仕様）
  3. 実装コードおよび `tests/` の fixture
- **備考**: 本書は検討用のドラフトであり、以降の議論で加筆・修正される前提。

`file:line` は該当コード箇所。仕様の行番号は docx をテキスト抽出した際の段落番号（おおよその位置の目安）。

---

## 0. 結論(最重要の 4 点)

1. **MSS 変換器が `gene→mRNA→CDS` 階層しか処理できない** — v0.5 の canonical example は CDS を `gene` 直下に置いており、原核生物 GFF も `gene→CDS`（mRNA なし）。どちらも現状は変換で捨てられる。→ A-1
2. **`is_circular` の大文字小文字不一致** — 実装は SO 流の `Is_circular` を読むが、v0.5 は小文字 `is_circular=true` を規定。conformant な環状フラグが無視され、origin-spanning feature に誤 ERROR が出る。→ A-2
3. **`location=join(...)` を誰も読んでいない** — trans-splicing / remote location の実座標はこの属性に入るが、パーサも変換器も col4/5 だけを見る。rice_cp（明示的な対象）で座標が壊れる。→ A-3
4. **仕様書そのものに正規表現・記述の不整合が複数あり**、実装以前に著者側へフィードバックすべき。→ E

---

## A. 実装と仕様の乖離(バグ / 要修正)

### A-1. CDS の親が gene のとき変換できない 〔重大〕

`mss/convert.py` の `build_gene_features`（277–321）は CDS を `_STRUCTURAL`（241 行）に含め、`collect_spans(mrna, "CDS")`（42–47）で **mRNA の子** から CDS を集める。ところが:

- v0.5 canonical example は `CDS ... Parent=gene_FABP5` = **gene の子**（さらに UTR は `exon` の子）。
- 原核生物 NCBI GFF3 は `gene→CDS`（mRNA なし）が標準。

前者では `collect_spans(mrna,"CDS")` が空 → `build_cds_feature` が `no-cds` を出して CDS を落とす。後者は `transcripts` も `rna_children` も空になり `no-rna` で gene ごとスキップ（280–284 行）。テスト（`test_mss_cds.py` の `mrna_with_cds`）も CDS を必ず mRNA 直下に付けており、この階層前提が固定化されている。**仕様の canonical example すら通らない**点で優先度が高い。

### A-2. `is_circular` の casing 〔重大〕

`model.py:102` は `self._first("Is_circular")`（大文字 I、SO 流）。v0.5 は小文字 `is_circular=true` を明記。結果:

- INSDC 準拠の環状フラグが `False` 扱い。
- `rules.py:52` の `rule_seqid_bounds` が `circular` 分岐に入らず、`end > seqlen` の origin-spanning feature に `feature-outside-region`（ERROR）を誤発報。

対処は両綴りを受理すること。なお仕様側も SO（`Is_circular`、landmark に付与）から逸脱し、INSDC では個々の origin-spanning feature に小文字で付与する点も注記が必要。

### A-3. `location=join(...)` / remote location 非対応 〔重大〕

trans-splicing の正準表現は `exception=trans_splicing` + `location=join(...INSDC座標...)` で、remote seqid（例 `PV366312.1:...`）まで含む。しかし:

- パーサは `location` を通常属性として素通し。
- `convert.build_insdc_location`（55–74）は col4/5 の span だけで座標を組む。

→ trans-spliced / remote feature は座標が誤変換される。`rule_special_case`（rules.py:144–159）は「trans_spliced なのに `location=` が無い」ことを INFO で指摘するだけで、`location=` の中身は使わない。

### A-4. pseudogene CDS に transl_table を必須化 〔中〕

`rules.py:106–120` は **全 CDS** に transl_table を要求（なければ ERROR）。仕様は「pseudogene CDS には省略すべき（should be omitted）」。`feature-mapping.tsv` の `pseudogenic_CDS→CDS /pseudo` を CDS として書いた擬遺伝子で誤 ERROR。`pseudo`/`pseudogene` 属性がある CDS は除外すべき。

### A-5. `_canonical_sort` が親→子順を保証しない 〔中〕

`writer.py:42–48` は `(seqid, start, end)` のみでソート。仕様は「親 gene は子より前に（子が上流にあっても）」を要求。gene と mRNA が同じ start のとき、mRNA は end が小さいので **gene より前** に並び、canonical order 違反になる（canonical example も gene→mRNA の順）。階層を尊重したソート（親を先に emit）が必要。

### A-6. multi-parent の重複行マージで 2 つ目の Parent が消える 〔中 / 仕様も曖昧〕

`parser.py:105–131` の `_add_row` は同一 ID の行を 1 feature にまとめ span を append、属性/Parent は**最初の行のみ保持**（差異は `attr-mismatch` WARNING）。ところが仕様は「複数 Parent を持つ feature は Parent ごとに重複行（同一 ID）で表現せよ」。したがって「2 transcript が共有する exon（同 ID・別 Parent）」を書くと、1 feature に span が二重登録され 2 つ目の Parent が失われる。**discontinuous feature（同 ID = 複数座標）の規則と、multi-parent 重複行の規則を混同**している構造的問題（→ C-2 の仕様曖昧さと表裏）。

### A-7. `partial` / `start_range` / `end_range` を消費していない 〔中〕

`convert.py` は部分性を UTR/開始・終止コドンの有無から**再推定**（`mrna_partial_flags` 108–121、`build_cds_feature` 168–169）するが、仕様が定義する明示属性 `partial=true` / `start_range` / `end_range` を読まない。投稿者が明示した部分性が失われ、推定と食い違う恐れ。validate 側にも `partial` 系のルールなし。

### A-8. `exception=ribosomal_slippage` 未対応 〔中〜低〕

validate / normalize / MSS のいずれも未処理。翻訳時のフレームシフトも適用されない（`translate.py` は transl_except のみ）。fixture `transl_except.gff3` に `exception=ribosomal slippage`（※スペース区切り）が既に存在。加えて `model.is_trans_spliced`（model.py:109–111）は `_`→`-` 置換のみで、スペース区切りの `exception` 値は拾えない点も注意。

---

## B. 改善点(仕様準拠だが質を上げられる)

- **col3 正規化で細粒度 SO 語を捨てている** — `passes.py:87–114` の `pass_so_terms` は `f.type` をコア型へリネームするだけで、仕様「細粒度の SO feature は col9 に入れる」を実施していない。元 SO 語を属性（例 `so_type=`）として残すべき。
- **`##sequence-region` を空白区切りで生成** — `passes.py:54` は `##sequence-region {seqid} 1 {length}`（空白）。仕様の正規表現はタブ区切りを要求（ただし仕様の canonical example 自体は空白 → 仕様側も不整合、E-1 参照）。ENA validator にかける前提ならタブが安全。
- **`##species` を `https` で生成** — `passes.py:10` の `_SPECIES_URL` は `https://`。仕様例・正規表現は `http://`。厳格 validator で弾かれ得る。
- **score の往復で `5`→`5.0`** — `writer.py:8` は `repr(float)`。整数値スコアの表記が変わる。原文保持なら整数はそのまま出したい。
- **protein_id / transcript_id を MSS が出力しない** — `build_cds_feature`（196–208）。DDBJ は protein_id を自付番するので意図的なら可だが、仕様では投稿者提供とされており、往復性の観点で注記が必要。
- **`###`（resolution boundary）は raw 保持のみ** — 前方参照解決の最適化には未使用（現状は許容）。
- **代表 transcript 選択が脆い** — `_representative_mrna`（convert.py:220–223）は ID が `.1` で終わるかで判定。命名規約に依存。

---

## C. 仕様に無い / 曖昧で、方針決定が要る点

- **C-1. part=# の整合検証がない** — 仕様は out-of-order interval に `part=#` を要求。パーサは読む（parser.py:94–97）が、「多座標 feature の全行が part を持つ/持たない」「番号が連番」といった検証ルールが無い。
- **C-2. 共有 exon の表現方法** — 「同 ID の discontinuous feature」と「Parent ごとの重複行」の切り分けが仕様上あいまい（A-6）。どちらを正準とするか決定的ルールが必要。
- **C-3. 環状トポロジの真正性** — 仕様「配列が circular と明示されていなければ origin-spanning 禁止」。トポロジは別配列ファイル側の情報で、GFF 単体では判定不能。`is_circular=true`（feature）と配列トポロジの突き合わせは未実装。sequence ファイル取り込みの設計が要る。
- **C-4. peptide-only FASTA の検証がない** — 仕様は核酸 FASTA を GFF に入れず peptide のみ許可。パーサは核酸/peptide を区別せず取り込み、writer は核酸 FASTA も再出力し得る。「FASTA は peptide のみ」を検証/警告するルールが無い。
- **C-5. `#!transl_table` のキー名** — 仕様の directive-value 欄は "nuclear"、usage example とコードは "primary"（`passes.py:66`）。どちらが正準か未確定。`primary:1,mitochondrion:2,...` の語彙定義も要る。
- **C-6. `insdc-gff-version` の既定値** — `config.py` は `"1.0.0"`。INSDC が宣言する現行版に合わせるべき（仕様 doc の "v0.5" とは別物）。
- **C-7. recoded_codon/stop_codon/anticodon 子の境界検証** — 仕様は「親 CDS/exon の範囲内に収まる」「終止コドンは CDS 末端」を要求。normalize は生成時にチェック（passes.py:158, 202）するが、**既に正準形で投稿された入力に対する validate ルールが無い**。

---

## D. SO GFF3 仕様の観点(base 仕様として遵守すべき点)

- **D-1. Parent の part_of 検証とサイクル検出が無い** — SO 仕様は「SO の part_of に無い Parent 関係は parse exception」「循環も exception」と明記。`_resolve_graph`（parser.py:134–146）は children/parents を張るだけで、型間関係もサイクルも検査しない。INSDC profile で緩めるとしても方針を明示すべき。
- **D-2. multiple parents** — SO は許容、INSDC は 1 行 1 親に制限。`rule_parents`（rules.py:92–103）は複数を ERROR で正しく制限（INSDC 準拠）。ただし A-6 のマージ問題は別途。
- **D-3. discontinuous features(同 ID = 複数座標)** — パーサの span マージで対応済み（良好）。
- **D-4. Target / Gap / Derives_from** — SO 予約属性だが未モデル化・未検証。INSDC 遺伝子アノテーションでは範囲外の可能性が高いが、alignment 系入力を受ける可能性があるなら判断が要る。
- **D-5. CDS phase は全 CDS で必須(SO)** — `rule_cds`（rules.py:117）が phase∈{0,1,2} を検査、良好。

---

## E. 仕様書ドキュメント自体の不具合(著者へフィードバック推奨)

- **E-1.** `##sequence-region` 正規表現 `^##sequence-region\t[^\t]+\t[0-9]+ \t[0-9]+` — 余分な空白（`[0-9]+ \t`）があり、かつタブ必須なのに canonical example は空白区切り。SO 例も空白。整合を取るべき。
- **E-2.** `#!transl_table` 正規表現 `...[A-Za-z0-9_]:[0-9]+...` — `[A-Za-z0-9_]` に量化子（`+`/`*`）が無く 1 文字しか一致しない。加えて "nuclear" vs "primary" の語彙不一致（C-5）。
- **E-3.** `##species` 正規表現にプレースホルダ文字列 `NCBI_Taxonomy_URI` がそのまま混入、二重空白もあり。
- **E-4.** `transcript_id` と `protein_id` の説明が同一（"This is the annotator provided protein identifier"）= コピペ誤り。
- **E-5.** 翻訳例外の例や anticodon 例で、行が改行なしで連結・タブが潰れている箇所がある（docx 由来の整形崩れ）。stop_codon 例が recoded 例と同じ ID `codon1` を使っている点も要確認。
- **E-6.** canonical example の `##gff-version 3.1.26` — SO 仕様版（1.26）を GFF version に埋め込んでいるように見える。SO 仕様は `3.#.#` を許すが意図確認が望ましい。

---

## F. 対応の優先度(提案)

| 優先 | 項目 | 種別 |
|---|---|---|
| 高 | A-1 CDS→gene / 原核 gene→CDS 対応 | 変換の適用範囲 |
| 高 | A-2 `is_circular` 小文字受理 | 環状ゲノム(細菌・plasmid) |
| 高 | A-3 `location=join` 解釈 | trans-splicing/remote(rice_cp) |
| 中 | A-4 pseudogene CDS 除外, A-5 親子順ソート, A-7 partial 属性消費 | 検証/変換の正確さ |
| 中 | C-7 / D-1 境界・part_of・サイクル検証 | validate 強化 |
| 低 | B 各項, A-8 ribosomal_slippage | 品質・網羅性 |

---

## G. 検証方法に関する注記

A-1 / A-2 / A-3 / A-6 は**コードの直接読解で確定**した挙動（`_resolve_graph`・`collect_spans`・`model.is_circular` の読解）。実行実証はこの環境に biopython が無く `pip install` が拒否されたため未実施だが、結論は確定している。着手時は、仕様の canonical example と rice_cp / 原核 fixture を回帰テストに追加してから修正する方針を推奨。
