# 設計書: 正規化器 特殊ケース①（recoded_codon / anticodon 子feature化 ＋ transl_except 翻訳）

- 日付: 2026-07-01
- 対象: フェーズ3-B-full の第1サブプロジェクト。CDS の `transl_except` 属性 → `recoded_codon`/`stop_codon` 子feature、tRNA の `anticodon` 属性 → `anticodon` 子feature（GFF→GFF, 3-B）。加えて Phase 2 CDS 翻訳が翻訳例外を適用するよう修正（テスト済み関数を vendoring）。
- 前提: Phase 3-B（common-case MVP）＋ feature-type fix が `main` にある。INSDC canonical 形は INSDC GFF3 Specification v0.5 に準拠。
- スコープ外（次サブプロジェクト）: trans-splicing（`location=join`）・circular 座標（end>seqlen）。プレースホルダ qualifier 能動処理も外。

---

## 1. 背景・INSDC canonical 形（spec v0.5 準拠）

3-A は `noncanonical-special-case`(INFO) として transl_except 属性 / anticodon 属性 / 非canonical trans-splicing を検出する（detect-only）。本サブプロジェクトは前2者の canonical 化を担う（1:1 対応）。spec v0.5 の canonical 形:

- **transl_except**（§218-234）: 非canonical = CDS の `/transl_except=(pos:<location>,aa:<amino_acid>)` 属性。canonical = CDS の**子feature**。col3 は `recoded_codon`(SO:0000145, 別アミノ酸へ) または `stop_codon`(SO:0000319, 終止へ・位置は CDS 末端)。recoded_codon には `codon_redefined`(SO:0000882) 属性=アミノ酸（例: `codon_redefined=selenocysteine`）。子は CDS 境界内。codon が2 CDS exon を跨ぐ稀ケースは同ID 2行。
- **anticodon**（§236-247）: 非canonical = tRNA の `anticodon=(pos:<location>,aa:<amino_acid>,seq:<text>)` 属性。canonical = tRNA の**子feature** col3=`anticodon`、`amino_acid=<aa>`・`sequence=<seq>` 属性、tRNA exon 境界内。跨ぎは同ID 2行。
- pos の座標は `139..141` または `complement(...)`（例: ecoli `transl_except=(pos:complement(4085235..4085237),aa:Sec)`）。URL-encode（`%2C`/`%28`/`%29`/`%3A`）され得る。
- aa は略号（`Sec`/`Glu`/`Term`）。canonical 例は full name（selenocysteine/glutamine）→ 略号→full name の小マップが要る。

Phase 2 は現状 CDS を plain `Seq(coding).translate(table)` で翻訳（convert.py:154）し、transl_except を扱わないため ecoli のセレノシステイン遺伝子（fdnG/fdoG/fdhF）で `translation-internal-stop` を誤発火する。テスト済み関数 `translate_cds_with_transl_except`（`nigyta/translate_with_exception` @ d3c3822、NIG 自前コード・再利用承認済み）で修正する。

---

## 2. 確定した方針（意思決定ログ）

| # | 論点 | 決定 |
|---|---|---|
| R-D1 | スコープ | recoded_codon/anticodon 子feature化（3-B）＋ transl_except 翻訳（Phase2）を **end-to-end** で。trans-splicing/circular は次サブプロジェクト |
| R-D2 | 3-B パス | `pass_transl_except`・`pass_anticodon` を `passes.py` に追加、`ALL_PASSES` 末尾へ。stdlib のみ（pos は正規表現、子feature 生成も stdlib） |
| R-D3 | recoded vs stop | aa が終止（Term/`*`）→ `stop_codon` 子（CDS末端）、それ以外 → `recoded_codon` 子（`codon_redefined=<full name>`） |
| R-D4 | 多part/境界外 | pos が単一range前提。CDS/tRNA 境界外、または join/多part な pos（2 exon 跨ぎ）は誤変換せず `needs-manual` report |
| R-D5 | 翻訳関数 | `translate_with_transl_except.py` を `src/ddbj_gff/mss/translate.py` として verbatim vendoring（provenance ヘッダ付）。Phase2 build_cds_feature が呼ぶ |
| R-D6 | 翻訳例外ソース | Phase2 は transl_except 属性（raw）＋ recoded_codon/stop_codon 子（3-B後）の**両方**から例外を収集し関数へ供給 |
| R-D7 | aa マッピング | 略号↔full name↔1文字 を集約（vendored `_AA_3TO1` 3→1文字を土台、full name 小 dict 補完）。3-B の codon_redefined 生成と Phase2 再構成が共有 |

依存: normalize パスは stdlib のみ維持。Biopython 翻訳は Phase2（`mss`、既に Biopython 使用）に閉じる。

---

## 3. アーキテクチャ・コンポーネント

**変更ファイル:**

| ファイル | 変更 |
|---|---|
| `src/ddbj_gff/normalize/passes.py` | `_parse_pos_spec` ヘルパ、`_AA_NAME` マップ、`pass_transl_except`・`pass_anticodon` 追加 |
| `src/ddbj_gff/normalize/normalize.py` | `ALL_PASSES` に2パス追加。`_APPLIED` に `add-child-feature` を追加 |
| `src/ddbj_gff/mss/translate.py` | **新規** — vendored `translate_cds_with_transl_except`（provenance ヘッダ） |
| `src/ddbj_gff/mss/convert.py` | `build_cds_feature` の plain translate を vendored 関数呼出に置換、例外を属性＋子から収集 |
| `tests/test_normalize_pass_transl_except.py` / `_anticodon.py` / `test_mss_translate.py` / `test_mss_cds`（追記）/ `test_normalize_integration.py`（ecoli slow） | テスト |

**共有ヘルパ（passes.py, stdlib）:**
- `_parse_pos_spec(spec: str) -> dict` — `(pos:LOC,aa:AA[,seq:SEQ])`（URL-decode 込）から `start,end,strand,aa,seq` を正規表現抽出。`complement(...)` → strand `-`、それ以外 `+`。単一 range 前提、join/複数part は `None`（→呼び手が needs-manual）。
- `_AA_NAME: dict[str,str]` — 略号→full name（`Sec`→`selenocysteine`, `Pyl`→`pyrrolysine`, `Glu`→`glutamine`, … ＋ `Term`/`*` は終止マーカ）。未知は入力値そのまま。

**`pass_transl_except(doc, ctx) -> list[Change]`:** 各 `f.type=="CDS"` で `transl_except` 属性がある場合、各値を `_parse_pos_spec`。
- パース不能/多part → `Change("needs-manual", …)`、属性は残す。
- 終止 aa → `stop_codon` 子、それ以外 → `recoded_codon` 子（`codon_redefined=_AA_NAME[aa]`）。
- 子 `Feature(id=f"{cds.id}_recoded_{n}", source=cds.source, type=<recoded_codon|stop_codon>, spans=[Span(pos)], attributes={...}, parent_ids=[cds.id])`。CDS span 範囲内チェック（外なら needs-manual・子なし）。
- 成功時 `transl_except` 属性を CDS から除去、`Change("add-child-feature", …)`。子は `doc.features` 追加＋`doc.feature_index` 登録＋親 `children` へ追加。

**`pass_anticodon(doc, ctx) -> list[Change]`:** 各 `f.type=="tRNA"` で `anticodon` 属性 → `anticodon` 子（`amino_acid=<full name or aa>`, `sequence=<seq>`, `parent_ids=[trna.id]`, span=pos）、tRNA span 境界内、属性除去。境界外/多part → needs-manual。

**`normalize.py`:** `ALL_PASSES = [pass_directives, pass_so_terms, pass_transl_except, pass_anticodon]`。`_APPLIED` に `add-child-feature` を追加（applied 扱い）。deepcopy・非破壊は不変。

---

## 4. Phase 2 翻訳統合

- **vendoring:** `src/ddbj_gff/mss/translate.py` に `translate_with_transl_except.py` を verbatim 取込。冒頭に provenance（source URL・commit `d3c3822`・"NIG own code, reused with authorization"）。公開 API `translate_cds_with_transl_except(feature: SeqFeature, parent_seq, stop_symbol="*") -> Seq`。
- **`build_cds_feature`（convert.py）:** 現行 `str(Seq(coding_full).translate(table=table_id))`（line 154）を置換:
  1. 翻訳例外を収集 — CDS の `transl_except` 属性（あれば）＋ mRNA/CDS の `recoded_codon`/`stop_codon` 子feature（子 span→pos、`codon_redefined`→aa を `_AA_NAME` 逆引きで略号化）。
  2. Bio `SeqFeature`（location=CDS `to_biopython_location()`、qualifiers=transl_table/codon_start/transl_except[再構成list]）を組み、vendored 関数で翻訳。
  3. 得た protein で既存検証（internal-stop / no-start）を実施。関数が initiator を `M` 強制・末尾 stop 除去する挙動に整合（no-start 判定は生 first codon 基準に調整、internal-stop 判定は変更なし）。
- **効果:** セレノシステイン等が `U` になり `translation-internal-stop` 誤発火が解消。パイプライン（3B→Phase2）でも raw 直接でも正翻訳。

---

## 5. テスト戦略・スコープ境界

**テスト（TDD・コンテナ内）:**
1. `test_normalize_pass_transl_except.py` — `(pos:139..141,aa:Sec)`→recoded_codon 子（codon_redefined=selenocysteine, span 139..141）＋属性除去; `complement(100..102)`→strand `-`; `aa=Term`→stop_codon 子（CDS末端）; 境界外→needs-manual・子なし; 属性なし→no-op。
2. `test_normalize_pass_anticodon.py` — `(pos:complement(14710..14712),aa:Glu,seq:ttc)`→anticodon 子（amino_acid=glutamine, sequence=ttc, strand -）＋属性除去; なし→no-op。
3. `test_mss_translate.py` — vendored 関数 smoke（transl_except 付 CDS が Sec→`U`、internal stop なし）。
4. `test_mss_cds`（追記）— recoded_codon 子を持つ CDS を build_cds_feature で翻訳 → `translation-internal-stop` 不発火。
5. **slow 統合 ecoli** — normalize → (a) validate で transl_except 由来 `noncanonical-special-case` 消滅（recoded_codon 子化）、(b) MSS で fdnG/fdoG/fdhF に `translation-internal-stop` 不発火＋recoded_codon 子存在。anticodon は crafted フィクスチャ（ecoli に anticodon 属性なし）。

**受け入れ基準:** 新ユニット＋Phase2 翻訳テスト green、全体スイート回帰なし（既存 ecoli validate slow は **raw** ecoli 検証につき transl_except noncanonical を期待 → 不変）。

**スコープ境界:**
- 内: `pass_transl_except`・`pass_anticodon`・vendored 翻訳・Phase2 統合・aa マッピング集約。
- 外（次）: trans-splicing（`location=join(...)`）・circular 座標（end>seqlen・is_circular）／プレースホルダ qualifier 能動処理。
- 非責務: 3-A 検証・SO-term 正規化パスは変更しない。
