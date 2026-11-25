# 研究方針：非可換 del Pezzo 曲面の非AS性・標準超曲面表示障害
2026-11-25 / non-as-v75（数学的内容の版）。正本は main。最初に本書を読む。

## 現行目標と到達点
現行の研究範囲はAOU三ブロック箙型。標識を忘れた包含 AS_3^(AOU) ∩ I_3 ⊆ H_hyp_3 は主稿 thm:ou-aou-all-unconditional に全支持・全特殊化の証明がある。標準noetherian OU比較を追加仮定としない。これは実現可能性や同型類の完全分類ではない。一般OU型への拡張は未解決であり、現行の限定目標と区別する。

## 主要結果と証明の所在
- 次数3の真の包含 H_hyp∩AS_OU ⊊ H_hyp ⊊ I_3：量子Fermat・非CY符号型の奇数次HH、次数3 Clifford orderの二次包絡の三次部分19（標準AS4は20）。構造対象を忘れた非表示も主稿で証明。非周期的Z-algebraへの拡張は未解決。
- 次数1・2 Clifford order：標準最小生成次数の重み付きambientではunpointed非表示、別のambientでは表示可能。主稿 weighted_* 節。次数は反標準Hilbert次数でありorderの標準因子の自己交点とは異なる。
- literal OU三次型の逆包含は弱del Pezzo反標準収縮による反例で否定（thm:ou-weak-cubic、cor:ou-weak-unpointed）。この例は I_3 外。
- AOU全係数認識：cubic_zero_supports.tex の cor:ou-only-one-zero-codim-five、thm:ou-one-zero-second-complete、prop:ou-one-zero-first-seven-orbits から thm:ou-aou-all-unconditional。九線束の全平面三次曲線認識（thm:ou-all-plane-cubic-recognition）、経路商rank、周期単純対象の厚化、標準比較の構成に依存。
- 一般OUへの部分結果：既約種数一曲線の球面的窓（thm:ou-any-quiver-elliptic-window）、固有DMスタックとの導来同値（thm:ou-derived-stack-recognition）、素・半素中心切断（thm:ou-prime-central-section、cor:ou-semiprime-central-section）、比較なしの半素PI/Goldie切断（cor:ou-literal-semiprime-pi-hypersurface）、平方零厚化（prop:ou-square-zero-depth）、中心有限な非被約切断（thm:ou-periodic-finite-curve-algebra、cor:ou-finite-central-nonreduced-cut）。各定理の仮定を保持する。
- 二次包絡の必要十分認識（既存分類稿と共有）、商スタック [P²/μ3] のAZ-ample torsion-free対象のrank下界3。主稿とノート参照。

## 次の課題（優先順）
1. AOU完全認識の独立監査：二族への帰着、第二族の全点判定、第一族七軌道、最後の二代表を順に照合する。周期対象の単純性、固定始域からの任意長商、AZ1/AZ2、Euler rank障害、右noetherian性、特異曲線のPerfとDbの区別を重点監査する。
2. 公表準備：Clifford非表示・弱del Pezzo反例・九線束認識の一次文献と新規性を確認。主稿のabstract・Introduction・結論を現在の定理に合わせて整理する。論文分割は別途判断し、現時点では一本を維持する。
3. 一般OU型への拡張を検討する場合は、中心一次元の作用・正則性、標準比較、非半素切断の根基の右有限生成性（cor:ou-pi-nil-finite-generation）を優先する。AOU外の全称命題は未証明。現行範囲の変更は明示する。
4. 非周期的Z-algebra ambient、符号型の相対変形とOU locus、Clifford二次包絡の高次関係は独立の残課題。

## 検証状態と再現
原稿に証明があることと独立監査済みであることを区別する。v75の記録では最後の十四グラフ、共有曲線外50閉点、追加ラベル30閉点、一変数表示29、例外既約因子延べ366の検算および主稿・ノート組版が報告されている。全runnerの一括成功、全証明の独立監査、網羅的新規性調査は未確認。今回の文書整理では再実行していない。
詳細な数学的知識・計算の参照は research/notes.tex と主稿、計算の正本は computations/、scripts/、certificates/。文献は references.bib。量子Fermat・符号型のH3認識の一部は nc-cubic-classification に依存する。OU・AOU・Artin–Zhang・Artin–Stafford・Levasseur・Opper・Clifford/order関連の引用の適用範囲を監査する。

```sh
python3 -m pip install -r requirements.txt
python3 computations/verify_ou_one_zero_first_square_difference.py
python3 computations/verify_ou_one_zero_first_final.py
python3 scripts/verify_all.py
python3 scripts/build.py --all
```
Python/SymPy、C++17、Singular 4.4（primdec.lib）、LuaLaTeX・BibTeX・latexmkが必要。PDFとSyncTeXはbuild/。旧会話・旧bundleは再開に不要。

## 恒常的研究原則と保存許可
1. 研究の価値を最大化するため、重要性・新規性・理解・応用・正確性・実現可能性に従い方針を柔軟に更新する。変更理由と元目標の解決・残存部分を現在形で記す。時系列履歴を残さない。
2. 証明済み、追加仮定付き、未完成、予想、計算的観察を区別し、仮定・依存する結果・証明箇所を明示する。原稿に証明があること、今回精査したこと、計算実行したことを混同しない。誤りと依存主張を訂正し、有用な反例・障害のみノートへ残す。
3. 公表論文では重要主張の一次文献を調査し、仮定と既存結果・本研究の関係をIntroductionで具体的に説明する。未確認の書誌情報を補完しない。必要な識別子とリンクを保持し、外部PDFの大量収録をしない。
4. 論文は日本語LuaLaTeXを基本とし、公表価値がある現在有効な成果・証明を削除せず、論理的に一本へ統合する。研究ノートは研究日誌とせず、再取得困難な着想・反例・計算・文献比較を保存する。
5. この方針を再開入口とし、現行の論文、ノート、証明書、計算・再生成手順だけから研究を再開できるようにする。旧版・会話・別保存のノートを必要としない。
6. bundleは現在のソースと必要なPDF・データ・コード・証明書のみを含む。過去ZIP、編集履歴、ログ、重複、キャッシュを除去する。数学的証明・検証可能性を犠牲にした容量制限は設けない。
7. **GitHub mainを唯一の正本とし、追加承諾なく通常コミットで保存する。** 定常的なLibrary二重保存・ZIP作成は不要。push不能時の緊急保存または明示指示時にLibraryを用いる。対象repoの原稿・ノート・方針・READMEを保存してからresearch-hubの該当項目を最新内容に更新し、別々に成否を確認する。既存Libraryファイルは明示指示なしに削除しない。生成、検証、GitHub push、Library保存を別状態として報告する。PR・Actions・force pushは用いずmainへfast-forwardで反映する。
