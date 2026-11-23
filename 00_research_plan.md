# 研究方針：非可換 del Pezzo 曲面の非AS性・標準超曲面表示障害
2026-11-24 / non-as-v8。**研究再開時は本ファイルを最初に読む**。正本はこのGitHubリポジトリのmain。公表用論文は **`papers/nc-cubic-non-as.tex` 一本**、計算は `computations/` と `scripts/`、証明書は `certificates/`、研究ノートは `research/notes.tex`、文献は `references.bib`。旧リポジトリ・旧会話・旧bundleは再開に不要。

## 現行目標・価値基準
rank-one torsion-free構造対象、反SerreのAZ-ample性、反標準Hilbert次数を持つ内在的候補と、(I) 標準AS四次元ambientの三次超曲面表示、(II) Okawa--UedaのAS surface性との境界を明確にする。真の包含は `H_hyp∩AS_OU ⊊ H_hyp ⊊ I` として証明される。**最優先課題は内在的次数三候補にも属するOU三次型の包含 `AS_OU_3 ∩ I_3 ⊆ H_hyp_3` の決着である（標識を忘れたクラス）。** OUのDef. 1.1(2)そのままの型（箙の型）と、さらに内在的候補I_3に属するOU型を区別する。前者については弱del Pezzo境界構成による否定の証明を主論文 `thm:ou-weak-cubic`・`cor:ou-weak-unpointed` に収録した。後者は未解決であり引き続き最優先とする。数学的重要性、新規性、理解の深化、応用可能性、正確性、実現可能性を総合して最も価値ある主張に集中する。

## 証明を保持した主要結果
- **滑らかな楕円OU型の内在的認識**（`thm:ou-smooth-elliptic-recognition`）：admissible smooth elliptic decupleの標準noetherian qgr/導来helix比較を持つOU型では、内在的三次条件から、与えられた標識の反標準環そのものが標準AS4の正規三次商となる。OU spherical functorの反標準三角、曲線側自己同値のt-exact性、曲線heartのGK次元≤1を示し、AZ1/AZ2とHilbert差3nから制限対象が線束と従う。正則中心一次元の商は次数3 THCRであり、VdB §12の中心拡大判定を適用する。六点一般位置・非torsion性・係数の全非消失を仮定しない。admissibilityだけからI_3を主張しない。証明を主稿に収録、独立監査・新規性調査未完了。
- **六標識の一般位置境界**（`prop:ou-marking-walls`）：第一ブロックの6順序について各30条件、99 primitive行（81行は2標識、18行は1標識を妨げる）。3975の有理rank二行空間のうち全6標識を覆う30空間は全てδを含む。δ非torsionなら全失敗は余次元≥3。整数計算 `computations/verify_ou_decuple_markings.py` 実行済み。この結果はVdBの移動後六点一般位置のみで、AOU比較の全仮定とは区別する。
- **可換弱del Pezzo導来型の認識**（`thm:recognition-weak-commutative`）：I_3 の圏が滑らかな次数3のweak del Pezzo Xと導来同値なら、Xはdel Pezzoであり標識付きでcoh Xとなる。反標準環の中心上有限性、一般点のrank=1、torsion-free性・中心の正規性、特異点での無限Extによる証明。OU仮定は不要。原稿に証明あり、新規性・独立監査未確認。
- **AOUの既知の爆発範囲のheart認識**（`prop:helix-heart-upgrade`, `cor:aou-abelian-upgrade`）：helixの一対象が標準三次商FのπFへ移る導来同値は、Serre移動・very strong性・高次Ext消滅・負の生成対象を用いてt-exact。AOU Thm.9.4の比較が適用される範囲に標準三次超曲面のアーベル圏表示を得る。与えられた任意の標識の一致までは主張しない。
- **楕円データの逆変換と境界**（`prop:decuple-inverse`, `prop:decuple-zero-coefficients`, `rem:decuple-distinct-boundaries`）：AOU (9.53) の L0,L1,六点を正規化decupleから復元。27の関係係数消失条件を明示し、非torsionのL1/L0ではVdBの爆発一般位置条件と一致しないことを示す。指数ベクトルの代入検算済み。非可換境界での内在条件は未解決。
- **OU三次型の逆包含反例（新規性未確認）**：弱del Pezzoの反標準収縮上のheart `coh(Y,f_*End(T))` はOU三次型。全反標準AZ-ample対象のvirtual rank r≥2、反標準Hilbert多項式の先頭係数3r²/2≥6。従って標準AS4三次超曲面とunpointed非同値。ただしI_3にも属さない。原稿に証明あり、全文の独立監査・組版・新規性確認は未了。
- **次数3・量子Fermat**（元`nc-cubic-non-as`）：中心正規三次商の反標準圏について、離散torsion付き`(μ_3)^2`商への導来Morita同値と `HH_{-1,0,1}=(2,5,2)`。OUのAS surfaceではないが、AS四次元ambientの超曲面にはなっている。計算：`computations/verify_quantum_fermat*.py`。
- **次数3・非CY符号型二次元族**（元`nc-cubic-non-as`）：非中心正規三次商、退化twisted group algebra `M_2 ⊕ M_2` のMorita還元、`HH=(2,4,2)` とOU非AS性。fixed ambientでのpointed同型類は一般に二次元。H3認識は既存分類稿の結果に依存。計算：`computations/verify_sign_family_sectors.py`。
- **次数3・Clifford order**（元`nonas-ncdp`）：`Z=F_1, L=C+2f, V=O^2⊕O(f)`。反標準環は次数一生成の左右noetherian domainで `(1-t^3)/(1-t)^4`、二次包絡の三次成分は **19** 次元。connected gradedな標準quadratic AS4 ambientの三次正規商とは、構造対象を忘れたアーベル圏としても同値でない。任意の非周期的Z-algebraへの拡張は未証明。証明書：`certificates/degree3_quadratic.json`。
- **次数2・Clifford order**：Hilbert級数 `(1-t^4)/((1-t)^3(1-t^2))`、標準最小生成次数(1,1,1,2)の四次商とのunpointed非同値。他方、最小生成次数(1,1,1)の別のAS4 ambientの中心四次商になる。`HH=(3,4,3)`。
- **次数1・Clifford order**：Hilbert級数 `(1-t^6)/((1-t)^2(1-t^2)(1-t^3))`、標準最小生成次数(1,1,2,3)の六次商とのunpointed非同値。他方、最小生成次数(1,1,3)の別のAS4 ambientの中心六次商になる。`HH=(5,3,5)`。

ここで次数は `dim R_n=1+en(n+1)/2` による。Cliffordの有理標準因子では `K_A^2=e/4` であり、混同しない。奇数次HHが非零であることはOU型への障害であり、任意のAS4超曲面の非存在を意味しない。

## 定義の境界とトーリック商スタック
- 内在的候補 `I_d` はExt-finite・Serre functor・s=S^{-1}[2]のt-exact性・GK rank一torsion-freeなO・AZ ampleness・反標準Hilbert級数の条件からなる。exceptionalityやsmooth proper性を加えた強い条件は次数三Clifford節で証明。
- `H_hyp` は標準四次元二次Koszul AS正則代数の正則正規三次商として反標準環が表示されるクラス。三段階の真の包含は `H_hyp∩AS_OU ⊊ H_hyp ⊊ I_3`。rank一とAS4超曲面表示を同時に証明した符号型族が第一の真性、次数三Cliffordが第二の真性の証人。量子FermatのGK rank一性はこの包含の証明に利用しない。OU型をさらにI_3に制限した場合の逆包含は未証明。
- 境界例 `[P^2/μ_3]` の対角作用は `(x:y:z)→(ζx:ζ^{-1}y:z)`。構造層はAZ-ampleでないが、その反標準環は `C[U,V,W,T]/(UVW−T^3)`。3つのA2特異点を持つ。正則表現で作るrank三対象EはAZ-ampleだがEnd(E)=C^3、Hilbert級数は異なる。証明は主論文のRemark `rem:mu3-stack` に記載。

## 二つの精密化定理（公表論文に証明を収録）
- **二次包絡の必要十分認識**：Rが一次生成、Hilbert級数が標準三次型なら、Q(R)がnoetherianな14641型AS4正則環であることと、R=Q(R)/(g)（正則正規三次元）が同値。次数3と4の核の次元1,4からgB1=B1gが従う。整域性にはLevasseur (1992), Thm. 6.2(ii) を使用。これは既存分類稿と共有する結果であり、初出の新規主張として扱わない。
- **商スタックのrank下界**：X=[P²/μ3]、s=⊗ω_X^{-1}では、AZ-ampleなtorsion-free Fのrankは必ず3以上。反射的双対F**に延長されるHomと、固定点での慣性群の三つの指標で証明。rank3の正則表現が下界を達成。構造層のみの非amplenessより強い。

## 主要な未解決問題と優先課題
0. **最優先：`AS_OU_3 ∩ I_3 ⊆ H_hyp_3` の決着。** 可換弱del Pezzo導来型、標準helix比較を持つ滑らかな楕円decuple型、およびAOU Thm.9.4の既知の比較範囲では肯定的認識を得た。次の中心は滑らかな楕円decupleを持たないOU型である。(a) OUの三次型の定義とAOU Mrel/Mχ上の楕円表示の像・補集合を調べ、存在する退化型を分類する。(b) integral nodal/cuspidal cubicおよび可約・非被約反標準曲線について、spherical functor、heartへの制限、正則一次元、THCRの二次包絡判定を延長できるか調べる。特異曲線のPerfとDb(coh)を混同せず、点対象によるt-exact性の証明を再構築する。(c) 曲線表示自体を持たない型について内在的ample条件から反標準divisorまたは正規一次元が得られるか検討する。標識付き認識とunpointedの反例を区別する。literal OU型の弱del Pezzo例はI_3外である。滑らかな楕円認識の標準導来比較、GK次元、中心拡大の外部入力の独立監査を並行する。
1. **公表前監査と先行研究**：OU (arXiv:2508.04986v1) Def. 1.1, §3 と Levasseur (1992) Thm. 6.2(ii) を一次資料で確認済み。他の引用と独立検証を続け、OU一般とAS4超曲面の関係、rank-three stack型境界の他の例を整理し、主論文の全引用について一次資料の指定定理・仮定・版を確認する。特に離散torsion付きMorita同値とorbifold HKRの指標、Artin--Zhangの原著、Clifford orderの非表示の初出、新規性は独立検証未完了。
2. **次数3の境界拡張**：一般の非周期的AS正則Z-algebraにおいて、自由対象のSerre helixと三次合成空間が20次元となるための条件を調べる。Mori--Nymanの訂正版の仮定を保持し、非存在を一般に断定しない。
3. **符号型の精密化**：H3認識と固定ambientの二次元分類の既存分類稿との仮定を再照合し、versal変形・全intrinsicモジュライ成分性・OU locusの開閉性を検討する。
4. **Clifford反標準環の高次関係**：次数3の二次包絡について次数4以降の核を計算し、超曲面でない現象の構造的説明を求める。現行の非表示定理の成立には不要。
新しい証明が得られた場合は一本の論文の適切な位置へ組み込む。経緯や版別補遺を蓄積しない。

## 主張と検証の状態
現行の滑らかな楕円認識の証明と六標識の境界命題を主論文の `papers/sections/cubic_elliptic_recognition.tex`、`cubic_elliptic_markings.tex` に統合した。OU arXiv:2007.07620v1 (3)–(5), Prop.2, Thm.5、Anno–Logvinenko arXiv:1309.5035 (5.12), Def.5.2、AOU §7, §9、VdB arXiv:math/9809116v1 §11, §12の該当箇所を一次資料で確認した。ATV1990の楕円三つ組とLBSVdB1996 Cor.2.7の適用はVdB §12の明示した定式化を参照し、両原著の全文照合は未実施。六標識の99行・3975行空間・30被覆空間の整数検算をPython標準ライブラリで再実行した。新規の圏論的証明の独立監査・形式化・網羅的新規性確認、旧計算の一括再実行は未実施。以下は保持した従来の検証範囲。
現行の可換認識、heart認識、逆変換および関係係数条件は主論文に証明を収録した。OU arXiv:2508.04986v1、AOU arXiv:2404.00175v1 Lemma 9.2, Theorem 9.4と証明、VdB arXiv:math/9809116v1 §§11–12、CVGL arXiv:math/0702442 §1 の関係箇所を参照。楕円逆変換の8つの形式的Picard群生成元と27条件を整数指数ベクトルで代入し一致を確認した。今回の組版、全文独立監査、形式化、網羅的新規性調査、既存計算の全再実行は未実施。以下の従来の計算・検証記述は従来の確認範囲を示す。
既存結果に加え、OU三次型の弱del Pezzo境界例の証明を主論文に収録した。OU arXiv:2508.04986v1 Definition 1.1, §§3,5 とAOU arXiv:2404.00175v1 §§4,7,9 の該当箇所を今回確認した。新しい証明は独立監査・形式化・組版未実施。 今回の整数交点計算では L²=3, C²=-2, L·C=0, T_i·L=(0,0,0,1,1,1,2,2,2), T_i·C=(0,1,1,1,1,0,0,1,1) をPythonで再計算・一致確認。TeXの新規ラベルの一意性、proof環境の対応、参照する既存双対性補題を確認した。実環境にはLuaLaTeX/latexmkがあるがGitHubへのDNS接続不能で全repoを取得できず、主稿全体の組版・既存有限計算の再実行は未実施。有限行列、Groebner基底、Hochschild sector のコードは計算部分の検証に限る。前稿記載の数値と文献記述を無断で新規性確定に読み替えない。文献の未確認箇所や証明の追加仮定は現行ノートに明記する。

## 恒常的研究原則と保存許可
1. 研究の価値を最大化するため、重要性・新規性・理解・応用・正確性・実現可能性に従い方針を柔軟に更新する。変更理由と元目標の解決・残存部分を現在形で記す。時系列履歴を残さない。
2. 証明済み、追加仮定付き、未完成、予想、計算的観察を区別し、仮定・依存する結果・証明箇所を明示する。原稿に証明があること、今回精査したこと、計算実行したことを混同しない。誤りと依存主張を訂正し、有用な反例・障害のみノートへ残す。
3. 公表論文では重要主張の一次文献を調査し、仮定と既存結果・本研究の関係をIntroductionで具体的に説明する。未確認の書誌情報を補完しない。必要な識別子とリンクを保持し、外部PDFの大量収録をしない。
4. 論文は日本語LuaLaTeXを基本とし、公表価値がある現在有効な成果・証明を削除せず、論理的に一本へ統合する。研究ノートは研究日誌とせず、再取得困難な着想・反例・計算・文献比較を保存する。
5. この方針を再開入口とし、現行の論文、ノート、証明書、計算・再生成手順だけから研究を再開できるようにする。旧版・会話・別保存のノートを必要としない。
6. bundleは現在のソースと必要なPDF・データ・コード・証明書のみを含む。過去ZIP、編集履歴、ログ、重複、キャッシュを除去する。数学的証明・検証可能性を犠牲にした容量制限は設けない。
7. **GitHub mainを唯一の正本とし、追加承諾なく通常コミットで保存する。** 定常的なLibrary二重保存・ZIP作成は不要。push不能時の緊急保存または明示指示時にLibraryを用いる。対象repoの原稿・ノート・方針・READMEを保存してからresearch-hubの該当項目を最新内容に更新し、別々に成否を確認する。既存Libraryファイルは明示指示なしに削除しない。生成、検証、GitHub push、Library保存を別状態として報告する。PR・Actions・force pushは用いずmainへfast-forwardで反映する。

## 再現
```sh
python3 -m pip install -r requirements.txt
python3 scripts/verify_all.py
latexmk papers/nc-cubic-non-as.tex
python3 scripts/make_bundle.py
```
原稿はLuaLaTeX、BibTeX、latexmkで組版する。PDFとSyncTeXは`build/`に出力しGitで追跡しない。実際に実行した検証だけを報告する。
