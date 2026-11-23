# 研究方針：非可換 del Pezzo 曲面の非AS性・標準超曲面表示障害
2026-11-24 / non-as-v9。**研究再開時は本ファイルを最初に読む**。正本はこのGitHubリポジトリのmain。公表用論文は **`papers/nc-cubic-non-as.tex` 一本**、計算は `computations/` と `scripts/`、証明書は `certificates/`、研究ノートは `research/notes.tex`、文献は `references.bib`。旧リポジトリ・旧会話・旧bundleは再開に不要。

## 現行目標・価値基準
rank-one torsion-free構造対象、反SerreのAZ-ample性、反標準Hilbert次数を持つ内在的候補と、(I) 標準AS四次元ambientの三次超曲面表示、(II) Okawa--UedaのAS surface性との境界を明確にする。真の包含は `H_hyp∩AS_OU ⊊ H_hyp ⊊ I` として証明される。**最優先課題は内在的次数三候補にも属するOU三次型の包含 `AS_OU_3 ∩ I_3 ⊆ H_hyp_3` の決着である（標識を忘れたクラス）。** OUのDef. 1.1(2)そのままの型（箙の型）と、さらに内在的候補I_3に属するOU型を区別する。前者については弱del Pezzo境界構成による否定の証明を主論文 `thm:ou-weak-cubic`・`cor:ou-weak-unpointed` に収録した。後者は未解決であり引き続き最優先とする。数学的重要性、新規性、理解の深化、応用可能性、正確性、実現可能性を総合して最も価値ある主張に集中する。

## 証明を保持した主要結果
- **既約種数一曲線型の内在的認識**（`thm:ou-integral-elliptic-recognition`）：滑らか・節点・尖点三次曲線上の次数0/1/2各3本の線束表示と標準noetherian qgr/導来helix比較を持つOU型では、I_3の与えられた標識の反標準環そのものが標準AS4正規三次商。Perf上のOU構成とDb(coh)への制限右随伴の拡張を区別し、Burban–Kreußler Cor.5.8と二段の評価核でhの偶数シフトを排除。K_0(Perf)のdeterminantでPic^0への作用を認識し、ATV regular triple条件を得る。六点一般位置・非torsion性は不要。曲線表示だけからI_3を主張しない。
- **AZ2の特異点局所自由性と被約可約曲線の条件付き認識**（`lem:ou-az2-singular-free`, `prop:ou-reduced-curve-rank`, `prop:ou-reduced-cubic-geometric`）：幾何的hに対するAZ2をO_E→k(p)に適用して特異点での自由性を証明。被約連結種数1曲線で幾何的h（全次数3）、曲線への導来随伴・反標準三角・exact conservative J（像GK次元≤1）があれば、制限層は線束。被約平面三次曲線でさらにhが指定twist積に由来すれば、全次数0線束のdeterminantにより成分置換を排除、AZ1で正多重次数を得てvery ample性・regular triple性を示し標準AS4三次商を認識する。可約の場合にこれらのh・Jの仮定を一般に導くことは未解決。
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
0. **最優先：`AS_OU_3 ∩ I_3 ⊆ H_hyp_3` の決着。** 可換弱del Pezzo導来型、標準helix比較を持つ既約種数一曲線型（滑らか・節点・尖点）およびAOU Thm.9.4の既知比較範囲では肯定的認識を得た。次は被約可約三次曲線（line+conicの横断/接触、3直線の三角/共点）の場合を優先する。(a) 9線束のspherical twist積から幾何的hを認識し、全成分の正次数を導くか、非幾何的/半ampleの場合を排除する。(b) 幾何的でampleなhと通常のhelix比較からJのheart認識を得る（既存補題の拡張）。幾何的h・exact conservative J・GK次元≤1が得られれば被約可約曲線の環論的認識は既に証明済み。(c) AOU Mrel/Mχで既約曲線表示の像と可約・非被約・非曲線型を特定し、残る型では反標準自然変換の存在から調べる。非被約曲線では特異点有限性を用いたAZ2補題をそのまま適用しない。標準qgr/導来比較は仮定として維持。外部入力の適用と新証明の独立監査も行う。literal OU型の弱del Pezzo例はI_3外である。
1. **公表前監査と先行研究**：OU (arXiv:2508.04986v1) Def. 1.1, §3 と Levasseur (1992) Thm. 6.2(ii) を一次資料で確認済み。他の引用と独立検証を続け、OU一般とAS4超曲面の関係、rank-three stack型境界の他の例を整理し、主論文の全引用について一次資料の指定定理・仮定・版を確認する。特に離散torsion付きMorita同値とorbifold HKRの指標、Artin--Zhangの原著、Clifford orderの非表示の初出、新規性は独立検証未完了。
2. **次数3の境界拡張**：一般の非周期的AS正則Z-algebraにおいて、自由対象のSerre helixと三次合成空間が20次元となるための条件を調べる。Mori--Nymanの訂正版の仮定を保持し、非存在を一般に断定しない。
3. **符号型の精密化**：H3認識と固定ambientの二次元分類の既存分類稿との仮定を再照合し、versal変形・全intrinsicモジュライ成分性・OU locusの開閉性を検討する。
4. **Clifford反標準環の高次関係**：次数3の二次包絡について次数4以降の核を計算し、超曲面でない現象の構造的説明を求める。現行の非表示定理の成立には不要。
新しい証明が得られた場合は一本の論文の適切な位置へ組み込む。経緯や版別補遺を蓄積しない。

## 主張と検証の状態
現行版の既約種数一曲線認識・AZ2局所自由性・被約可約曲線の条件付き認識は主論文 `papers/sections/cubic_elliptic_recognition.tex` に証明を収録。OU arXiv:2007.07620v1 pp.1–2, (3)–(5), Prop.2, Thm.5 と Burban–Kreußler arXiv:math/0503496v2 Cor.5.8（p.37、原PDFを目視）、Prop.4.13、VdB arXiv:math/9809116v1 §12 p.121（中心拡大判定の明示的定式化）を今回再確認した。ATV (4.5), (4.8), Thms.6.6/6.8 は原著の公開テキスト抜粋と出版社書誌を確認し、原版PDF全体は未照合。特異曲線のPerfとDb(coh)を区別して証明したが、全圏論的証明の独立監査・網羅的新規性調査・形式化は未実施。OU 2025論文の今回のPDF再取得は失敗し、通常のnoetherian qgr/導来helix比較は引き続き仮定とする。
六標識の99 primitive行・3975有理二次元行空間・30被覆空間の整数計算と、逆変換・27係数条件の検算は従来の実施記録を保持する。今回の全計算再実行とは区別する。弱del Pezzoの整数交点計算、Clifford関係・判別曲線・Hochschild sector の実施記録は研究ノート・再現コードを参照。LuaLaTeX/latexmkは実行環境にあるが、今回もGitHubへのDNS接続ができず全repoの取得と現行主稿の組版は未実施。形式証明器による確認は行っていない。

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
