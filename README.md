# 非可換del Pezzo曲面の非AS性と超曲面表示の障害

**現行版 2026-11-25 / non-as-v45。** `shinnosukeokawa/nonas-ncdp` の次数1・2・3のClifford成果を統合し、主論文を **[papers/nc-cubic-non-as.tex](papers/nc-cubic-non-as.tex) 一本** にまとめています。

研究を再開する際は、最初に **[研究方針](00_research_plan.md)**、次に **[研究ノート](research/notes.tex)** を読み、必要な箇所を主論文の節ファイルで確認してください。旧リポジトリや旧版のbundleは必要ありません。

## Cloud Codexの作業順序

2026-11-25の実行方針：標準OU比較の仮定除去 → 1–3零型の完全認識・排除 → 共通依存証明の独立監査 → 独立公表可能な成果の確定。上位課題の障害で停止せず、次順位へ進み適宜戻る。詳細と保存手順は [研究方針](00_research_plan.md) および [AGENTS.md](AGENTS.md) を参照。15–30分を目安にmainへ通常コミット・pushし、research-hubを更新する。

## 箙を限定しない認識

内在的三次候補の圏が有限生成連結可換次数付き環の `qgr` に同値なら、任意の内在的標識は滑らかな三次曲面の標準標識に同値である。次数一生成性は不要である。より一般に、射影粗空間を持つ任意の複素固有Deligne–Mumfordスタックについて、内在的三次条件は滑らかさと全安定化群の自明性を強制する。次数一生成の左右noetherian整域でGK次元三、斉次商環の次数零部が曲面の可換関数体となる表示も、Rogalski–Staffordの分類から同じ結論に至る。主稿の `thm:ou-effective-orbifold-recognition` と続く二つの系に証明がある。さらに、有限型・固有で射影粗空間を持つDMスタックとの**導来同値**だけで、任意の内在的標識を認識できる。元の同値の t 完全性や反標準束の豊富性も仮定しない（主稿 `thm:ou-derived-stack-recognition`）。一般のOU圏にこれらの表示が得られるかは未解決である。

整な射影曲面上の連接代数の加群圏でも、一般点の代数が有限拡大体上の行列環なら同じ結論を得る。したがって、反標準環が有限生成の斉次中心部分代数上有限で、斉次商環の零次部が可換体なら、生成次数を問わず滑らかな三次曲面である（主稿 `thm:ou-split-order-recognition`、`cor:ou-finite-central-birational-recognition`）。一般点が非分裂な除法代数である場合はこの条件の外にある。

## AOU三ブロック箙型の全係数に対する条件付き判定

**現在の到達点。** AOU三ブロック箙型では、全27係数が非零の場合と三個以上が零の場合に加えて、二零台134217455で両方の零の補集合が二項式を満たす場合の全特殊化を処理した。指定零011の正規化では二番目の条件は `d=1` である。二つの零の交換により四族を三つの三変数表示にまとめ、23曲面、81曲線、340件の閉点因子を検証した。二零の残る四族でも60曲面のそれぞれで開判定を検証し、六次元係数空間内の残余を次元一以下に絞った（`cor:ou-two-zero-codim-five`）。一零では二族の内部に余次元四以上の残余がある。三零については十二の二変数族・62曲線・199件の閉点因子を既に検証した。主稿 `thm:ou-two-zero-double-binomial` と [再現コード](computations/verify_ou_two_zero_double_families.py)。

主稿の `cor:ou-aou-all-coefficients` で、標準noetherian OU比較が成立するAOU三ブロック箙型では、関係係数を問わず I₃ から標準AS四次元三次商の表示が従う。一・二・三零の補集合二項式外では追加の比較仮定を要しない。一・二零の五グラフ閉包の全特殊化の証明は `thm:ou-one-two-zero-all-boundaries`、検算は `python3 computations/verify_ou_one_two_zero_components.py` による。比較仮定の一般的除去と、AOU三ブロック箙以外を含むOU三次型の被覆は未解決である。

補集合二項式内の五台・十通りの指定零に現れる全43因子について、29件は三直線、14件は四成分の節点cycle上に、各因子の稠密開集合で九線束を構成した。四成分の場合は平面三角形の一節点とその例外曲線上の無限近点を含む六点爆発を参照し、標準射影の右noetherian性を導いた。1–3零型で比較仮定なしに残る係数は余次元三以上である。主稿 `thm:ou-complement-expanded-cycles`、検算 `python3 computations/verify_ou_complement_cycles.py`。

六個以上の零を持つAOU三ブロック型でも、標準比較を仮定せずに包含を証明した。残っていた七零の一点（台116309742、$u=1$）で、次数3と6の中心元が次数10の全594次元を生成し、巡回関係代数が中心上有限となる。主稿 `cor:ou-six-zero-automatic`、検算 `python3 computations/verify_ou_seven_zero_center.py`。

五零型でも全域で比較仮定を除去した。最後の正規化点134135534、$(u,v,w)=(1,1,1)$で、標数2の中心元との積が次数12の全819次元を張ることと、付値環を用いた右noetherian性の持上げを証明した。主稿 `thm:ou-five-zero-automatic`、検算 `python3 computations/verify_ou_seven_zero_center.py`。台132054766の補集合族は三直線（$a=-1$で共点）、134135534の補集合族は$a\ne1$で四成分の節点曲線を用いる。証明 `prop:ou-five-zero-cycle-pencils`、検算 `python3 computations/verify_ou_five_zero_cycle_pencils.py`。

四零の境界九行でも22件の残余因子（三成分6、四成分16）の一般点を処理した。これによりAOUの全係数台を通じて、比較仮定なしの残余が正規化係数空間の余次元三以上となる。主稿 `cor:ou-aou-codim-three`、検算 `python3 computations/verify_ou_four_zero_cycle_pencils.py`。

全AOU型のvery strong窓で矢印長3の中心が非零であり、一般点では次元1である。三次・四次経路空間が束となる開集合で、閉路と矢印の可換条件の核次元を特殊化した。特殊点での中心元の正則性は今後の課題である（`prop:ou-aou-central-specialization`）。

四零型も全特殊化で比較を仮定せずに包含を証明した。九つの境界内の曲線表示を二変数の局所化環上で確認し、そこから残る十三の一係数族には三次中心で割った後の頂点代数を使う。各頂点に三生成元のOre拡大の濾過商を構成し、右加群としての有限性と全ての例外根を検証した。AOU型で残るのは1–3零の余次元三以上の特殊化。主稿 `thm:ou-four-zero-automatic`、検算 `python3 computations/verify_ou_corner_comparison.py`。

一〜三零の未判定集合を、各正規化係数torusの余次元四以上に絞った。補集合上の三・四成分の曲線表示を係数条件の全域に拡張し、共点三直線も含めた。43因子内に残る余次元一成分162件には、三次中心商の九頂点それぞれで三生成元二次正則環を使う開条件を検証した。主稿 `prop:ou-complement-uniform-cycles`、`lem:ou-atv-open-corners`、`thm:ou-small-zero-codim-four`。再検算は `python3 computations/verify_ou_small_zero_corners.py`。

三零の台134184687は全特殊化で包含を証明した。この台の補集合二項式に残る三つの二変数族に対し、頂点代数の判定を曲線・閉点まで反復する。十五曲線の全例外根を検証した。当時残った三零台134151919と134217454に対しては、各42件を矢印スカラーと添字対称性で八族と四族にまとめ、全残余を覆う十二の二変数族を明示した（`prop:ou-three-zero-twelve-surfaces`）。主稿 `thm:ou-three-zero-double-pencils`、検算 `python3 computations/verify_ou_three_zero_double_pencils.py`。

## 標準比較と少数零型の進展（2026-11-25）

AOU三ブロック箙型のOU代数で、指定された九線束による平面三次曲線表示を持つ全場合と、全27係数が非零の全場合について、標準OU比較を仮定せずに認識する。有限次元加群を消した完全導来圏で窓を構成し、零次数成分の縮約とcornerの切断加群、対角三零の色交換関係を用いて標準射影の右noetherian性を証明する。二添字を共有する二零があると八次経路商が17以上（必要値15）となるため、追加の零を許して排除される。1–3零の13軌道中、五台がこの全域排除、完全matching一台も色交換と周期障害により全域でI₃外、七台は有理数上の21勾配idealにより各非空稠密開集合で比較仮定なしの認識に至る。補集合二項式の特異外側面も六次経路商の下界12（必要値10）により排除し、標準比較付きの二項式障害を超曲面全体に延長した。AOU三ブロック箙で関係係数の零がちょうど三個なら、標準OU比較の下で全九台・全特殊化について I₃ から標準AS四次元三次商の表示が従う。追加の比較仮定が残るのは三台の補集合二項式上だけである。一・二零では、補集合二項式の外に残る格子障害は正規化係数空間で余次元三以上となる。 一・二零の残余、曲線表示外の比較、AOU三ブロック以外を含む被覆は未解決。

証明は主稿の `prop:ou-perfect-tail-comparison`、`cor:ou-integral-automatic-comparison`、`thm:ou-plane-automatic-comparison`、`cor:ou-dense-automatic-comparison`、`thm:ou-two-collinear-zeros`、`prop:ou-few-zero-smooth-open` にある。再計算は `computations/verify_ou_collinear_zeros.py`（python-flint 0.9.0で高速化可能）と `computations/verify_ou_few_zero_cuts.py` による。

三零型の証明は主稿の `thm:ou-three-zero-complete-comparison`、七台の次元表は `prop:ou-few-zero-cyclic-codimension` にある。三巡回の同時失敗の素イデアル分解は `python3 computations/verify_ou_few_zero_rotations.py`、五つの経路商と二つの三直線族は `python3 computations/verify_ou_three_zero_boundaries.py` で再計算できる。前者には Singular 4.4 の `Singular` コマンドと標準の `primdec.lib` を用いる。実行ファイルを指定するには環境変数 `SINGULAR` を設定する。

補集合二項式を持つ五台・十通りの指定零のそれぞれで、その二項式の稠密開集合にも円錐曲線と直線上の九線束表示を構成した（主稿 `thm:ou-complement-line-conic`）。ここでは標準比較が自動となり、I₃から外れる。従って1–3零の各支持軌道で、比較仮定なしの包含が未判定なのは係数torusの余次元二以上の閉集合の中に限られる。開条件の非空性は `python3 computations/verify_ou_complement_line_conic.py` により検算できる。

一・二零の補集合二項式外では、三巡回の同時障害の75素成分を対称性で16軌道に分けた。61成分で経路商次元が必要値を超え、残り14成分の五つの代表族でも各稠密開集合に三直線上の九線束を構成した。二項式の外で未判定の余次元は四以上である。主稿 `thm:ou-one-two-zero-graphs` と `computations/verify_ou_one_two_zero_components.py`、`computations/verify_ou_one_two_zero_triangles.py` にグラフ、証明、再現手順を収録した。

## 内在的定義の三段階と境界例

次数3の内在的候補 `I_3`、標準AS四次元三次超曲面型 `H_hyp`、OUのAS surfaceでもある超曲面型 `H_hyp∩AS_OU` を区別し、真の包含 `H_hyp∩AS_OU ⊊ H_hyp ⊊ I_3` をIntroductionで証明しています（rank一の符号型族と次数三Clifford例を証人とします）。**OUの箙型による定義そのままでは逆包含は偽です。** 弱del Pezzoの反標準収縮からOU三次型のheartを構成し、全反標準ample対象の導来rankが2以上となるため標準三次超曲面と同値でないことを主論文 `thm:ou-weak-cubic` と `cor:ou-weak-unpointed` で証明します。ただしこの例は内在的候補I_3にも属さず、I_3に属するOU型での逆包含は未解決です。これを最優先課題とします。rank一・Hilbert級数を緩めた場合の境界例として `[P²/μ₃]` を計算し、rank三のAZ-ample対象と、構造層のnon-ampleness、特異な反標準三次模型 `UVW=T³` を区別しました。

## OU三次型の逆包含と反標準rank障害

滑らかな平面三次曲線上の6点を選び、p1,p2,p4だけを共線とする。爆発Xの反標準収縮f:X→YとAOUの9線束の三ブロックfoundation Tを用いる。対応するOU型は `coh(Y,f_*End(T))` であり、`coh X` 自体ではありません。型は逆Serre bimoduleのtensor代数表示から直接確認します。収縮曲線C上の `O_C`, `O_C(-1)[1]` はこのheartに入り、点対象を拡張として分解します。AZ2は反標準ample対象Pの局所射影性、AZ1は両対象への非零Homを保証し、導来rank r≥2、Hilbert多項式の先頭係数3r²/2≥6が従います。

今回、OU Definition 1.1, 3.1, 5.2、AOU §§4,7,9 の関連箇所を一次資料で確認し、新証明を主稿に収録しました。組版・全証明の独立監査・新規性確認は未実施です。

## 最優先課題：内在的条件を満たすOU型の認識

`AS_OU_3 ∩ I_3 ⊆ H_hyp_3` の決着を最優先とします。現行主稿では以下の範囲を明確にしました。

- 可換な次数3のweak del Pezzo Xと導来同値なI_3の圏は、標識付きで滑らかな三次曲面の連接層圏となる（`thm:recognition-weak-commutative`）。
- AOU Thm.9.4の指定された非可換爆発の比較が適用される範囲では、helixの一対象と自由対象との対応により導来同値がt-exactとなり、標識を忘れた標準AS4三次商の表示を得る（`prop:helix-heart-upgrade`, `cor:aou-abelian-upgrade`）。
- AOU (9.53) の楕円データを逆変換し、27の二次関係の消失係数を線束の等式で記述した（`prop:decuple-inverse`, `prop:decuple-zero-coefficients`）。非可換の場合、関係係数の全非消失とVdBの六点一般位置を混同できない。

- **既約種数一曲線型の内在的認識**（`thm:ou-integral-elliptic-recognition`）：滑らか・節点・尖点三次曲線上の指定9線束と標準OU helix比較を持つ型に対し、内在的標識の反標準環を標準AS4三次商と認識。PerfとDb(coh)を区別したspherical functor・Burban–Kreußler Cor.5.8・評価核の計算により幾何的hを得る。六点一般位置・非torsion性は不要。
- **AZ2と可約被約曲線への条件付き認識**（`lem:ou-az2-singular-free`, `prop:ou-reduced-curve-rank`, `prop:ou-reduced-cubic-geometric`）：幾何的hに対するAZ2から全特異点で局所自由、Hilbert差から制限線束を得る。被約平面三次曲線で幾何的h・exact conservative J（像GK次元≤1）・指定twist積を仮定すると、全成分の正多重次数とvery ample性、regular triple性により標準三次商を認識する。
- **可約被約三次曲線の多重次数認識**（`thm:ou-reduced-md-criterion`, `cor:ou-reduced-md-recognition`）：3ブロックの多重次数の和A,B,Cに対してδ_md=A+C−2B、w=5B−3A−2C。幾何的かつampleな反標準移動hの必要十分条件はδ_md=0, w_i>0。直線+円錐曲線（横断・接触）、三直線（三角・共点）を含む。2/3成分の6/9個の正かつδ_md=0の必要条件候補を列挙し、全ての古典参照型への巡回帰着を検算した。この条件と標準OU helix比較、内在的条件を持つ標識について標準AS4三次商表示を得る。
- **可約残余型の内在的排除**（`prop:ou-residual-delta-obstruction`, `thm:ou-residual-geometric-ample`）：標準OU helix比較と内在条件の下で、反標準自然変換εがheart上で非零ならδ_md=0（曲線側の幾何性・Jの完全性は不要）。幾何的hならwの負成分はvery strong性、零成分は有限長対象に無限個の相異なる単純商が生じることから排除。従って幾何的hは自動的にampleとなり標準三次商表示を得る。δ_md=0候補の16/168型と内訳を整数列挙で再確認。
- **節点cycle型の線形成長と反標準写像**（`thm:ou-nodal-cycle-growth`, `cor:ou-nodal-cycle-residual`）：横断する直線+円錐曲線、三角形の三直線では、導来helixのvery strong性だけから、noetherian比較・内在的標識を仮定せずhの写像類が有限冪でmultitwist、全ての固定Mへの全HomがO(n)、δ_md=0、全helix対象上でε≠0と示す。Opperのgentle algebra模型とNielsen–Thurston分解を用いる。I_3があれば構造対象上でεは単射。接触・共点型にも変形によりδ_md=0を移すが、εの非零性は未解決。
- **全被約節点cycle型の完全認識**（`thm:ou-nodal-complete-recognition`, `prop:ou-nodal-complete-word-certificate`）：三ブロック九線束の曲線を有限個の具体的写像類で同時に線束の曲線へ移す。Opperの曲線模型と導来自己同値群、自由群の共役類の等式、共通shiftの除去、共通零成分の縮約により、二成分16候補は12条件付き認識・4排除、三成分168候補は90条件付き認識・78排除、未判定0。標準OU helix比較と内在的条件の下で与えられた標識の反標準環そのものが標準AS4三次商。候補数は実現可能性や同型類数ではない。従来のgrading・multitwist精密化、幾何的共役、負指数排除も独立の結果として保持する。
- **加法型の代数的認識**（`thm:ou-additive-partial-recognition`）：単純rank二ベクトル束の評価・双対評価、spherical twistの九線束同時変換、一次jetによる共通零成分の縮約、一般位置を仮定しない古典参照型で、III型の16候補は12条件付き認識・4排除・未判定0、IV型の168候補は72条件付き認識・72排除・24未判定。未判定wは(-1,2,2)（18候補）と(0,0,3)（6候補）。条件付き認識は標準OU比較と内在的条件の下での帰着であり、曲面の存在や同型類数を意味しない。
- **IV型残余の関係代数の比較**（`prop:ou-iv-negative-relation-matches`, `prop:ou-iv-residual-relation-family`, `cor:ou-iv-six-residual`）：一次jetで切断・合成を明示し、負の多重次数の18候補の関係式を、許容される全gluing parameterでw=(1,1,1)の別表示と同定した（矢印の符号変更のみ）。IV型の最後の6候補を同じ五パラメータ関係代数族（六係数に一線形関係）へ帰着した。曲線側の自己同値による共役は主張せず、directed algebraの同型を用いる。係数空間と同型類moduli、線束表示とI_3の存在は区別する。
- **全被約平面三次曲線型の完全認識**（`thm:ou-iv-periodic-exclusion`, `thm:ou-reduced-complete-recognition`）：IV型の最後の関係代数族を、三つの反標準周期単純対象と固定helix対象の任意長の有限長商から排除した。AZ1・AZ2とEuler形式は内在的標識の数値rankを1とする一方、三周期対象へのEuler数から3以上を強制する。IV型の168必要多重次数候補は90条件付き認識・78排除、未判定0。指定の九線束表示と標準noetherian OU比較を持つ全被約平面三次曲線型について、内在的三次条件から標準AS4三次商表示を得る。全OU型の被覆と残存関係代数族のvery strong性自体は未解決。
- **全平面三次曲線表示の認識と行列式による復元**（`thm:ou-all-plane-cubic-recognition`, `thm:ou-smooth-determinant-reconstruction`, `thm:ou-integral-determinant-criterion`）：三重直線は全次数が3の倍数なので排除される。二重直線と直線の和は前向きExt消滅と箙の合成の生成性から九線束が二重直線上で共通捩りを除いて自明となり、尖点三次曲線へ縮約して同じdirected algebraを得る。従って任意の平面三次曲線の指定九線束表示、AOU三ブロック箙型の標準noetherian OU比較と内在的条件の下で標準AS4三次商表示を得る。一般の関係テンソルから三つの行列式の線形式により P⁵ の固定二項三次超曲面を切断する。3線形式が独立で切断が滑らかなら九線束表示を直接復元。既約特異切断では偶奇座標直線の九交点が滑らかで、特異点における三つの線束行列のrankが二であることが同じ頂点ラベルの既約曲線表示の必要十分条件。可約・非被約切断やrank落ち等の条件外は別問題。
- **全係数非零のrank2可換locusとrank3可約切断の認識**（`cor:ou-cyclic-support-obstruction`, `cor:ou-dense-commutative-rank-criterion`, `thm:ou-reduced-determinant-grid`）：三次syzygyから各slice比のrank1を排除（係数零の場合はsupport graph条件）。全係数非零・very strong導来helixの三つの巡回表示で、行列式rank2が一つでも現れることと滑らかな可換三次曲面のAOU九線束型であることは同値。rank2の切断は滑らかな三次曲面で、九直線と行列核を六点爆発の標準線束列に同定し、標準helix比較ならcoh Xと同値。その他は全巡回rank3、零平面は二項三次内に含まれない。rank3の係数による復元条件は(i)零平面上の偶奇座標形式の独立、(ii)係数行列の偶奇2列の独立、(iii)各kに対してあるk'≠kでdet(C_k'/C_k)≠0（成分ごとの比）。この三条件で被約曲線の九線束を復元（可約も含む、既約では必要十分）。格子点の成分は行だけまたは列だけで決まるため多重次数・核・直交性を制御できる。巡回で既約rank障害が可約表示に移り解消する具体例も記号検算したが、その例のvery strong性・I_3は未確認。
- **四零係数型の認識・排除の完了（標準比較付き）**（`thm:ou-four-zero-path-exclusion`, `prop:ou-four-zero-smooth-open`, `thm:ou-four-zero-grid-recognition`, `prop:ou-four-zero-rank-two-recognition`, `prop:ou-four-zero-mixed-grid`）：全30台の指数行列rank19、四自由係数。22台を長さ4/5/6/7/8の厳密rank証明書1/3/3/7/8台で全特殊化排除（全行列11、再帰商2、Laurent単項式・原始多項式消去9）。残余133136111,133660399,133889775,134135535,134147823,134151407,134151918,134183663。133136111を除く七台は正規化(2,3,5,7)で元の巡回の行列式rank3、滑らかな三次切断（21勾配idealの厳密計算）を持つため、各族の非空Zariski開集合上で標準比較付き認識。`thm:ou-four-zero-grid-recognition` で90局所化勾配idealと巡回を用い133660399,134183663を全域認識、133889775はv=q=1,u+w=1以外、134147823は(v−1)(q−1)≠0かつ(u+w,v−q)≠(0,0)、134151407はv≠1、134151918は(u−1)(q−1)≠0で認識。`prop:ou-four-zero-rank-two-recognition` で133136111のrank2 locus v=1,q=−1,u+w=1（零空間の三次曲面の特異点有限、Bertini）、134135535の全係数1（明示平面）も認識。133136111のrank3の全巡回は零列に対応する座標点で必ず特異。`lem:ou-two-zero-slice-grid` は空の部分matchingも逆数点が零平面外なら許す。`prop:ou-four-zero-mixed-grid` は134135535の各係数≠1で27勾配minor idealとtorus rank条件を判定、三matchingのうち二組等値・二組和1の両方を避ける巡回があれば認識。全巡回障害は同値四係数または相補的二値各二個に限る。`thm:ou-four-zero-triangle` により133889775の残余も三直線上（w=−1の共点を含む）で認識、全域認識完了。`thm:ou-four-zero-mixed-triangle` により134135535の等値族は共通値1/2以外を認識、相補値の三配置はδ_md=(0,1,−1)によりvery strong性外、全係数1/2も五次商8>7で排除。134135535の係数1を含む境界は以下の補集合障害で処理する。一般定理 `thm:ou-isolated-matching-obstruction` と `cor:ou-isolated-matching-supports` により孤立零の完全matchingを持つ10台を標準比較付きI_3外とした（既排除八台と新規の四零133136111・三零133947135）。奇置換三列Oは対角零・非対角非零でrank2または3、rank2でも右核の全座標非零と一変数形式陰関数定理で位数1の巡回を構成。133136111のrank2九線束表示の存在とは両立する。補集合の形式曲線は一般に行列式零平面の点(1,1,0,0,0,0)の局所枝を正規化し、rank1の2×2行列をt²基底変換後に両因子が正位数となるよう分解して構成できる（`thm:ou-complement-general-arc`）。四零の九境界を `thm:ou-four-zero-complement-boundaries` で標準比較付きI_3外：134135535の任意の係数1（全1も含む）、134147823のv=1又はq=1、134151407のv=1、134151918のu=1又はq=1かつv≠1。`thm:ou-four-zero-boundary-paths` は134147823のu=v−1,w=1−v,q=vを長さ6商10>9、134151918のu=v=1を長さ8商19>15、v=q=1を長さ7商14>12で排除。残る134147823 w=−u,q=v,v≠1,u≠v−1は `thm:ou-four-zero-last-conic` で直線と滑らかな円錐曲線の長さ二交叉上の九線束を明示復元（δ_md=0,w_md=(2,1)）。`cor:ou-four-zero-complete-recognition` により四個以上の零係数型を標準比較の下で認識・排除完了。残る1–3零・標準比較は未解決。台の実現可能性は未主張。
- **五零係数型の認識・排除の完了（標準比較付き）**（`thm:ou-five-zero-path-exclusion`, `lem:ou-diagonal-zero-arc-lifting`, `thm:ou-five-zero-periodic-exclusion`, `thm:ou-five-zero-grid-recognition`, `thm:ou-five-zero-triangle`, `thm:ou-five-zero-expanded-curves`, `thm:ou-five-zero-117374702-boundary`, `lem:ou-complement-periodic-simple`, `prop:ou-complement-arc-obstruction`, `thm:ou-five-zero-complement-exclusion`, `cor:ou-five-zero-complete-recognition`）：80台は全て三自由係数。71台を長さ4/5/6/7/8の厳密経路商証明書6/20/12/19/14台で排除（59台は全関係行列、12台は右矢印作用を保持する再帰商計算）。九台残余。125763310,132086511は対角三零と列1,2,5のminor（−uw−v,−w−1）が非零なら形式陰関数定理、零でも明示的な一変数形式方程式により全係数で標準比較付きI_3外。各sliceの零が空でない部分matchingなら核成分の大域非零性で被約格子復元できるよう `lem:ou-two-zero-slice-grid` を一般化。90格子勾配idealを厳密検算し132054255,133099247を全域、133873391は(1,1,1)以外、117374702はu≠1、132054766,134147822はw≠1、134135534はu,w≠1で認識。133873391の(1,1,1)も巡回c'_{ijk}=c_{kij}の三角形(x+y)(x−z)(y−z)=0上の明示的な九線束の貼り合わせ・核・合成比較で認識（δ_md=0,w=(1,1,1)）。従って同台は全域認識。四台の境界も直線・円錐曲線の長さ二交叉D上の明示的な六線束貼り合わせと核・合成比較で認識（δ_md=0,w=(3,0)、全平面三次曲線認識定理を適用）。117374702はu=1,v≠1、132054766はw=1,v≠−1、134147822はw=1全域、134135534はu=1,w≠1とw=1,u≠1（係数正規形の対称性(u,v,w)↦(w,1/v,u)）を追加。117374702 u=v=1は長さ七・巡回0・始点0終点2の1458×729関係行列のQ(w)上の厳密rank716により全域でvery strong性から排除（商次元13以上、必要値12）。残る132054766 v=−1,w=1、134135534 u=w=1も孤立零(011)と補集合の六頂点に支持を持つ周期単純対象の同時形式厚化で標準比較付きI_3外。補集合八係数が1、三つの外側2×2面が可逆なら補集合表現RはS(R)=R[2]で単純。明示A(t),B(t)の特殊化Q⊕Rは各頂点一次元、両headが有界なので固定始域の任意長商と数値的rank1≥2の矛盾を得る。全五零80台の認識・排除が完了（実現可能性は未主張）。残余九台は一有限体特殊化で長さ12まで必要値と一致するが実現可能性は未主張。標準比較とI_3を認識に仮定する。
- **六個以上の零係数型の認識・排除の完了（標準比較付き）**（`thm:ou-six-zero-path-exclusion`, `thm:ou-six-zero-periodic-exclusion`, `lem:ou-two-zero-slice-grid`, `thm:ou-six-zero-curve-recognition`, `thm:ou-six-zero-expanded-curve`, `cor:ou-six-zero-complete-recognition`）：201台中53/76/30/23/12台を長さ4/5/6/7/8の厳密経路商rankで排除（計194、有限体は候補端点のscreeningのみ）。残る七台は117358318,117370606,125759214,131005679,131037935,133099246,134145774（131037935だけ三係数、他は二係数）。131037935は既存周期障害、125759214は形式曲線 R=-v+t^3, D'=1−(v−1)R/v, h=(u+D')(R+D')/(RD')（u+v≠0、u=−v≠−1、(-1,1)の別表現を全て扱う）で標準比較付きI_3外。五族の二項三次切断をxの二次式として原始性・判別式平方条件・45偶奇対の勾配idealで厳密判定。117370606,133099246,131005679は全非零係数で既約、134145774は(-1,-1)だけ直線+円錐曲線だが各sliceの二零が異なる行列位置にある場合の新しい九線束復元補題で全域認識。117358318のu≠1は既約切断で認識、u=1もN=V(XZ−Y²)とL=V(X−Y+vZ)の長さ二の交叉D=s²−st+vt²での明示的な線束貼り合わせ・核・合成比較により全v≠0で認識（v=1/4の接触型を含む）。δ_md=0,w=(3,0)なので全平面三次曲線認識定理を適用。六零の未判定は零。五零の現状は上記の部分認識を参照。残る1–4零と標準比較。認識は標準比較・I_3付き、実現可能性は未主張。
- **七個以上の零係数型の認識・排除（標準比較付き）**（`prop:ou-zero-support-enumeration`, `thm:ou-large-zero-path-exclusion`, `thm:ou-seven-eight-zero-path-exclusion`, `thm:ou-latin-square-periodic-exclusion`, `cor:ou-intrinsic-at-most-seven-zeros`, `prop:ou-seven-zero-curve-recognition`）：9 sliceの連結性と異slice両順序の強制零を含まない木成分の排除で全台を列挙（零係数数0–12の同値類数1,1,3,9,30,80,201,344,387,217,53,4,1）。torus正規化・厳密経路rankにより10–12零の58台全部、9零の217台中216台、8零の387台中386台、7零の344台中336台を排除。8零の唯一の残余108979950は9零のLatin square型の零を一つ除いた台。対角三零・一方向巡回三零・ちょうど二添字が等しい18係数非零に周期単純対象と任意長有限長商による障害を一般化し、9零、8零の残余と7零の125234926をI_3から排除（標準比較付き）。従ってI_3内のvery strong OU型の零係数は高々7個。七零の116309742は自由係数u=c_211≠1で滑らかな切断から認識、116834030はu=c_222≠0の全域で認識（u=1/16の唯一の通常節点は六座標非零で核rank2）。定理 `thm:ou-seven-zero-line-conic-exclusion` は128957047,129972855,133118583をvery strong性から全係数で排除。128957047のu=1は六次商次元11（必要値10）、他は直線+滑らかな円錐曲線の切断と明示的な九線束の貼り合わせ、核、直交性・合成比較により δ_md=(-1,1)≠0 を得て既存の多重次数必要条件に反する（標準比較・I_3不要）。さらに `thm:ou-seven-zero-formal-arc-exclusion` で117368558を形式表現と有界head（u≠−1で位数1、u=−1で3）による固定対象の任意長有限長商・周期rank障害から全係数でI_3外とした（標準比較付き）。`thm:ou-seven-zero-eighth-path-exclusion` は125531767を八次経路商の有理関数体rank2170、商次元17>15から全係数でvery strong性外とした。最後の116309742のu=1も `thm:ou-seven-zero-six-cycle-exclusion` により標準比較の下でI_3外とした。六孤立零に対応するs固定単純対象の同時厚化を各頂点rank2の形式表現で構成。q=t^6、ラベル差の次数と三周期対称性から18方程式・30変数に帰着、対角式をqで割った初期値でJacobian小行列式−8を厳密確認。形式陰関数定理と位数1の二次minorにより固定始域の任意長有限長商を得て 2r([O])=2≥6 の矛盾。`cor:ou-seven-zero-complete-recognition` により七個以上の零係数型の認識・排除を完了（標準比較付き）。六零以下の現状は上の帰着定理を参照。台の実現可能性は未主張。
- **全係数非零型の認識（標準比較付き）**（`thm:ou-dense-residual-three-lines`, `thm:ou-dense-complete-recognition`）：残余五係数族 F1=0 の零平面はz0=z1=t、x+w=y+z、(D+S)t−Sx+(b−d)p+(c−d)q=0（p=y−x,q=z−x）となり三次式はtpq=0。S≠0は三角形、S=0は共点三直線。三成分上の射行列を特異点を含むxw≠0の共通開集合上の行列と単元の行列rescalingで貼り合わせ、共点の場合も線束を直接構成。核・余核は全点rank2から局所自由、九線束のmdegはU0=(1,0,−1), U1,U2=(0,0,0), V0=(0,1,0), V1,V2=(1,0,0), W0=(1,0,1), W1,W2=(1,1,0)。δ_md=0,w=(1,1,1)、前向きExt消滅・ブロック内完全直交・合成関係を証明。既存のrank2認識・rank3復元・残余への帰着と合わせて全係数非零かつvery strong型を尽くし、標準noetherian OU比較と内在的三次条件の下で標準AS4正規三次商として認識する。新規性・全証明の独立監査は未確認。一般の零係数型・標準比較は未解決。
- **座標従属の五係数族への帰着**（`lem:ou-dense-constant-face-obstruction`, `thm:ou-dense-five-coefficient-residual`）：rank3全係数非零の正規形M=[[a,b],[c,d]], N=λ1+μM（λ≠0）でvery strongならμ≠0, τ=1+λ−μ≠0。μ=0の定数面族の七次1458×729多項式行列（始点・終点添字1）の有理関数体rank716を実行確認、全特殊化で商次元13以上（必要値12）。τ=0は基準slice変更で同じ族。元の巡回が(ii),(iii)を満たすと他の二巡回も(ii),(iii)を満たし、(i)はF1 F2 F3≠0と同値。D=ad−bc, S=a+d−b−c, T=abc−abd−acd+bcd, F1=(μ−1)D+λS, F3=μT−(1+λ)D, F2=F3−F1。全巡回の失敗はkの添字変更でF1=0, D≠0, μ=1−λS/Dに帰着。標準比較・I_3不要。残余族は上の三直線表示により認識済み。
- **核rank障害の高々一巡回性・格子条件との同時成立・残る座標従属**（`thm:ou-dense-at-most-one-torus-failure`, `thm:ou-dense-common-grid-torus`, `prop:ou-dense-dependent-coordinate-normal-form`）：全係数非零・very strongなら(iii)に失敗する巡回は高々一つ。二方向のrank1面から非零辺ベクトルの二部graphを作り、非連結なら四次次元超過、連結ならd_ijk=m_ij z_kで五次双対の次元が8以上（必要値7）。混合障害は一方向の(iii)失敗と他二方向の(ii)失敗に帰着し、障害長方形の比較からdの二面 M=[[0,0],[0,q]], N=u[[1,w],[v,vw]] の族のみ残る。六次405×243多項式行列の有理関数体rank231から全特殊化で商次元12以上（必要値10）となり排除。従って(ii),(iii)の共通成功巡回が存在する。標準比較・I_3不要。rank3でなお(i)が失敗する場合は、従属する座標対の共有辺を(0,0)とし正規化して二面N=λ1+μM（M非定数、λ≠0）に帰着。この六係数族は上の帰着定理によりさらに絞り込まれる。
- **四次経路次元公式と三巡回同時格子障害の排除**（`prop:ou-dense-fourth-path`, `thm:ou-dense-no-triple-grid-failure`）：全係数非零の巡回二次関係代数で、端点(i,j)を含み交叉比が第3添字に依らない長方形数bに対し、長さ4の経路商次元はmax(5,4+b)。very strongなら各巡回の障害長方形は高々1。三巡回同時の(ii)失敗は立方体三辺の係数1への正規化で、二つの四係数族に帰着。食い違う軸が2以上の型は四次公式で排除。残りは有理関数体上の405×243（長さ6）・1458×729（長さ7）行列の厳密rank 231・716から、全特殊化で商次元12以上・13以上となりEulerの必要値10・12に反する。標準比較・I_3不要。`computations/verify_ou_grid_higher_syzygies.py` を実行し8立方体型、四次公式の有理数例と有理関数体rankを確認した。この格子障害排除は上の共通成功巡回定理の入力。
- **高次syzygyによる三巡回同時の核rank障害の排除**（`thm:ou-dense-no-triple-torus-failure`）：全係数非零かつvery strong導来helixなら少なくとも一つの巡回で(iii)が成立する。三方向の障害sliceを0とし係数をrescaleして全座標零面を1にすると、残る2×2×2立方体のc−1の各面は非零rank1。台は全立方体または対頂点のみ。積型では長さ5の経路商次元8（必要値7）、対頂点型では長さ4の経路商次元6以上（必要値5）となる。一般証明は主稿、標準比較・I_3・行列式rank3は不要。`computations/verify_ou_dense_higher_syzygies.py` で全256台と二つの有理数係数例を検算。上の高々一巡回性により強化される。
- **六標識の境界**（`prop:ou-marking-walls`）：第一ブロックの6順序でVdB一般位置を比較し、99の整数等式と3975のrank二行空間を得る。δ非torsionの場合に全標識が失敗する集合は余次元3以上。これはAOUの全仮定を回復する主張ではない。

一般の包含は未解決ですが、全係数非零と六個以上の零係数型は標準noetherian OU比較と内在的三次条件の下で認識・排除済みです。今回、五零80台中71台を厳密経路rankで排除し、二台を形式厚化による周期障害で排除、七族を開集合上で認識し、さらに133873391の(1,1,1)も三角形上の九線束表示により認識しました。今回、残る四台の境界を直線・円錐曲線の長さ二交叉上の九線束表示により認識し、未判定を三つの一係数族へ縮小しました。今回、四零の30台中22台を厳密rankで排除し、残余八族中七族で滑らかな切断の非空開集合を確認しました。今回さらに四零の二族を全域、133889775を一係数残余以外、他三族を明示開集合で認識し、二台のrank2 locusと134135535の混合sliceの認識範囲を得ました。五零の三残余のうち117374702 u=v=1を七次経路商の厳密rank716で排除しました。残る二つの一係数族も補集合の周期単純対象と明示形式厚化により標準比較付きI_3外とし、五零の認識・排除を完了しました。今回133889775の一係数残余を三直線上で認識し同台全域の認識を完了、134135535の等値族・相補値三配置を認識または多重次数障害・五次経路商で排除しました。今回、孤立零の完全matchingから内在的条件を排除する一般定理を得て、四零133136111と三零133947135を全域排除しました。以前の八台の周期障害も統一されます。`verify_ou_isolated_matching.py` を全実行し、全台列挙・10台のmatching・rank2奇置換行列と形式方程式の恒等式を厳密検算しました。四零の最後の境界も補集合の一般形式曲線、三つの厳密経路商と直線・円錐曲線の九線束表示で処理し、四個以上の零係数型の認識・排除を完了しました。次は1–3零と標準比較を優先します。今回 `verify_ou_four_zero_triangle.py` と `verify_ou_four_zero_mixed_triangle.py` の全実行で貼り合わせ・核・rank・合成・多重次数と全係数1/2の五次商を検証しました。今回 `verify_ou_four_zero_grid.py` を全実行し90格子勾配ideal、rank2 locus・滑らかな平面切断・27勾配minor ideal・torus逆数点を検証しました。前回の `verify_ou_four_zero_supports.py` の3区間全22証明書・30台照合・21勾配idealの検証は既存の記録です。前回の `verify_ou_five_zero_expanded.py` も全実行し単元条件・核・rank・全合成・多重次数・係数対称性を再確認しました。前回の `verify_ou_five_zero_triangle.py` の貼り合わせ、核、rank、合成、多重次数は既存の検証記録です。統合 `verify_ou_five_zero_supports.py` を4区間に分け全71排除証明書を厳密に再実行しました。`verify_ou_five_zero_curves.py` の全実行で90格子idealと形式方程式を検証しました。前回の最後の六零残余の貼り合わせ・核・合成rankと多重次数の `verify_ou_six_zero_residual.py` は既存の検証記録です。以前の全201台の分割経路rank判定と統合 `verify_ou_six_zero_supports.py` の先頭20台の試行、`verify_ou_six_zero_curves.py` の全実行（平方条件、45格子点、形式表現と核・位数）は既存の検証記録です。前回の `verify_ou_seven_zero_six_cycle.py` のJacobian小行列式とhead条件は既存の検証記録です。前回の `verify_ou_seven_zero_last.py` の形式曲線と核・位数、125531767の八次商、特殊値の八次九次全端点の計算は既存の検証記録です。前回の `verify_ou_seven_zero_line_conics.py` の三族の記号恒等式・特殊値の六次経路商は既存の検証記録です。全族排除は既存の多重次数必要条件に依存します。前回の `verify_ou_zero_supports.py` の7–12零全範囲の軌道数・正規化・厳密経路rank・周期条件と `verify_ou_seven_zero_curves.py` の二曲線の平面切断・特異パラメータ・Hessianは既存の検証記録です。前回の `verify_ou_residual_three_lines.py` の三直線分解、貼り合わせ、核、多重次数、二つの係数例の検算は既存の記録です。前回の `verify_ou_coordinate_dependence.py` の小行列式恒等式と七次rank716は既存の検証記録です。以前の `verify_ou_grid_higher_syzygies.py` の8立方体型・四次公式・六次七次rank231/716・五次六次rank73/231、全256台と核rank障害の二例の検算は既存の記録です。標識付き認識には標準noetherian qgr/導来比較を仮定します。既存の検証では `python3 computations/verify_ou_determinant_cyclic.py` を実行し、Euler行列、偶奇置換、双対Gauss恒等式、古典九直線の巡回類、巡回復元の例を検証しました。可換認識・可約曲線の復元・今回の高次syzygy排除は主論文の証明によります。周期対象の検算は既存の検証記録として保持します。節点型の全184候補の自由群検算は既存の検証記録として保持します。他の計算の過去の実施記録は今回の全計算再実行とは区別します。今回 `verify_ou_five_zero_boundary.py` を実行し、117374702の七次関係行列の有理関数体rank716を確認しました。当時の残余二族の十次までの有限体検算は観察であり、今回の排除とは別です。今回 `verify_ou_five_zero_complement.py` を全実行し、両五零残余の係数・三面行列式・形式表現と核・両向きのhead条件を厳密に検算しました。周期性・単純性と圏論的排除は主稿の証明によります。今回 `verify_ou_four_zero_boundary.py` を全実行し、九境界の三面行列式、三つの有理関数体経路商、円錐曲線のHessianと交叉上の単元・貼り合わせ・核・全合成を厳密に検算しました。一般形式枝と圏論的排除は原稿の証明によります。今回の全runner・新規証明の独立監査・新規性確認・組版は未実施です。

## 二次包絡と商スタックの精密化
標準三次Hilbert級数と一次生成性の下で、二次包絡Q(R)が14641型AS正則環であることと、Rがその正則正規三次商になることは同値です（Levasseurの整域性定理に依存）。トーリック商スタック[P²/μ₃]では反標準AZ-ampleなtorsion-free対象はrank 3以上であり、正則表現が最小rankを達成します。詳細な証明は主論文のIntroductionにあります。

## 公表用主論文の構成

量子Fermat点、非Calabi--Yau符号型二次元族、`\mathbb F_1`上の次数3 Clifford order、次数2・1のClifford orderを同一原稿で扱います。数学的に別の二種類の障害を区別します。

| ケース | 証明する性質 | 関係する制限 |
|---|---|---|
| 次数3・量子Fermat | HH_{-1,0,1}=(2,5,2)、OUのAS surfaceでない | AS四次元ambientの三次超曲面**ではある** |
| 次数3・非CY符号型二次元族 | HH=(2,4,2)、OU非AS性 | AS四次元ambientの非中心正規三次商**ではある** |
| 次数3・Clifford | 反標準環は次数1生成、二次包絡の三次部分が19次元、構造対象を忘れた非表示 | connected gradedな標準quadratic AS4 ambientの正規三次商について |
| 次数2・Clifford | (1,1,1,2)型ambientの正規四次商とのunpointed非同値、HH=(3,4,3) | 別の(1,1,1)型AS4 ambientの中心四次商にはなる |
| 次数1・Clifford | (1,1,2,3)型ambientの正規六次商とのunpointed非同値、HH=(5,3,5) | 別の(1,1,3)型AS4 ambientの中心六次商にはなる |

次数`e`は `dim R_n=1+e n(n+1)/2` で規格化した反標準Hilbert次数です。Clifford orderでの有理標準因子の自己交点は`e/4`です。**奇数次HHの非零性だけでは、AS四次元ambient内の超曲面表示の非存在を証明できません。** 本文では両者を別々に扱います。

主論文から `papers/sections/cubic_*.tex`、`weighted_*.tex` を読み込みます。これらは別の論文ではなく、同じ原稿の証明節です。著者欄は空欄、日本語LuaLaTeXです。

## 再生成

TeX Live/MacTeX（LuaLaTeX、LuaTeX-ja、BibTeX、latexmk）、Python 3、SymPy 1.14.0、Singular 4.4を利用します。

```sh
python3 -m pip install -r requirements.txt
python3 scripts/verify_all.py
latexmk papers/nc-cubic-non-as.tex
# または
python3 scripts/build.py --all
python3 scripts/make_bundle.py
```

`scripts/verify_all.py` は六標識の整数行空間・可約三次曲線の多重次数の計算、元の有限群sector計算と、次数1・2・3のClifford次数・判別曲線・二次包絡の検証をまとめて実行します。証明書 `certificates/degree3_quadratic.json` は45次小行列式と核ベクトルを保持します。有限計算は圏論的なMorita同値や一般の分類定理を形式証明するものではありません。

PDFとSyncTeX、auxiliary filesは`build/`に生成し、Gitでは追跡しません。主論文のPDFは`build/nc-cubic-non-as.pdf`、研究ノートは`build/notes.pdf`です。`python3 scripts/make_bundle.py` は必要なソース・証明書・再現コード・二つのPDFを含む現行bundleを`output/bundles/`に作り、CRCを確認します。PDFなしなら`--source-only`を指定します。

### VS Code / LaTeX Workshop

リポジトリのルートをVS Codeで開き、主論文`papers/nc-cubic-non-as.tex`を表示してください。`Cmd+Option+B`でビルド、`Cmd+Option+V`でPDF表示、`Cmd+Option+J`でSyncTeXによる位置対応です。節ファイルにも主原稿へのroot指定があります。`.vscode/settings.json`のレシピと`.latexmkrc`はいずれも`build/`を出力先とします。変更をpullした後は必要なら **Developer: Reload Window** を実行してください。

## 文献と検証上の注意

書誌情報は [references.bib](references.bib) に統合してあります。Clifford/conicとramification型は既知の枠組みです。量子Fermat・符号型のH3認識の一部は [nc-cubic-classification](https://github.com/uedakazushi/nc-cubic-classification) の結果に依存します。初出・新規性の網羅的調査、非周期的`\mathbb Z`-algebra ambientへの拡張、相対変形の成分性は未解決です。追加仮定と証明の所在は研究方針・研究ノートで確認してください。

少数零型について補集合八係数の正規化可能性を原始二項式で特徴付け、正規化後の三面可逆性が残余スカラー作用に不変であることを証明。零集合222、211/222、200/211/222、121/211/222、120/211/222の五台でこの二項式超曲面の非空開集合を標準比較付きI_3外とした。台全域の判定ではない。
