# 非可換del Pezzo曲面の非AS性と超曲面表示の障害

**現行版 2026-11-24 / non-as-v31。** `shinnosukeokawa/nonas-ncdp` の次数1・2・3のClifford成果を統合し、主論文を **[papers/nc-cubic-non-as.tex](papers/nc-cubic-non-as.tex) 一本** にまとめています。

研究を再開する際は、最初に **[研究方針](00_research_plan.md)**、次に **[研究ノート](research/notes.tex)** を読み、必要な箇所を主論文の節ファイルで確認してください。旧リポジトリや旧版のbundleは必要ありません。

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
- **六個以上の零係数型の認識・排除の完了（標準比較付き）**（`thm:ou-six-zero-path-exclusion`, `thm:ou-six-zero-periodic-exclusion`, `lem:ou-two-zero-slice-grid`, `thm:ou-six-zero-curve-recognition`, `thm:ou-six-zero-expanded-curve`, `cor:ou-six-zero-complete-recognition`）：201台中53/76/30/23/12台を長さ4/5/6/7/8の厳密経路商rankで排除（計194、有限体は候補端点のscreeningのみ）。残る七台は117358318,117370606,125759214,131005679,131037935,133099246,134145774（131037935だけ三係数、他は二係数）。131037935は既存周期障害、125759214は形式曲線 R=-v+t^3, D'=1−(v−1)R/v, h=(u+D')(R+D')/(RD')（u+v≠0、u=−v≠−1、(-1,1)の別表現を全て扱う）で標準比較付きI_3外。五族の二項三次切断をxの二次式として原始性・判別式平方条件・45偶奇対の勾配idealで厳密判定。117370606,133099246,131005679は全非零係数で既約、134145774は(-1,-1)だけ直線+円錐曲線だが各sliceの二零が異なる行列位置にある場合の新しい九線束復元補題で全域認識。117358318のu≠1は既約切断で認識、u=1もN=V(XZ−Y²)とL=V(X−Y+vZ)の長さ二の交叉D=s²−st+vt²での明示的な線束貼り合わせ・核・合成比較により全v≠0で認識（v=1/4の接触型を含む）。δ_md=0,w=(3,0)なので全平面三次曲線認識定理を適用。六零の未判定は零。残る1–5零と標準比較。認識は標準比較・I_3付き、実現可能性は未主張。
- **七個以上の零係数型の認識・排除（標準比較付き）**（`prop:ou-zero-support-enumeration`, `thm:ou-large-zero-path-exclusion`, `thm:ou-seven-eight-zero-path-exclusion`, `thm:ou-latin-square-periodic-exclusion`, `cor:ou-intrinsic-at-most-seven-zeros`, `prop:ou-seven-zero-curve-recognition`）：9 sliceの連結性と異slice両順序の強制零を含まない木成分の排除で全台を列挙（零係数数0–12の同値類数1,1,3,9,30,80,201,344,387,217,53,4,1）。torus正規化・厳密経路rankにより10–12零の58台全部、9零の217台中216台、8零の387台中386台、7零の344台中336台を排除。8零の唯一の残余108979950は9零のLatin square型の零を一つ除いた台。対角三零・一方向巡回三零・ちょうど二添字が等しい18係数非零に周期単純対象と任意長有限長商による障害を一般化し、9零、8零の残余と7零の125234926をI_3から排除（標準比較付き）。従ってI_3内のvery strong OU型の零係数は高々7個。七零の116309742は自由係数u=c_211≠1で滑らかな切断から認識、116834030はu=c_222≠0の全域で認識（u=1/16の唯一の通常節点は六座標非零で核rank2）。定理 `thm:ou-seven-zero-line-conic-exclusion` は128957047,129972855,133118583をvery strong性から全係数で排除。128957047のu=1は六次商次元11（必要値10）、他は直線+滑らかな円錐曲線の切断と明示的な九線束の貼り合わせ、核、直交性・合成比較により δ_md=(-1,1)≠0 を得て既存の多重次数必要条件に反する（標準比較・I_3不要）。さらに `thm:ou-seven-zero-formal-arc-exclusion` で117368558を形式表現と有界head（u≠−1で位数1、u=−1で3）による固定対象の任意長有限長商・周期rank障害から全係数でI_3外とした（標準比較付き）。`thm:ou-seven-zero-eighth-path-exclusion` は125531767を八次経路商の有理関数体rank2170、商次元17>15から全係数でvery strong性外とした。最後の116309742のu=1も `thm:ou-seven-zero-six-cycle-exclusion` により標準比較の下でI_3外とした。六孤立零に対応するs固定単純対象の同時厚化を各頂点rank2の形式表現で構成。q=t^6、ラベル差の次数と三周期対称性から18方程式・30変数に帰着、対角式をqで割った初期値でJacobian小行列式−8を厳密確認。形式陰関数定理と位数1の二次minorにより固定始域の任意長有限長商を得て 2r([O])=2≥6 の矛盾。`cor:ou-seven-zero-complete-recognition` により七個以上の零係数型の認識・排除を完了（標準比較付き）。六零以下の現状は上の帰着定理を参照。台の実現可能性は未主張。
- **全係数非零型の認識（標準比較付き）**（`thm:ou-dense-residual-three-lines`, `thm:ou-dense-complete-recognition`）：残余五係数族 F1=0 の零平面はz0=z1=t、x+w=y+z、(D+S)t−Sx+(b−d)p+(c−d)q=0（p=y−x,q=z−x）となり三次式はtpq=0。S≠0は三角形、S=0は共点三直線。三成分上の射行列を特異点を含むxw≠0の共通開集合上の行列と単元の行列rescalingで貼り合わせ、共点の場合も線束を直接構成。核・余核は全点rank2から局所自由、九線束のmdegはU0=(1,0,−1), U1,U2=(0,0,0), V0=(0,1,0), V1,V2=(1,0,0), W0=(1,0,1), W1,W2=(1,1,0)。δ_md=0,w=(1,1,1)、前向きExt消滅・ブロック内完全直交・合成関係を証明。既存のrank2認識・rank3復元・残余への帰着と合わせて全係数非零かつvery strong型を尽くし、標準noetherian OU比較と内在的三次条件の下で標準AS4正規三次商として認識する。新規性・全証明の独立監査は未確認。一般の零係数型・標準比較は未解決。
- **座標従属の五係数族への帰着**（`lem:ou-dense-constant-face-obstruction`, `thm:ou-dense-five-coefficient-residual`）：rank3全係数非零の正規形M=[[a,b],[c,d]], N=λ1+μM（λ≠0）でvery strongならμ≠0, τ=1+λ−μ≠0。μ=0の定数面族の七次1458×729多項式行列（始点・終点添字1）の有理関数体rank716を実行確認、全特殊化で商次元13以上（必要値12）。τ=0は基準slice変更で同じ族。元の巡回が(ii),(iii)を満たすと他の二巡回も(ii),(iii)を満たし、(i)はF1 F2 F3≠0と同値。D=ad−bc, S=a+d−b−c, T=abc−abd−acd+bcd, F1=(μ−1)D+λS, F3=μT−(1+λ)D, F2=F3−F1。全巡回の失敗はkの添字変更でF1=0, D≠0, μ=1−λS/Dに帰着。標準比較・I_3不要。残余族は上の三直線表示により認識済み。
- **核rank障害の高々一巡回性・格子条件との同時成立・残る座標従属**（`thm:ou-dense-at-most-one-torus-failure`, `thm:ou-dense-common-grid-torus`, `prop:ou-dense-dependent-coordinate-normal-form`）：全係数非零・very strongなら(iii)に失敗する巡回は高々一つ。二方向のrank1面から非零辺ベクトルの二部graphを作り、非連結なら四次次元超過、連結ならd_ijk=m_ij z_kで五次双対の次元が8以上（必要値7）。混合障害は一方向の(iii)失敗と他二方向の(ii)失敗に帰着し、障害長方形の比較からdの二面 M=[[0,0],[0,q]], N=u[[1,w],[v,vw]] の族のみ残る。六次405×243多項式行列の有理関数体rank231から全特殊化で商次元12以上（必要値10）となり排除。従って(ii),(iii)の共通成功巡回が存在する。標準比較・I_3不要。rank3でなお(i)が失敗する場合は、従属する座標対の共有辺を(0,0)とし正規化して二面N=λ1+μM（M非定数、λ≠0）に帰着。この六係数族は上の帰着定理によりさらに絞り込まれる。
- **四次経路次元公式と三巡回同時格子障害の排除**（`prop:ou-dense-fourth-path`, `thm:ou-dense-no-triple-grid-failure`）：全係数非零の巡回二次関係代数で、端点(i,j)を含み交叉比が第3添字に依らない長方形数bに対し、長さ4の経路商次元はmax(5,4+b)。very strongなら各巡回の障害長方形は高々1。三巡回同時の(ii)失敗は立方体三辺の係数1への正規化で、二つの四係数族に帰着。食い違う軸が2以上の型は四次公式で排除。残りは有理関数体上の405×243（長さ6）・1458×729（長さ7）行列の厳密rank 231・716から、全特殊化で商次元12以上・13以上となりEulerの必要値10・12に反する。標準比較・I_3不要。`computations/verify_ou_grid_higher_syzygies.py` を実行し8立方体型、四次公式の有理数例と有理関数体rankを確認した。この格子障害排除は上の共通成功巡回定理の入力。
- **高次syzygyによる三巡回同時の核rank障害の排除**（`thm:ou-dense-no-triple-torus-failure`）：全係数非零かつvery strong導来helixなら少なくとも一つの巡回で(iii)が成立する。三方向の障害sliceを0とし係数をrescaleして全座標零面を1にすると、残る2×2×2立方体のc−1の各面は非零rank1。台は全立方体または対頂点のみ。積型では長さ5の経路商次元8（必要値7）、対頂点型では長さ4の経路商次元6以上（必要値5）となる。一般証明は主稿、標準比較・I_3・行列式rank3は不要。`computations/verify_ou_dense_higher_syzygies.py` で全256台と二つの有理数係数例を検算。上の高々一巡回性により強化される。
- **六標識の境界**（`prop:ou-marking-walls`）：第一ブロックの6順序でVdB一般位置を比較し、99の整数等式と3975のrank二行空間を得る。δ非torsionの場合に全標識が失敗する集合は余次元3以上。これはAOUの全仮定を回復する主張ではない。

一般の包含は未解決ですが、全係数非零と六個以上の零係数型は標準noetherian OU比較と内在的三次条件の下で認識・排除済みです。今回、最後の六零残余117358318のu=1を直線・円錐曲線の長さ二の交叉上の明示的な貼り合わせで認識しました。次は五零の80台から1–5零、標準比較を優先します。今回 `verify_ou_six_zero_residual.py` を実行し、交叉上の等式・核・九端点の合成rankと多重次数を検証しました。以前の全201台の分割経路rank判定と統合 `verify_ou_six_zero_supports.py` の先頭20台の試行、`verify_ou_six_zero_curves.py` の全実行（平方条件、45格子点、形式表現と核・位数）は既存の検証記録です。前回の `verify_ou_seven_zero_six_cycle.py` のJacobian小行列式とhead条件は既存の検証記録です。前回の `verify_ou_seven_zero_last.py` の形式曲線と核・位数、125531767の八次商、特殊値の八次九次全端点の計算は既存の検証記録です。前回の `verify_ou_seven_zero_line_conics.py` の三族の記号恒等式・特殊値の六次経路商は既存の検証記録です。全族排除は既存の多重次数必要条件に依存します。前回の `verify_ou_zero_supports.py` の7–12零全範囲の軌道数・正規化・厳密経路rank・周期条件と `verify_ou_seven_zero_curves.py` の二曲線の平面切断・特異パラメータ・Hessianは既存の検証記録です。前回の `verify_ou_residual_three_lines.py` の三直線分解、貼り合わせ、核、多重次数、二つの係数例の検算は既存の記録です。前回の `verify_ou_coordinate_dependence.py` の小行列式恒等式と七次rank716は既存の検証記録です。以前の `verify_ou_grid_higher_syzygies.py` の8立方体型・四次公式・六次七次rank231/716・五次六次rank73/231、全256台と核rank障害の二例の検算は既存の記録です。標識付き認識には標準noetherian qgr/導来比較を仮定します。既存の検証では `python3 computations/verify_ou_determinant_cyclic.py` を実行し、Euler行列、偶奇置換、双対Gauss恒等式、古典九直線の巡回類、巡回復元の例を検証しました。可換認識・可約曲線の復元・今回の高次syzygy排除は主論文の証明によります。周期対象の検算は既存の検証記録として保持します。節点型の全184候補の自由群検算は既存の検証記録として保持します。他の計算の過去の実施記録は今回の全計算再実行とは区別します。新規証明の独立監査・新規性確認・組版は未実施です。

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

TeX Live/MacTeX（LuaLaTeX、LuaTeX-ja、BibTeX、latexmk）、Python 3、SymPy 1.14.0を利用します。

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
