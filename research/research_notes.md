# 非AS型の内在的三次曲面：現行研究ノート
2026-11-17 / non-as-v2。研究方針は 00_research_plan.md、本文証明は papers/nc-cubic-non-as.tex。

## 二種類の成果と既存入力
量子Fermat: CY ambientと中心三次。L=(Z/3)^2、H=L^vee、T=C_mu L=M3。Hの楕円固定sector4つがHH_{±1}に計2ずつ寄与し、HH0=5。smooth central切断のpointed同型類は一点。これまでの証明を主論文に保持する。

符号型: q01=1、他の5組は-1。G=f3(x0,x1)+x0*x2^2+x1*x3^2、ad*Delta(f3)≠0。非CY ambientと非中心三次。H3認識とfixed sign ambientのpointed二次元分類は nc-cubic-moduli の papers/sections/cubic_noncy_sign_example.tex に既存証明がある。今回の新しいHochschild計算は HH_{-1,0,1}=(2,4,2)。

## 符号型の本質的な注意点
Gamma=(Z/2)^3、重み(x0,x1,x2,x3)=(e0,e0,e1,e2)、cocycle nu(u,v)=(-1)^(u1*v0+u2*v0+u2*v1)。交換形式のradicalは e0+e1+e2。従ってT=C_nu Gamma=M2⊕M2であり、M4と誤記してはいけない。H=Gamma^veeは二つの行列成分を交換し、その安定化群H0={(a,b,c):a+b+c=0}へMorita還元する。非退化Schur multiplierはH0上の一意の非自明類で、commutatorは(-1)^(a_g b_h+b_g a_h)。

H0非自明元gA=(0,1,1)の固定集合はP1+3点(Euler5)。gB=(1,0,1)はE2={x2=0}と一点(Euler1)。gC=(1,1,0)はE3={x3=0}と一点(Euler1)。固定楕円のH0作用はH^0(Omega1)上それぞれ(-1)^b、(-1)^cでありtwisted characterに一致する。Picard H^2のtraceは7,3,-1,-1、invariant rankは2。従ってHH0=4、HH±1=2。Pythonのfinite sector確認は computations/verify_sign_family_sectors.py。

## 注意・未解決
実際のMorita同値とHKRは本文の代数論・引用文献に依存し、Python有限列挙だけでは証明されない。符号型のwhole family H3認識は既存分類の精密な仮定に依存する。全intrinsicモジュライの成分性・相対降下・versal変形・OU locusの開閉性・新規性の網羅調査は未解決。量子Fermatと符号型のHH0が異なるのはH0の不変Picard rankとゼロ次sectorの違いによる。

再現: python3 computations/verify_all.py; latexmk -lualatex papers/nc-cubic-non-as.tex。
