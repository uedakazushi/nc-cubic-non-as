# 非AS型非可換三次曲面 — nc-cubic-non-as

**まず [研究方針](00_research_plan.md) を読む。** 公表用LuaLaTeX原稿は [nc-cubic-non-as.tex](papers/nc-cubic-non-as.tex)。補足は [研究ノート](research/research_notes.md)。

Okawa–UedaのAS surfaceより内在的な三次del Pezzo型曲面の条件が広いことを、二種類の例で示す。

- **量子Fermat点**：CY ambient、中心三次、Hochschild homologyの次元（−1,0,1）=(2,5,2)。
- **符号型二次元族**：非CY ambient、非中心正規三次。滑らかさ ad*Delta(f3)≠0 のもとでHH=(2,4,2)。固定符号ambientの一般pointed同型類の次元は2。

どちらもdiscrete torsion付き有限商stackへの導来Morita同値と、固定楕円曲線上のtwisted orbifold HKRによる。符号型のH3認識とpointed分類は[ H3主研究](https://github.com/uedakazushi/nc-cubic-moduli)の既存結果に依存する。全intrinsicモジュライの成分性・新規性の網羅的調査は未完成。

再現:
- python3 -m pip install -r requirements.txt
- python3 computations/verify_all.py
- latexmk -lualatex papers/nc-cubic-non-as.tex

GitHub mainを正本とし、PDF・ZIP・中間生成物は追跡しない。関連研究は [ambient比較](https://github.com/uedakazushi/nc-cubic-ambient-moduli) および [H3分類](https://github.com/uedakazushi/nc-cubic-moduli)。
