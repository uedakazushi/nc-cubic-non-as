# 非可換del Pezzo曲面の非AS性と超曲面表示の障害

**現行版 2026-11-19 / non-as-v5。** `shinnosukeokawa/nonas-ncdp` の次数1・2・3のClifford成果を統合し、主論文を **[papers/nc-cubic-non-as.tex](papers/nc-cubic-non-as.tex) 一本** にまとめています。

研究を再開する際は、最初に **[研究方針](00_research_plan.md)**、次に **[研究ノート](research/notes.tex)** を読み、必要な箇所を主論文の節ファイルで確認してください。旧リポジトリや旧版のbundleは必要ありません。

## 内在的定義の三段階と境界例

次数3の内在的候補 `I_3`、標準AS四次元三次超曲面型 `H_hyp`、OUのAS surfaceでもある超曲面型 `H_hyp∩AS_OU` を区別し、真の包含 `H_hyp∩AS_OU ⊊ H_hyp ⊊ I_3` をIntroductionで証明しています（rank一の符号型族と次数三Clifford例を証人とします）。**OUのAS surface全体が標準AS4超曲面になるとは主張しません**。rank一・Hilbert級数を緩めた場合の境界例として `[P²/μ₃]` を計算し、rank三のAZ-ample対象と、構造層のnon-ampleness、特異な反標準三次模型 `UVW=T³` を区別しました。

## 新たに証明した精密化
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

`scripts/verify_all.py` は元の有限群sector計算と、次数1・2・3のClifford次数・判別曲線・二次包絡の検証をまとめて実行します。証明書 `certificates/degree3_quadratic.json` は45次小行列式と核ベクトルを保持します。有限計算は圏論的なMorita同値や一般の分類定理を形式証明するものではありません。

PDFとSyncTeX、auxiliary filesは`build/`に生成し、Gitでは追跡しません。主論文のPDFは`build/nc-cubic-non-as.pdf`、研究ノートは`build/notes.pdf`です。`python3 scripts/make_bundle.py` は必要なソース・証明書・再現コード・二つのPDFを含む現行bundleを`output/bundles/`に作り、CRCを確認します。PDFなしなら`--source-only`を指定します。

### VS Code / LaTeX Workshop

リポジトリのルートをVS Codeで開き、主論文`papers/nc-cubic-non-as.tex`を表示してください。`Cmd+Option+B`でビルド、`Cmd+Option+V`でPDF表示、`Cmd+Option+J`でSyncTeXによる位置対応です。節ファイルにも主原稿へのroot指定があります。`.vscode/settings.json`のレシピと`.latexmkrc`はいずれも`build/`を出力先とします。変更をpullした後は必要なら **Developer: Reload Window** を実行してください。

## 文献と検証上の注意

書誌情報は [references.bib](references.bib) に統合してあります。Clifford/conicとramification型は既知の枠組みです。量子Fermat・符号型のH3認識の一部は [nc-cubic-classification](https://github.com/uedakazushi/nc-cubic-classification) の結果に依存します。初出・新規性の網羅的調査、非周期的`\mathbb Z`-algebra ambientへの拡張、相対変形の成分性は未解決です。追加仮定と証明の所在は研究方針・研究ノートで確認してください。
