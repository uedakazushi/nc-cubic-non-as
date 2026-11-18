# 研究方針：非AS型の内在的非可換三次曲面
2026-11-17 / non-as-v2 / **研究再開時に最初に読む**。次に主論文 papers/nc-cubic-non-as.tex と research/research_notes.md を読む。正本はGitHub main。

## I. 恒常的な研究原則・保存許可
数学的重要性、新規性、理解の深化、応用可能性、正確性、実現可能性を総合して価値を最大化する。必要に応じ方針を現行版として書き直し、変更履歴を蓄積しない。主張を証明済み・仮定付き・未完成・予想・計算的観察に区分し、原稿の証明／今回の監査／計算の実行を混同しない。先行研究の一次資料を確認し、未確認事項や新規性は未確認と記す。論文は日本語LuaLaTeXを基本とし、十分な証明を残す。ノートは現行有効な着想・障害を整理し、作業履歴・旧ZIP・重複・キャッシュは収録しない。

研究区切りの承諾不要Library保存は一般許可されるが、この研究ではGitHub mainを正本とし、push成功時はLibrary保存せず、push失敗時のみ緊急避難に利用する。既存Libraryファイルは明示指示なしに削除しない。PDF・ZIP・中間生成物はGit追跡せずTeXソースから再生する。PR、Actions、force pushを使わずmainへfast-forwardする。生成・検証・push・Library保存を別状態として報告する。

## II. 目標と到達点
Okawa--Ueda Definition 1.1 のAS surfaceと、内在的H3条件の差を明確な導来Morita不変量で公表する。主論文は papers/nc-cubic-non-as.tex。
- 中心量子Fermat（三乗根、CY ambient、中心三次）：捩れ付き(mu3)^2商、HH_{-1,0,1}=(2,5,2)。
- 符号型（二位数、非CY ambient、非中心正規三次）：G=f3(x0,x1)+x0*x2^2+x1*x3^2、ad*disc(f3)非零。捩れ付きH0=(Z/2)^2商、HH=(2,4,2)。一般pointed同型類次元2。
- OUのAS surfaceなら完全強例外列があり、HHは次数0に集中。従ってどちらもOU型でない。

符号型のH3認識およびfixed ambientでの二次元分類は nc-cubic-classification/papers/sections/cubic_noncy_sign_example.tex に依拠する既存結果。非AS性とHH計算は本原稿の新しい部分。全intrinsic moduliの成分性は未証明。

## III. 優先課題
1. 符号型の退化twisted group algebra T=M2⊕M2に対するH作用・Reynolds full性・H0 corner縮約とprojective表現を独立に再監査する。
2. 二本の楕円曲線のPoincare residueの指標、Hodge次数とdiscrete torsionの符号、HH0固定点寄与を点検する。
3. OU Definition 1.1、Arinkin--Caldararu--Hablicsek、Baranovsky、Caldararu--Giaquinto--Witherspoonの一次資料と本稿の適用条件を精査し、新規性を確認する。
4. 符号型二次元族のversal変形と全H3モジュライの既約成分性、OU locusの開閉性を調べる。

## IV. 再現
- python3 -m pip install -r requirements.txt
- python3 computations/verify_all.py
- latexmk -lualatex papers/nc-cubic-non-as.tex

finite-sectorプログラムは群と指標の有限計算のみを検査し、dg Morita同値やorbifold HKRの数学的証明を代替しない。
