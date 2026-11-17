# 研究方針：量子Fermat三次曲面の非AS性

2026-11-17 / ou-obstruction-v1 / 独立研究正本。**再開時に最初に読む。** 次に `research/research_notes.md` と `papers/quantum_fermat_not_as_surface.tex`、必要な検算へ進む。

## I. 恒常的な研究原則と保存許可

数学的価値・正確性・理解の深化・新規性・公表可能性を総合して研究方針を変更する。重要な変更は現行方針として理由と残る課題を記し、履歴の追記にしない。証明済み・仮定付き・未完成・予想・計算的観察を区別し、他論文の既存成果と本稿の成果を分ける。一次文献を確認し、未確認の新規性や書誌を補わない。原稿は原則日本語LuaLaTeX、十分な証明を保持する。ノートは現在有効な知識と未解決障害を整理し、過去の会話・旧ZIP・ログ・キャッシュ・PDFをGitに含めない。

研究区切りの承諾不要Library保存は一般許可されるが、この研究ではGitHub mainを正本として、push成功時はLibraryに保存せず、GitHub保存に失敗した場合だけ緊急避難に使用する。既存Libraryファイルは明示指示なしに削除しない。

## II. 目標と達成基準

量子Fermat三次曲面の非可換射影幾何と、Okawa--UedaのDefinition 1.1のAS surfaceとの差を、独立した短報として完全に証明し公表する。導来Morita同値、離散torsionのtwisted orbifold HKR、固定楕円曲線の指標とHH次数、OUの例外列の帰結を全て正確に示す。反例が持つ内在的H3条件は既存の結果と区別する。

## III. 到達点

短報の中心計算は`dim HH_{-1}=dim HH_1=2`、`dim HH_0=5`。OUのAS surfaceには完全強例外列があり、そのHochschild homologyは次数0に集中するため反例となる。原稿には捩れ付き群商stackへの導来Morita比較とHodge固有指標の証明を記載。計算コードは`computations/verify_quantum_fermat.py`と`computations/verify_quantum_fermat_ou_sectors.py`、H3例の別検算として`verify_boundary_fermat.py`を保持する。有限計算だけで導来同値や圏の正則性を証明したとは主張しない。

## IV. 最優先課題

(1) 捩れのcocycleとequivariant化・Morita同値の全仮定を独立検査する。(2) twisted HKRのsector別符号と固定曲線の指標を一次資料と突き合わせる。(3) OU Definition 1.1の適用範囲と完全強例外列の帰結を精査し、既存の非AS反例との差をIntroductionで明確にする。(4) 論文の説明を簡潔に整理し、主分類論文への依存を最小限にする。

## V. 再現

```sh
python3 -m pip install -r requirements.txt
python3 computations/verify_all.py
latexmk -lualatex papers/quantum_fermat_not_as_surface.tex
```

生成・検算・commit・pushは別状態として報告し、独立査読や網羅的な文献調査を未実施なら明示する。
