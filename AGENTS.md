# Codex研究指示 — nc-cubic-non-as

2026-11-17。まず`00_research_plan.md`、次に`research/research_notes.md`と`papers/nc-cubic-non-as.tex`を読む。

このリポジトリは量子Fermat点と非CY符号型二次元族のOU AS surface障害の論文の独立正本。全分類、ambient比較、楕円E6の各研究は別リポジトリのmainを正本とし、ここから誤って編集しない。

数学的価値・正確性・理解の深化・新規性・公表可能性を総合して研究方針を変更する。重要な変更は現行方針として理由と残る課題を記し、履歴の追記にしない。証明済み・仮定付き・未完成・予想・計算的観察を区別し、他論文の既存成果と本稿の成果を分ける。一次文献を確認し、未確認の新規性や書誌を補わない。原稿は原則日本語LuaLaTeX、十分な証明を保持する。ノートは現在有効な知識と未解決障害を整理し、過去の会話・旧ZIP・ログ・キャッシュ・PDFをGitに含めない。

研究区切りの承諾不要Library保存は一般許可されるが、この研究ではGitHub mainを正本として、push成功時はLibraryに保存せず、GitHub保存に失敗した場合だけ緊急避難に使用する。既存Libraryファイルは明示指示なしに削除しない。

約15分ごと・論理的な区切りで、証明・方針・検算コードの現行版を整理し、PR、Actions、force pushを使わずmainへ通常のfast-forwardでpushする。共有mainが進んだらfetchして両立させる。PDFは追跡せず、VSCode/latexmkにより再生成する。主要検算は`python3 computations/verify_all.py`。新規性や全圏論的同値の独立監査を計算PASSだけから主張しない。
