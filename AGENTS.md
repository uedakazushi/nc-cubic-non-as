# Research workflow: nc-cubic-non-as

## 現行範囲
最初に 00_research_plan.md を読む。AOU三ブロック箙型の全支持・全特殊化に対する包含は主稿 thm:ou-aou-all-unconditional に証明がある。一般OU型は未解決。現行の限定範囲を無断で変更しない。次は独立監査と公表準備を優先する。
主原稿は papers/nc-cubic-non-as.tex、分割節は同一論文の一部。研究ノートは research/notes.tex。量子Fermat、符号型、次数1・2・3 Clifford、OU認識の有効な証明を保持する。日本語LuaLaTeX、著者空欄。
認識・排除・実現可能性を区別し、標準比較、内在的条件、very strong性の仮定を落とさない。有限体計算は探索用、排除には厳密証明を要求する。特異曲線のPerfとDbを混同しない。補集合の外側面行列式は補集合cubeの正規化後に判定する。周期対象・AZ1/AZ2・Euler rank・noetherian比較の依存関係を独立監査する。
再現手順は研究方針とREADMEを参照。全runner・組版・独立監査・新規性確認を区別する。一次文献の定理と仮定を明記する。

Maintain a self-contained lightweight current set: plan, notes, paper, bibliography, necessary reproducibility code and certificates. No logs, old ZIPs, backups or redundant drafts. No routine Library duplicate or ZIP for GitHub-managed research; use Library only for emergency preservation if push fails or when expressly requested. Never delete Library files without explicit instruction.

GitHub main is source of truth. Push fast-forward to main after a coherent piece of work (target 15–30 minutes); do not use PRs, Actions or force-push. Do not commit generated PDF, SyncTeX, auxiliary files or ZIP. Maintain `.latexmkrc` and `.vscode/settings.json` so both `latexmk` from the repository root and LaTeX Workshop use `build/`. If a repository state changes remotely, reread the head before updating.

After a research milestone, first update and save the current paper, notes, research plan and README in this repository. Then read the latest `uedakazushi/research-hub` README and AGENTS.md and update the relevant concise entry (actual achievements, certainty, unresolved core, next step, dependencies, canonical paths, actual inspection date), preserving unrelated/concurrent changes. Use normal main commits without force push or Actions, and verify repository and hub updates separately. Report any failed or unverified step. Do not duplicate proofs or detailed todos in the hub.

Original source projects: `uedakazushi/nc-cubic-non-as` and `shinnosukeokawa/nonas-ncdp`; this repository now contains their current mathematical results, so do not require retrieving either source for resumption. Do not delete or mutate the source repository without explicit instructions.
