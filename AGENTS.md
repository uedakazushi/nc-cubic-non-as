# Research workflow: nc-cubic-non-as

Read `00_research_plan.md` first. The only standalone publication paper is `papers/nc-cubic-non-as.tex`; its `papers/sections/` inputs are not separate papers. Preserve every currently valid proof (quantum Fermat, sign family, degree 3/2/1 Clifford orders), and integrate new results into the logical position, not historical addenda. Research notes are `research/notes.tex`.

Use Japanese LuaLaTeX. Keep the author blank unless explicitly instructed otherwise. Separate proved theorems, assumptions, gaps, computations, and literature claims. Cite primary sources precisely and label unverified novelty. Run `python3 scripts/verify_all.py` and `python3 scripts/build.py --all` when the environment provides dependencies. Do not claim those checks passed unless they actually ran.

One current self-contained lightweight bundle only: plan, notes, paper, bibliography, necessary reproducibility code, certificates and available compiled PDF(s). No logs, old ZIPs, backups or redundant drafts. May save the updated bundle to Library at research milestones without additional permission; never delete Library files without explicit instruction.

GitHub main is source of truth. Push fast-forward to main after a coherent piece of work; do not use PRs, Actions or force-push. Do not commit generated PDF, SyncTeX, auxiliary files or ZIP. Maintain `.latexmkrc` and `.vscode/settings.json` so both `latexmk` from the repository root and LaTeX Workshop use `build/`. If a repository state changes remotely, reread the head before updating.

Original source projects: `uedakazushi/nc-cubic-non-as` and `shinnosukeokawa/nonas-ncdp`; this repository now contains their current mathematical results, so do not require retrieving either source for resumption. Do not delete or mutate the source repository without explicit instructions.
