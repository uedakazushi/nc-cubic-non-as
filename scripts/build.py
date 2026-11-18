#!/usr/bin/env python3
"""Build the single publication paper (and optionally the current research notes)."""
from pathlib import Path
import argparse
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
PAPER = "papers/nc-cubic-non-as.tex"
NOTES = "research/notes.tex"

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--all", action="store_true", help="also compile research notes")
    a = p.parse_args()
    if not shutil.which("latexmk") or not shutil.which("lualatex") or not shutil.which("bibtex"):
        raise SystemExit("TeX Live/MacTeX with latexmk, lualatex and bibtex is required.")
    for source in ((PAPER, NOTES) if a.all else (PAPER,)):
        cmd = ["latexmk","-r",str(ROOT / ".latexmkrc"),"-lualatex",
               "-synctex=1","-interaction=nonstopmode","-halt-on-error",
               "-file-line-error","-outdir=build",source]
        print(" ".join(cmd),flush=True)
        subprocess.run(cmd,cwd=ROOT,check=True)
        stem = Path(source).stem
        for ext in (".pdf",".synctex.gz"):
            f = ROOT / "build" / (stem + ext)
            if not f.is_file():
                raise SystemExit("Missing expected output: " + str(f))
        print("Built",ROOT / "build" / (stem + ".pdf"))

if __name__ == "__main__":
    main()
