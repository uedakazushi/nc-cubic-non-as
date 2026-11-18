#!/usr/bin/env python3
"""Create a current, self-contained research bundle without old versions or logs."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import argparse

ROOT=Path(__file__).resolve().parents[1]
FILE_NAMES=[
    ".latexmkrc",".gitignore","README.md","AGENTS.md",
    "00_research_plan.md","requirements.txt","references.bib",
    ".vscode/settings.json",".vscode/tasks.json",
]
SOURCE_DIRS={"papers":{".tex"},"research":{".tex",".md"},
             "tex":{".tex"},"computations":{".py",".json"},
             "scripts":{".py"},"certificates":{".json"}}
EXCLUDE={".build","build","output","__pycache__",".git",".venv","archive","archives"}
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source-only",action="store_true")
    a=parser.parse_args()
    files=[ROOT/name for name in FILE_NAMES]
    for d,exts in SOURCE_DIRS.items():
        directory=ROOT/d
        if directory.exists():
            files.extend(p for p in directory.rglob("*") if p.is_file()
                         and p.suffix in exts
                         and not EXCLUDE.intersection(p.relative_to(ROOT).parts)
                         and not any(x.startswith(".") for x in p.relative_to(ROOT).parts))
    if not a.source_only:
        for stem in ("nc-cubic-non-as","notes"):
            files.append(ROOT/"build"/(stem+".pdf"))
    for f in files:
        if not f.is_file(): raise SystemExit("Missing: "+str(f))
    dest=ROOT/"output"/"bundles"/("nc-cubic-non-as-2026-11-19-v5"+("-source" if a.source_only else "")+".zip")
    dest.parent.mkdir(parents=True,exist_ok=True)
    with ZipFile(dest,"w",ZIP_DEFLATED) as z:
        for f in sorted(set(files)):
            z.write(f,"nc-cubic-non-as/"+f.relative_to(ROOT).as_posix())
    with ZipFile(dest) as z:
        assert z.testzip() is None
    print("Validated bundle:",dest)
if __name__=="__main__":main()
