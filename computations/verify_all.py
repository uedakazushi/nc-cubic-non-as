#!/usr/bin/env python3
"""Reproduce finite quantum Fermat/O-U obstruction certificates (not a proof of derived equivalence)."""
import json,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
for name in ("verify_quantum_fermat.py","verify_quantum_fermat_ou_sectors.py"):
    subprocess.run([sys.executable,str(root/name)],cwd=root,check=True)
cert=json.loads((root/"quantum_fermat_ou_sector_certificate.json").read_text())
counts=cert["hh_dimensions"]
assert counts["-1"]==2 and counts["1"]==2 and counts["0"]==5
print("PASS finite algebra data; twisted sectors HH_-1=HH_1=2, HH_0=5")
