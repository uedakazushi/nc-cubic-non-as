#!/usr/bin/env python3
"""Reproduce finite quantum Fermat/O-U obstruction certificates (not a proof of derived equivalence)."""
import json,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
for name in ("verify_quantum_fermat.py","verify_quantum_fermat_ou_sectors.py","verify_sign_family_sectors.py","verify_toric_mu3.py"):
    subprocess.run([sys.executable,str(root/name)],cwd=root,check=True)
cert=json.loads((root/"quantum_fermat_ou_sector_certificate.json").read_text())
counts=cert["hh_dimensions"]
assert counts["-1"]==2 and counts["1"]==2 and counts["0"]==5
sign=json.loads((root/"sign_family_sector_certificate.json").read_text())
assert sign["hh_dimensions"]=={"-1":2,"0":4,"1":2}
print("PASS finite sector checks: Fermat HH=(2,5,2); sign family HH=(2,4,2)")
