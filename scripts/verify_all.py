#!/usr/bin/env python3
"""Run all independently reproducible finite checks from both source programs."""
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
PROGRAMS=(
    "computations/verify_all.py",
    "computations/verify_ou_decuple_markings.py",
    "computations/verify_ou_reduced_multidegrees.py",
    "computations/verify_ou_nodal_mapping_classes.py",
    "computations/verify_ou_reduced_algebraic_moves.py",
    "computations/verify_ou_iv_relation_matches.py",
    "computations/verify_ou_iv_periodic_objects.py",
    "computations/verify_ou_plane_curve_reduction.py",
    "computations/verify_ou_determinant_cyclic.py",
    "computations/verify_ou_dense_higher_syzygies.py",
    "computations/verify_ou_grid_higher_syzygies.py",
    "computations/verify_ou_coordinate_dependence.py",
    "computations/verify_ou_residual_three_lines.py",
    "computations/verify_boundary_fermat.py",
    "scripts/del_pezzo_degree2_verification.py",
    "scripts/nc_dp1_quintic_checks.py",
    "scripts/hochschild_verification.py",
    "scripts/verify_degree3_smooth.py",
    "scripts/verify_degree3_relations.py",
)
def main():
    for file in PROGRAMS:
        print("VERIFY",file,flush=True)
        subprocess.run([sys.executable,str(ROOT / file)],cwd=ROOT,check=True)
    print("All requested finite calculation scripts completed.")
if __name__=="__main__": main()
