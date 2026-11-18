#!/usr/bin/env python3
"""Verify the saved exact degree-three quadratic-obstruction certificate.

The default run reads the certificate without changing it. Use --regenerate
explicitly to recompute and replace the minor and kernel witnesses. All
arithmetic is exact; neither floating-point nor finite-field ranks are used.
"""
import argparse
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates" / "degree3_quadratic.json"
RELATION_NAMES = [
    "CD-DC", "BC+CB", "BD+DB", "2AA-2CC-AD-DA",
    "AB+BA-2DD", "2BB-2DD-AC-CA",
]


def construct():
    e2 = sp.eye(16)
    w = lambda a, b: e2[:, 4 * a + b]
    # q11 = c s^2 + d t, q12 = c t^2, q22 = c t^2 + d s,
    # q13 = d, q23 = 0, q33 = c.
    relations = [
        w(2, 3) - w(3, 2),
        w(1, 2) + w(2, 1),
        w(1, 3) + w(3, 1),
        2*w(0, 0) - 2*w(2, 2) - w(0, 3) - w(3, 0),
        w(0, 1) + w(1, 0) - 2*w(3, 3),
        2*w(1, 1) - 2*w(3, 3) - w(0, 2) - w(2, 0),
    ]
    e1 = sp.eye(4)
    # First 24 columns r_j x_i; final 24 columns x_i r_j.
    matrix = sp.Matrix.hstack(
        *[sp.kronecker_product(r, e1[:, i])
          for r in relations for i in range(4)],
        *[sp.kronecker_product(e1[:, i], r)
          for r in relations for i in range(4)],
    )
    return sp.Matrix.hstack(*relations), matrix


def check_euler():
    """Recompute the Euler matrix and exceptional-object numerical formulas."""
    def chi_a(a, b):
        return 1 + 2*(-a*a + 2*a*b) + a + b

    def chi_p(a, b):
        return 4 + 2*(-a*a + 2*a*b) + 3*(a+b)

    objects = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    euler = sp.zeros(4)
    for i, (ei, ai, bi) in enumerate(objects):
        for j, (ej, aj, bj) in enumerate(objects):
            e, a, b = ej-ei, aj-ai, bj-bi
            if e == -1:
                e, a, b = 1, a-1, b-2
            euler[i, j] = chi_a(a, b) if e == 0 else chi_p(a, b)
    assert euler.det() == 100
    alpha, beta, y = sp.symbols("alpha beta y")
    coord = sp.Matrix([1-alpha/2-beta/2+y/2, y,
                       (alpha-y)/2, (beta-2*y)/2])
    self_euler = sp.expand((coord.T * euler * coord)[0])
    structure_euler = sp.expand((euler * coord)[0])
    y_exceptional = (-alpha**2 + 2*alpha*beta + 2*alpha)/5
    assert sp.simplify(self_euler.subs(y, y_exceptional)-1) == 0
    assert sp.simplify(structure_euler.subs(y, y_exceptional)
        - (1+(-alpha**2+2*alpha*beta+alpha+beta)/2)) == 0
    return euler


def metadata(euler):
    return {
        "word_order": "lexicographic on A,B,C,D; zero-based indices",
        "column_order": "r_j x_i (j=0..5,i=0..3), then x_i r_j in same loop order",
        "quadratic_relations": RELATION_NAMES,
        "relations_rank": 6,
        "cubic_matrix_shape": [64, 48],
        "minor_determinant": 512,
        "cubic_rank": 45,
        "quadratic_envelope_degree3_dimension": 19,
        "AS4_polynomial_expected_degree3_dimension": 20,
        "euler_object_order": ["A", "P", "A(C0)", "A(f)"],
        "euler_matrix": [[int(v) for v in euler.row(i)] for i in range(4)],
        "euler_determinant": int(euler.det()),
    }


def verify_certificate(relations, matrix, certificate, euler):
    """Check the stored witnesses, without finding replacements for them."""
    assert relations.shape == (16, 6) and relations.rank() == 6
    assert matrix.shape == (64, 48)
    for key, expected in metadata(euler).items():
        assert certificate.get(key) == expected, f"Certificate field differs: {key}"
    rows, columns = certificate["minor_rows"], certificate["minor_columns"]
    for indices, bound in ((rows, 64), (columns, 48)):
        assert isinstance(indices, list) and len(indices) == 45
        assert all(type(index) is int and 0 <= index < bound for index in indices)
        assert len(set(indices)) == 45
    determinant = matrix.extract(rows, columns).det(method="domain-ge")
    assert determinant == certificate["minor_determinant"] == 512
    vectors = certificate["right_kernel_vectors"]
    assert isinstance(vectors, list) and len(vectors) == 3
    for vector in vectors:
        assert isinstance(vector, list) and len(vector) == 48
        assert all(type(value) is int for value in vector)
    kernel = sp.Matrix.hstack(*[sp.Matrix(vector) for vector in vectors])
    assert kernel.shape == (48, 3)
    assert matrix * kernel == sp.zeros(64, 3)
    assert kernel.rank() == 3
    # The nonzero minor gives rank >=45; the kernel gives rank <=48-3.
    assert len(rows) == matrix.cols - kernel.cols == certificate["cubic_rank"]
    assert matrix.rows - certificate["cubic_rank"] == 19


def regenerate_certificate(matrix, euler):
    """Find witnesses only when regeneration is explicitly requested."""
    columns = list(matrix.rref()[1])
    rows = list(matrix[:, columns].T.rref()[1])
    assert len(columns) == len(rows) == 45
    assert matrix.extract(rows, columns).det(method="domain-ge") == 512
    nullvectors = []
    for vector in matrix.nullspace():
        denominator = sp.ilcm(*[v.q for v in vector])
        vector = vector * denominator
        divisor = sp.igcd(*list(vector))
        nullvectors.append(vector / divisor)
    certificate = metadata(euler)
    certificate.update({
        "minor_rows": rows,
        "minor_columns": columns,
        "right_kernel_vectors": [[int(v) for v in vector]
                                  for vector in nullvectors],
    })
    return certificate


def main():
    if not __debug__:
        raise SystemExit("Verification requires assertions; do not use optimized Python.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regenerate", action="store_true",
                        help="explicitly recompute and replace the saved certificate")
    args = parser.parse_args()
    relations, matrix = construct()
    euler = check_euler()
    if args.regenerate:
        certificate = regenerate_certificate(matrix, euler)
    else:
        certificate = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    verify_certificate(relations, matrix, certificate, euler)
    if args.regenerate:
        CERTIFICATE.parent.mkdir(parents=True, exist_ok=True)
        CERTIFICATE.write_text(json.dumps(certificate, indent=2) + "\n", encoding="utf-8")
        print("Regenerated and verified:", CERTIFICATE.relative_to(ROOT))
    else:
        print("Verified without modification:", CERTIFICATE.relative_to(ROOT))
    print("Relations: rank 6; cubic matrix: 64 x 48; minor determinant: 512.")
    print("Independent integral kernel vectors: 3; cubic rank: 45; quotient dimension: 19.")
    print("Euler determinant: 100; exceptional-object numerical formulas: verified.")
    print("PASS: exact degree-three quadratic obstruction and Euler checks.")


if __name__ == "__main__":
    main()
