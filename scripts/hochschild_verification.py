#!/usr/bin/env python3
"""Exact arithmetic certificates for the HH calculation of two Clifford orders.

Run: python hochschild_verification.py
Requires Python 3 and SymPy. All computations take place over QQ.

The mathematical input (smoothness, the conic-bundle decomposition, HKR,
Lefschetz, and the blowup formula) is proved or cited in the accompanying note.
This script verifies the intersection-theoretic and numerical calculations;
it is not a formal verification of those geometric theorems.
"""

import sympy as sp

h, a, b = sp.symbols("h a b")


def homogeneous_part(poly, variables, degree):
    return sp.expand(sum(
        coefficient * sp.prod(v ** exponent for v, exponent in zip(variables, powers))
        for powers, coefficient in sp.Poly(sp.expand(poly), *variables, domain=sp.QQ).terms()
        if sum(powers) == degree
    ))


def reduce_product_ring(poly):
    """Reduce in QQ[h,a,b]/(h^3,a^2,b^2)."""
    return sp.expand(sum(
        coefficient * h ** powers[0] * a ** powers[1] * b ** powers[2]
        for powers, coefficient in sp.Poly(sp.expand(poly), h, a, b, domain=sp.QQ).terms()
        if powers[0] < 3 and powers[1] < 2 and powers[2] < 2
    ))


# Degree 2: X is a divisor Q=2h+a+b in P^2 x P^1 x P^1.
Q = 2*h + a + b
c_ambient = (1+h)**3 * (1+a)**2 * (1+b)**2
c3_x2 = reduce_product_ring(homogeneous_part(
    c_ambient * sum((-Q)**j for j in range(4)), (h,a,b), 3
))
e_x2_chern = sp.Poly(sp.expand(Q*c3_x2), h,a,b).coeff_monomial(h**2*a*b)
assert c3_x2 == 2*a*b*h - a*h**2 - b*h**2
assert e_x2_chern == 2
print("Degree 2 intersection ring: QQ[h,a,b]/(h^3,a^2,b^2), integral(h^2*a*b)=1")
print("c3(TX_2) =", c3_x2)
print("integral_X2 c3 =", e_x2_chern)

# Genus of a (3,3) curve: 2g-2 = D(D+K), K=-2a-2b.
D2 = 3*a+3*b
two_g2_minus_two = sp.Poly(sp.expand(D2*(D2-2*a-2*b)), a,b).coeff_monomial(a*b)
g2 = 1 + two_g2_minus_two/2
assert g2 == 4

# Degree 1: cubic threefold Y in P^4, then blow up one P^1.
c3_y = homogeneous_part((1+h)**5 * sum((-3*h)**j for j in range(4)), (h,), 3)
e_y = 3*sp.Poly(c3_y,h).coeff_monomial(h**3)
e_x1_blowup = e_y + 2  # e(Bl_{P1}Y)=e(Y)+e(P1).
assert c3_y == -2*h**3
assert e_y == -6 and e_x1_blowup == -4
print("Cubic threefold: c3(TY) =", c3_y, "; integral_Y h^3=3")
print("e(Y) =", e_y, "; e(Bl_P1 Y) =", e_x1_blowup)

# Genus of a plane quintic: 2g-2 = 5h(5h-3h).
g1 = 1 + sp.Rational(5*(5-3),2)
assert g1 == 6

cases = [
    {"degree":2, "genus":g2, "e_base":4, "b2_total":3, "e_chern":e_x2_chern},
    {"degree":1, "genus":g1, "e_base":3, "b2_total":2, "e_chern":e_x1_blowup},
]
for case in cases:
    e_total = 2*case["e_base"] + 2 - 2*case["genus"]
    assert e_total == case["e_chern"]
    b3 = 2 + 2*case["b2_total"] - e_total
    h21 = b3/2
    hh0 = (2 + 2*case["b2_total"]) - case["e_base"]
    assert h21 == case["genus"]-1
    assert hh0 == case["e_base"]
    hh = {-1:h21, 0:hh0, 1:h21}
    assert hh == ({-1:3,0:4,1:3} if case["degree"]==2 else {-1:5,0:3,1:5})
    print("degree =", case["degree"],
          "; genus(D) =", case["genus"],
          "; e(X) =", e_total,
          "; b3(X) =", b3,
          "; h21(X) =", h21,
          "; HH_i dimensions =", hh,
          "; all |i|>=2 vanish",
          "; HH Euler characteristic =", hh0-2*h21)

print("All exact arithmetic checks passed.")

