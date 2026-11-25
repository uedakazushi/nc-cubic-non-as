#!/usr/bin/env python3
"""Matrix-unit checks for the three-step cusp lattice argument over C[[t]]."""
from itertools import product

from sympy import Matrix, eye, expand, symbols, zeros


def ceil_third(n):
    return -((-n)//3)


def verify():
    t = symbols('t')
    indices = range(3)

    # L_i has generators t^ceil((i-p)/3) e_p. J_j maps every L_i into L_{i+j}.
    for p,q in product(indices,repeat=2):
        j4 = ceil_third(4+q-p)
        conductor_sum = min(2+ceil_third(i-p)-ceil_third(i-q) for i in indices)
        assert j4==conductor_sum

    def unit(j,q):
        n,p = divmod(q+j,3)
        matrix = zeros(3)
        matrix[p,q] = t**n
        return matrix

    for j in range(4):
        pairing = Matrix(3,3,lambda a,b:expand((unit(j,a)*unit(3-j,b)).trace()).coeff(t,1))
        assert abs(pairing.det())==1

    x = symbols('x0:3')
    element = sum((a*unit(1,i) for i,a in enumerate(x)),zeros(3))
    assert element**3==t*x[0]*x[1]*x[2]*eye(3)
    shift = sum((unit(2,i) for i in indices),zeros(3))
    for q in indices:
        assert shift*unit(4,q)==unit(6,q)
    print('PASS conductor lattices, four perfect pairings, cubic trace and degree-six products')


if __name__=='__main__':
    verify()
