#!/usr/bin/env python3
"""Exact certificates for the two five-zero complement obstructions."""
from itertools import product
from sympy import Matrix, symbols, factor, cancel
from verify_ou_zero_supports import normalize, I, PAIRS, TRIPLES


def main():
    a, t = symbols('a t')
    u, v, w = symbols('u:3')
    for mask, values, free_expected in (
        (132054766, (a, -1, 1), [(2,1,2),(2,2,1),(2,2,2)]),
        (134135534, (1, a, 1), [(2,1,1),(2,2,1),(2,2,2)]),
    ):
        coeff, K, free = normalize(mask)
        assert free == free_expected
        c = {q: K.to_sympy(z).subs(dict(zip((u,v,w),values)), simultaneous=True)
             for q,z in coeff.items()}
        p = (1,0,2)  # new-to-old indices for j,k
        c = {(i,j,k):c[i,p[j],p[k]] for i,j,k in TRIPLES}
        assert c[0,0,0] == 0
        assert all(c[q] == 1 for q in product((1,2), repeat=3))
        # Isolated zero and three complementary face determinants.
        assert all(c[q] != 0 for q in TRIPLES if sum(x != 0 for x in q) == 1)
        dets = []
        for h in I:
            ind = [(0,) if h == l else (1,2) for l in I]
            dets.append(factor(Matrix(2,2,[c[q] for q in product(*ind)]).det()))
        assert dets == ([-1,a,-1] if mask == 132054766 else [-1,1,a])
        C = [Matrix(3,3,lambda i,j:c[i,j,k]) for k in I]
        if mask == 132054766:
            A = Matrix([[1,t,a*t],[t,1,1],[(a+1)*t/a-t**3,1,1-a*t**2]])
        else:
            A = Matrix([[1,0,t],[t,1,1],[t,1,1]])
        A0 = Matrix([[1,0,0],[0,1,1],[0,1,1]])
        B0 = Matrix([[1,0,0],[0,1,1],[0,-1,-1]])
        assert A.subs(t,0) == A0
        B = []
        for k in I:
            M = C[k].multiply_elementwise(A)
            assert cancel(M.det()) == 0
            # Normalize the kernel by coordinate 0 for k=0, 1 otherwise.
            rows = (1,2) if k == 0 else (0,1)
            j = 0 if k == 0 else 1
            vker = M.row(rows[0]).cross(M.row(rows[1])).T
            assert cancel(vker[j].subs(t,0)) == (dets[2] if k == 0 else -1)
            vker = (vker/vker[j]).applyfunc(cancel)
            assert vker.subs(t,0) == B0[:,k]
            assert (M*vker).applyfunc(cancel) == Matrix.zeros(3,1)
            B.append(vker)
        B = Matrix.hstack(*B)
        assert all(cancel(sum(c[i,j,k]*A[i,j]*B[j,k] for j in I)) == 0
                   for i,k in PAIRS)
        assert A[0,2] in (a*t,t) and A[1,0] == t
        # Transposed reversed arrow matrices also have two cross entries of order one.
        assert cancel(B[1,0]/t).subs(t,0) == -1
        assert cancel(B[0,1]/t).subs(t,0) == (a if mask == 132054766 else 1)
        print('PASS', mask, 'face determinants', dets, 'formal arc and head bounds')


if __name__ == '__main__':
    main()
