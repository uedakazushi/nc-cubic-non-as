#!/usr/bin/env python3
"""Exact support enumeration and algebraic identities for the matching obstruction.
The universal formal lifting and categorical obstruction are proved in the paper.
"""
from collections import Counter
from itertools import product, permutations
from pathlib import Path
from tempfile import TemporaryDirectory
import os
import subprocess
from sympy import Matrix, symbols, cancel

I = range(3)
TRIPLES = list(product(I, repeat=3))
PERMS = list(permutations(I))


def isolated_matching(mask):
    zero = {a for n, a in enumerate(TRIPLES) if not (mask >> n & 1)}
    isolated = {a for a in zero if all(sum(x != y for x, y in zip(a, b)) != 1 for b in zero)}
    return next((tuple((i, p[i], q[i]) for i in I)
                 for p in PERMS for q in PERMS
                 if all((i, p[i], q[i]) in isolated for i in I)), None)


def check_supports():
    with TemporaryDirectory() as tmp:
        exe = Path(tmp) / 'orbits'
        src = Path(__file__).with_name('ou_zero_supports.cpp')
        subprocess.run([os.environ.get('CXX', 'c++'), '-O3', '-std=c++17',
                        str(src), '-o', str(exe)], check=True)
        rows = [tuple(map(int, line.split())) for line in
                subprocess.check_output([str(exe)], text=True).splitlines()]
    expected = {3: [133947135], 4: [133136111],
                5: [125763310, 132086511], 6: [125759214, 131037935],
                7: [117368558, 125234926], 8: [108979950], 9: [108455662]}
    matches = {}
    for z, mask in rows:
        assert z == 27 - mask.bit_count()
        match = isolated_matching(mask)
        if match is not None:
            # Relabel the old triples to (i,i,i); all exactly-two-equal
            # coefficients must be nonzero in this orientation.
            c = {(i,j,k): (mask >> (9*match[i][0]+3*match[j][1]+match[k][2])) & 1
                 for i,j,k in TRIPLES}
            assert all(c[i,i,i] == 0 for i in I)
            assert all(c[a] for a in TRIPLES if len(set(a)) == 2)
            matches.setdefault(z, []).append(mask)
    assert {z: sorted(ms) for z, ms in matches.items()} == expected
    assert Counter(z for z, _ in rows) == dict(enumerate((1,1,3,9,30,80,201,344,387,217,53,4,1)))
    print('PASS isolated matching supports:', {z: len(ms) for z, ms in sorted(matches.items())})


def check_odd_matrix():
    a,b,c,d,e,f = symbols('a b c d e f', nonzero=True)
    O = Matrix([[0,a,b],[c,0,d],[e,f,0]])
    assert O.det() == a*d*e+b*c*f
    O = O.subs(f, -a*d*e/(b*c))
    right = Matrix([-d/c,-b/a,1])
    left = Matrix([[-d/b,1,-c/e]])
    assert all(cancel(x)==0 for x in O*right)
    assert all(cancel(x)==0 for x in left*O)
    assert O[:2,:2].det() == -a*c
    # Solve O*p = y with p_2=0, when the left-kernel condition holds.
    y0,y1 = symbols('y0 y1')
    y = Matrix([y0,y1,e*y1/c-a*d*e/(b*c)*y0/a])
    p = Matrix([y1/c,y0/a,0])
    assert all(cancel(x)==0 for x in O*p-y)
    # The leading cubic equation on the rank-two locus has nonzero
    # constant term, because all coordinates of the right kernel are units.
    lam,t,h,p0,p1,p2 = symbols('lam t h p0 p1 p2')
    vec = right+t*h*Matrix([p0,p1,p2])
    H = lam*h*h-vec[0]*vec[1]*vec[2]
    assert cancel(H.subs(t,0)-(lam*h*h-b*d/(a*c)))==0
    assert cancel(H.diff(h).subs(t,0)-2*lam*h)==0
    print('PASS odd-column rank-two kernel and lifting identities')


if __name__ == '__main__':
    check_supports()
    check_odd_matrix()
