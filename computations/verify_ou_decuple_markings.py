#!/usr/bin/env python3
"""Exact integral check of all six VdB general-position wall sets."""
from itertools import combinations, permutations
from math import gcd


def add(*vs):
    return tuple(map(sum, zip(*vs)))


def mul(c, v):
    return tuple(c*x for x in v)


def sub(v, w):
    return add(v, mul(-1, w))


def canon(v):
    g = gcd(*v)
    assert g
    return tuple(x//(g if next(x for x in v if x) > 0 else -g) for x in v)


def wedge(v, w):
    return canon(tuple(v[i]*w[j]-v[j]*w[i]
                       for i,j in combinations(range(len(v)), 2)))


def main():
    basis = [tuple(int(i == j) for j in range(9)) for i in range(9)]
    A, m, n = basis[:3], basis[3:6], basis[6:]
    M, N = add(*m), add(*n)
    delta = add(*A, N, mul(-2, M))
    K = [sub(sub(M, mi), nj) for mi in m for nj in n]
    masks = {}
    for bit, (a,b,c) in enumerate(permutations(A)):
        D = add(M, mul(-2,a), mul(-1,c))
        L0 = add(b,mul(-1,a),N,mul(-1,M))
        L1 = sub(mul(2,L0),D)
        assert sub(L1,L0) == delta
        P = [sub(mi,a) for mi in m] + [add(M,mul(-1,a),mul(-1,c),mul(-1,nj)) for nj in n]
        direct = [sub(P[i],P[j]) for i in range(3) for j in range(3,6)]
        direct += [sub(add(*(P[i] for i in I)),L1) for I in combinations(range(6),3)]
        direct += [sub(add(*P),mul(2,L1))]
        walls = [sub(k,c) for k in K] + [sub(k,add(a,mul(2,delta))) for k in K]
        walls += [sub(k,sub(b,mul(3,delta))) for k in K]
        walls += [sub(sub(c,a),mul(2,delta)),sub(sub(b,c),mul(3,delta)),sub(sub(b,a),mul(5,delta))]
        assert len(walls) == 30
        assert all(gcd(*v)==1 and sum(v)==0 and sum(v[3:6])+2*sum(v[6:])==0 for v in walls)
        assert set(map(canon,direct)) == set(map(canon,walls))
        assert len(set(map(canon,walls))) == 30
        for v in map(canon,walls):
            masks[v] = masks.get(v,0) | (1<<bit)
    assert len(masks)==99
    assert sum(b.bit_count()==1 for b in masks.values())==18
    assert sum(b.bit_count()==2 for b in masks.values())==81
    planes={}
    for (v,b),(w,c) in combinations(masks.items(),2):
        p=wedge(v,w)
        planes[p]=planes.get(p,0)|b|c
    bad={p for p,b in planes.items() if b==63}
    expected={wedge(sub(k,a),delta) for k in K for a in A}
    expected|={wedge(sub(a,b),delta) for a,b in combinations(A,2)}
    assert len(planes)==3975 and len(bad)==30 and bad==expected
    print('6 markings; 30 walls per marking; 99 primitive walls; 3975 rank-two planes; 30 covering planes, all containing delta.')


if __name__=='__main__':
    main()
