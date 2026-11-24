#!/usr/bin/env python3
"""Exact finite-over-center certificate for support 116309742 at u=1."""
from functools import lru_cache
from itertools import product

from sympy import QQ
from sympy.polys.matrices import DomainMatrix

from verify_ou_zero_supports import I, TRIPLES, support


def row_matrix(rows, width):
    return DomainMatrix(dict(enumerate(rows)), (len(rows), width), QQ)


class CyclicPaths:
    def __init__(self, mask):
        self.support = set(support(mask))

    @lru_cache(None)
    def space(self, length, rotation, start, end):
        """Return basis paths and the reduction map for appending the last arrow."""
        if length == 0:
            return ([(start,)] if start == end else []), {}
        if length == 1:
            return [(start, end)], {(start,): {0: QQ.one}}

        words = [prefix + (end,) for middle in I
                 for prefix in self.space(length - 1, rotation, start, middle)[0]]
        columns = {word: col for col, word in enumerate(words)}
        rows = []
        turn = (rotation + length - 2) % 3
        for first in I:
            for prefix in self.space(length - 2, rotation, start, first)[0]:
                row = {}
                for middle in I:
                    triple = (first, middle, end)
                    if triple[-turn:] + triple[:-turn] not in self.support:
                        continue
                    basis, reductions = self.space(length - 1, rotation, start, middle)
                    for k, value in reductions[prefix].items():
                        row[columns[basis[k] + (end,)]] = value
                rows.append(row)

        reduced, pivots = row_matrix(rows, len(words)).rref()
        entries = reduced.to_dok()
        pivot_rows = {col: i for i, col in enumerate(pivots)}
        free = [col for col in range(len(words)) if col not in pivot_rows]
        indices = {col: i for i, col in enumerate(free)}
        reductions = {}
        for col, word in enumerate(words):
            if col in indices:
                reductions[word[:-1]] = {indices[col]: QQ.one}
            else:
                pivot = pivot_rows[col]
                reductions[word[:-1]] = {i: -entries[pivot, h] for i, h in enumerate(free)
                                         if (pivot, h) in entries}
        return [words[col] for col in free], reductions

    @lru_cache(None)
    def word(self, rotation, word):
        """Coordinates of an arbitrary path in its endpoint space."""
        if len(word) == 1:
            return {0: QQ.one}
        length = len(word) - 1
        basis, _ = self.space(length - 1, rotation, word[0], word[-2])
        _, reductions = self.space(length, rotation, word[0], word[-1])
        out = {}
        for k, coefficient in self.word(rotation, word[:-1]).items():
            for col, value in reductions[basis[k]].items():
                out[col] = out.get(col, QQ.zero) + coefficient * value
        return {col: value for col, value in out.items() if value}

    def center(self, length):
        """A central element is a sum of loops commuting with every arrow."""
        assert length % 3 == 0
        columns = [(r, w) for r, i in product(I, repeat=2)
                   for w in self.space(length, r, i, i)[0]]
        rows = {}
        for col, (rotation, loop) in enumerate(columns):
            i = loop[0]
            prev = (rotation - 1) % 3
            for j in I:
                for k, value in self.word(rotation, loop + (j,)).items():
                    rows.setdefault((rotation, i, j, k), {})[col] = value
                for k, value in self.word(prev, (j,) + loop).items():
                    row = rows.setdefault((prev, j, i, k), {})
                    row[col] = row.get(col, QQ.zero) - value
        kernel = row_matrix(list(rows.values()), len(columns)).nullspace().to_dod()
        return [{columns[col]: value for col, value in row.items()} for row in kernel.values()]

    def product_rank(self, length, rotation, start, end, generators):
        rows = []
        for degree, central in generators:
            for tail in self.space(length - degree, rotation, start, end)[0]:
                row = {}
                for (turn, loop), value in central.items():
                    if (turn, loop[0]) == (rotation, start):
                        for col, v in self.word(turn, loop + tail[1:]).items():
                            row[col] = row.get(col, QQ.zero) + value * v
                rows.append({col: value for col, value in row.items() if value})
        dimension = len(self.space(length, rotation, start, end)[0])
        return row_matrix(rows, dimension).rank(), dimension


def main():
    paths = CyclicPaths(116309742)
    generators = []
    for degree, expected in ((3, 1), (6, 4)):
        basis = paths.center(degree)
        assert len(basis) == expected
        generators.extend((degree, element) for element in basis)
        print(f"degree {degree}: center dimension {len(basis)}")
    dimensions = []
    for rotation, start, end in TRIPLES:
        rank, dim = paths.product_rank(10, rotation, start, end, generators)
        assert rank == dim, (rotation, start, end, rank, dim)
        dimensions.append(dim)
    print(f"PASS all 27 degree-ten spaces spanned; total dimension {sum(dimensions)}")


if __name__ == "__main__":
    main()
