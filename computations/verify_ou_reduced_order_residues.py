"""Residue constraints at the possible rank-(1,2) junction (characteristic zero)."""

from itertools import combinations_with_replacement, permutations, product

from sympy import Matrix, symbols


def transported(types, action):
    return tuple(sorted(action[i] for i in types))


def orthogonality_possible(first, other, count):
    """Test for a pair of nonzero residues. The common factor 1/2 is absorbed."""
    at_bad = (0, 1)
    first_difference = [
        int(first[0] == i) - int(first[-1] == i) for i in range(count)
    ]
    second_difference = [at_bad.count(i) - other.count(i) for i in range(count)]
    if not any(first_difference) or not any(second_difference):
        return not any(first_difference) and not any(second_difference)
    return Matrix([first_difference, second_difference]).rank() == 1


def split_second():
    # The second projection is faithful. A simple supported only on the other
    # branch with dimension two would contribute a nonzero unbalanced residue
    # on a primitive idempotent. All remaining types are one-dimensional.
    counts = dict.fromkeys(product((1, 2), repeat=2), 0)
    for count in range(2, 5):
        labels = tuple(range(count))
        for other in combinations_with_replacement(labels, 2):
            if set((0, 1) + other) != set(labels):
                continue
            for first in ((i,) for i in labels):
                check_first(first, other, labels, counts)
            for first in combinations_with_replacement(labels, 2):
                if first[0] != first[1]:
                    check_first(first, other, labels, counts)
    assert counts[1, 1] > 0 and counts[2, 1] > 0
    assert counts[1, 2] == counts[2, 2] == 0
    print('PASS incompatible residue multiplicities', sorted(counts.items()))


def check_first(first, other, labels, counts):
    for action in permutations(labels):
        if transported(first, action) != first:
            continue
        for period in (1, 2):
            at_bad = (0, 1)
            images = (transported(at_bad, action), transported(other, action))
            expected = (at_bad, other) if period == 1 else (other, at_bad)
            if images != expected:
                continue
            returning = action if period == 1 else tuple(action[i] for i in action)
            if returning[0:2] != (1, 0):
                continue
            key = (len(first), period)
            counts[key] += 1
            assert not orthogonality_possible(first, other, len(labels))


def double_pole():
    t, a, b = symbols('t a b')
    step = Matrix([[0, a], [t*b, 0]])
    assert step**2 == t*a*b*Matrix.eye(2)
    assert (step**2).trace().expand().coeff(t) == 2*a*b
    # In the dual numbers, the analogous square of a radical element vanishes.
    assert ((t*b)**2).as_poly(t).rem((t**2).as_poly(t)).is_zero
    print('PASS degree-two trace on the periodic two-step flag')


if __name__ == '__main__':
    split_second()
    double_pole()
