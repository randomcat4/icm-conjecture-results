#!/usr/bin/env python3
"""Exact certificate for the optimized rational five-point candidate."""

from fractions import Fraction as F
from itertools import permutations


N = 5
ZERO = (F(0), F(0))
ONE = (F(1), F(0))


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def multiply(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def scale(q, x):
    return (q * x[0], q * x[1])


def parity(p):
    return -1 if sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p))) % 2 else 1


def determinant(a):
    total = ZERO
    for p in permutations(range(len(a))):
        term = ONE
        for i, j in enumerate(p):
            term = multiply(term, a[i][j])
        total = add(total, scale(F(parity(p)), term))
    return total


def endpoint(K, B, sign=1):
    return [[(K[i][j], F(sign) * B[i][j]) for j in range(N)] for i in range(N)]


def identity_minus(Q):
    return [[(F(i == j) - Q[i][j][0], -Q[i][j][1]) for j in range(N)] for i in range(N)]


def leading_principal_minors(Q):
    values = []
    for size in range(1, N + 1):
        value = determinant([row[:size] for row in Q[:size]])
        assert value[1] == 0 and value[0] > 0
        values.append(value[0])
    return values


def event_probabilities(Q):
    probabilities = []
    for mask in range(1 << N):
        complement = {i for i in range(N) if not (mask & (1 << i))}
        matrix = [[add(Q[i][j], (F(-1 if i == j and i in complement else 0), F(0)))
                   for j in range(N)] for i in range(N)]
        value = scale(F(-1 if len(complement) % 2 else 1), determinant(matrix))
        assert value[1] == 0 and value[0] > 0
        probabilities.append(value[0])
    assert sum(probabilities, F(0)) == 1
    return probabilities


def interval_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def interval_scale(q, x):
    return (q * x[0], q * x[1]) if q >= 0 else (q * x[1], q * x[0])


def log_interval(x, terms=60):
    reduced, exponent = x, 0
    while reduced < 1:
        reduced *= 2
        exponent -= 1
    while reduced >= 2:
        reduced /= 2
        exponent += 1

    def bounds(y):
        z = (y - 1) / (y + 1)
        lower = 2 * sum((z ** (2 * j + 1) / (2 * j + 1) for j in range(terms)), F(0))
        upper = lower + 2 * z ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z * z))
        return lower, upper

    return interval_add(bounds(reduced), interval_scale(F(exponent), bounds(F(2))))


def entropy_interval(probabilities):
    total = (F(0), F(0))
    for p in probabilities:
        total = interval_add(total, interval_scale(-p, log_interval(p)))
    return total


def main():
    K_numbers = [
        [281472, -264453, 221468, 107027, -263963],
        [-264453, 313361, -297320, -221468, 70061],
        [221468, -297320, 313361, 264453, 70061],
        [107027, -221468, 264453, 281472, 263963],
        [-263963, 70061, 70061, 263963, 813761],
    ]
    B_numbers = [
        [0, -185, -98, -372, -73],
        [185, 0, -175, 98, 299],
        [98, 175, 0, 185, -299],
        [372, -98, -185, 0, -73],
        [73, -299, 299, 73, 0],
    ]
    K = [[F(x, 10**6) for x in row] for row in K_numbers]
    B = [[F(x, 10**6) for x in row] for row in B_numbers]
    assert all(K[i][j] == K[j][i] for i in range(N) for j in range(N))
    assert all(B[i][j] == -B[j][i] for i in range(N) for j in range(N))

    Q_plus = endpoint(K, B, 1)
    Q_minus = endpoint(K, B, -1)
    leading_principal_minors(Q_plus)
    leading_principal_minors(identity_minus(Q_plus))
    leading_principal_minors(Q_minus)
    leading_principal_minors(identity_minus(Q_minus))

    center = event_probabilities(endpoint(K, B, 0))
    plus = event_probabilities(Q_plus)
    minus = event_probabilities(Q_minus)
    assert plus == minus

    center_entropy = entropy_interval(center)
    endpoint_entropy = entropy_interval(plus)
    gap = (endpoint_entropy[0] - center_entropy[1],
           endpoint_entropy[1] - center_entropy[0])
    assert gap[0] > F(1, 130_000_000)
    assert gap[1] < F(1, 120_000_000)

    print("optimized five-point rational construction")
    print("proved: both endpoints are strict positive contractions")
    print("proved: endpoint laws agree event by event")
    print("proved: 1/130000000 < endpoint entropy - midpoint entropy < 1/120000000")
    print("minimum center event probability:", float(min(center)))
    print("minimum endpoint event probability:", float(min(plus)))
    print("directed gap interval:", float(gap[0]), float(gap[1]))


if __name__ == "__main__":
    main()
