#!/usr/bin/env python3
"""Exact certificate for a symmetric six-point DPP entropy counterexample."""

from fractions import Fraction as F
from itertools import permutations
from collections import Counter


N = 6
ZERO = (F(0), F(0))
ONE = (F(1), F(0))


def cadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def cmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def cscale(q, x):
    return (q * x[0], q * x[1])


def parity(p):
    inversions = sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))
    return -1 if inversions % 2 else 1


PERMUTATIONS = [(p, parity(p)) for p in permutations(range(N))]


def determinant(a):
    total = ZERO
    for p, sign in PERMUTATIONS:
        term = ONE
        for i, j in enumerate(p):
            term = cmul(term, a[i][j])
        total = cadd(total, cscale(F(sign), term))
    return total


def matrix_product(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(N)), F(0))
             for j in range(N)] for i in range(N)]


def transpose(a):
    return [[a[j][i] for j in range(N)] for i in range(N)]


def event_probabilities(K, B, t):
    probabilities = []
    for mask in range(1 << N):
        complement = {i for i in range(N) if not (mask & (1 << i))}
        matrix = [[
            (K[i][j] - F(i == j and i in complement), t * B[i][j])
            for j in range(N)] for i in range(N)]
        value = cscale(F(-1 if len(complement) % 2 else 1), determinant(matrix))
        assert value[1] == 0
        assert value[0] > 0
        probabilities.append(value[0])
    assert sum(probabilities, F(0)) == 1
    return probabilities


def interval_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def interval_scale(q, x):
    return (q * x[0], q * x[1]) if q >= 0 else (q * x[1], q * x[0])


def log_interval(x, terms=50):
    assert x > 0
    reduced = x
    exponent = 0
    while reduced < 1:
        reduced *= 2
        exponent -= 1
    while reduced >= 2:
        reduced /= 2
        exponent += 1

    def atanh_bounds(y):
        z = (y - 1) / (y + 1)
        lower = 2 * sum((z ** (2 * j + 1) / (2 * j + 1) for j in range(terms)), F(0))
        upper = lower + 2 * z ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z * z))
        return lower, upper

    base = atanh_bounds(reduced)
    log_two = atanh_bounds(F(2))
    return interval_add(base, interval_scale(F(exponent), log_two))


def entropy_interval(probabilities):
    total = (F(0), F(0))
    for p in probabilities:
        total = interval_add(total, interval_scale(-p, log_interval(p)))
    return total


def main():
    R = [
        [F(2, 3), F(-1, 3), F(2, 3)],
        [F(-1, 3), F(2, 3), F(2, 3)],
        [F(2, 3), F(2, 3), F(-1, 3)],
    ]
    P = [[F(0) for _ in range(N)] for _ in range(N)]
    for i in range(3):
        P[i][i] = P[i + 3][i + 3] = F(1, 2)
        for j in range(3):
            P[i][j + 3] = P[i + 3][j] = R[i][j] / 2

    S = [[F(0), F(1), F(-1)], [F(-1), F(0), F(1)], [F(1), F(-1), F(0)]]
    B = [[F(0) for _ in range(N)] for _ in range(N)]
    for i in range(3):
        for j in range(3):
            B[i][j] = S[i][j]
            B[i + 3][j + 3] = -S[i][j]

    I = [[F(i == j) for j in range(N)] for i in range(N)]
    assert transpose(P) == P
    assert matrix_product(P, P) == P
    assert transpose(B) == [[-x for x in row] for row in B]
    assert matrix_product(P, B) == matrix_product(B, P)

    epsilon = F(1, 475)
    t = F(1, 2300)
    K = [[epsilon * I[i][j] + (1 - 2 * epsilon) * P[i][j]
          for j in range(N)] for i in range(N)]

    # K and I-K have spectral margin epsilon.  The Frobenius bound gives
    # ||t iB||_op^2 <= t^2 ||B||_F^2 = 12 t^2 < epsilon^2.
    frobenius_squared = sum(x * x for row in B for x in row)
    assert frobenius_squared == 12
    assert frobenius_squared * t * t < epsilon * epsilon

    center = event_probabilities(K, B, F(0))
    plus = event_probabilities(K, B, t)
    minus = event_probabilities(K, B, -t)
    assert plus == minus

    # Exact eight-row table in the note.  Each triple is
    # (multiplicity, center probability, endpoint change).
    u = epsilon * (1 - epsilon)
    z = t * t
    rows = [
        (2, u**3, 3*z*u*(3*z + 2*u - 1)),
        (12, u**2*(1 - 2*u)/2,
         -z*(18*z*u - 3*z + 12*u**2 - 6*u + 1)/2),
        (12, u*(1 - 2*u)**2/4,
         z*(18*z*u - 6*z + 12*u**2 - 8*u + 1)/2),
        (2, (1 - 2*u)**3/8,
         3*z*(1 - 2*u)*(3*z + 2*u - 1)/2),
        (12, u*(36*u**2 - 20*u + 5)/36,
         z*(54*z*u - 12*z + 36*u**2 - 10*u + 1)/6),
        (6, u*(9*u**2 - 8*u + 2)/9,
         z*(27*z*u - 6*z + 18*u**2 - 11*u + 2)/3),
        (12, (1 - 2*u)*(9*u**2 - 4*u + 1)/18,
         -z*(54*z*u - 15*z + 36*u**2 - 20*u + 2)/6),
        (6, (1 - 2*u)*(36*u**2 - 4*u + 1)/72,
         -z*(54*z*u - 15*z + 36*u**2 - 8*u - 1)/6),
    ]
    expected_groups = Counter({(p, p + d): multiplicity
                               for multiplicity, p, d in rows})
    assert Counter(zip(center, plus)) == expected_groups
    representatives = [0, 1, 3, 7, 9, 10, 11, 13]
    for mask, (_, p, d) in zip(representatives, rows):
        assert center[mask] == p
        assert plus[mask] == p + d

    center_entropy = entropy_interval(center)
    endpoint_entropy = entropy_interval(plus)
    gap = (endpoint_entropy[0] - center_entropy[1],
           endpoint_entropy[1] - center_entropy[0])
    assert gap[0] > F(1, 40_000_000)
    assert gap[1] < F(1, 39_000_000)
    assert center_entropy[0] > F(2_877_726_131_217_012, 10**15)
    assert center_entropy[1] < F(2_877_726_131_217_013, 10**15)

    print("six-point symmetric construction")
    print("epsilon = 1/475; t = 1/2300")
    print("proved: Q_+ and Q_- are strict positive contractions")
    print("proved: endpoint laws agree event by event")
    print("proved: the 64 events agree with the eight-row table")
    print("proved: 1/40000000 < endpoint entropy - midpoint entropy < 1/39000000")
    print("proved: 2.877726131217012 < H(K) < 2.877726131217013")
    print("minimum center event probability:", float(min(center)))
    print("minimum endpoint event probability:", float(min(plus)))
    print("directed gap interval:", float(gap[0]), float(gap[1]))


if __name__ == "__main__":
    main()
