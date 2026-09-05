"""Dependency-free exact check, independent of the SymPy discovery code.

Determinants and their t^2 coefficients are expanded directly over permutations.
All inequalities use fractions.Fraction; floating point is used only for display.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path


N = 5
TERMS = 24


def parity(p):
    inversions = sum(p[i] > p[j] for i in range(N) for j in range(i + 1, N))
    return -1 if inversions % 2 else 1


PERMS = [(p, parity(p)) for p in permutations(range(N))]


def determinant(M):
    ans = F(0)
    for sigma, sign in PERMS:
        term = F(sign)
        for i in range(N):
            term *= M[i][sigma[i]]
        ans += term
    return ans


def determinant_small(M):
    n = len(M)
    if n == 0:
        return F(1)
    ans = F(0)
    for sigma in permutations(range(n)):
        inversions = sum(sigma[i] > sigma[j] for i in range(n) for j in range(i + 1, n))
        term = F(-1 if inversions % 2 else 1)
        for i in range(n):
            term *= M[i][sigma[i]]
        ans += term
    return ans


def determinant_second_jet(M, B):
    """Second derivative at zero of det(M+i*t*B), returned as a rational."""
    ans = F(0)
    for sigma, sign in PERMS:
        for a, b in combinations(range(N), 2):
            term = F(-sign) * B[a][sigma[a]] * B[b][sigma[b]]
            for i in range(N):
                if i != a and i != b:
                    term *= M[i][sigma[i]]
            ans += term
    return 2 * ans


def wedge(u, v):
    return [[u[i] * v[j] - v[i] * u[j] for j in range(N)] for i in range(N)]


def add(A, B, alpha=F(1)):
    return [[A[i][j] + alpha * B[i][j] for j in range(N)] for i in range(N)]


def outer(u, v):
    return [[u[i] * v[j] for j in range(N)] for i in range(N)]


def projection_and_direction():
    one = [F(1)] * N
    x = list(map(F, [-2, -1, 0, 1, 2]))
    q2 = list(map(F, [2, -1, -2, -1, 2]))
    q3 = list(map(F, [-1, 2, 0, -2, 1]))
    P = add([[z / 5 for z in row] for row in outer(one, one)],
            [[z / 10 for z in row] for row in outer(x, x)])
    B = add(wedge(one, x), wedge(q2, q3), F(-1, 2))
    return P, B


def construction():
    P, B = projection_and_direction()
    epsilon = F(1, 10**6)
    K = [[(1 - 2 * epsilon) * P[i][j] + (epsilon if i == j else 0)
          for j in range(N)] for i in range(N)]
    return K, B


def principal(M, indices):
    return [[M[i][j] for j in indices] for i in indices]


def factor_integer(n):
    factors = defaultdict(int)
    p = 2
    while p * p <= n:
        while n % p == 0:
            factors[p] += 1
            n //= p
        p += 1
    if n > 1:
        factors[n] += 1
    return factors


def rational_prime_exponents(q):
    ans = defaultdict(int)
    for p, e in factor_integer(q.numerator).items():
        ans[p] += e
    for p, e in factor_integer(q.denominator).items():
        ans[p] -= e
    return ans


def verify_boundary_table():
    P, B = projection_and_direction()
    Q = [[(1 if i == j else 0) - P[i][j] for j in range(N)] for i in range(N)]
    rows = list(csv.DictReader((Path(__file__).parent / "boundary_table.csv").open(encoding="utf-8")))
    assert len(rows) == 1 << N
    sum_b = F(0)
    sum_wb = F(0)
    log_exponents = defaultdict(F)
    for mask, row in enumerate(rows):
        assert int(row["mask"]) == mask
        S = [i for i in range(N) if mask & (1 << i)]
        C = [i for i in range(N) if i not in S]
        w = abs(len(S) - 2)
        a = determinant_small(principal(P, S)) if len(S) <= 2 else determinant_small(principal(Q, C))
        M = [[P[i][j] - (1 if i == j and i in C else 0) for j in range(N)] for i in range(N)]
        event_sign = -1 if len(C) % 2 else 1
        b = event_sign * determinant_second_jet(M, B)
        assert int(row["size"]) == len(S)
        assert int(row["order_w"]) == w
        assert F(row["leading_a"]) == a > 0
        assert F(row["second_jet_b"]) == b
        sum_b += b
        sum_wb += w * b
        for prime, exponent in rational_prime_exponents(a).items():
            log_exponents[prime] += -b * exponent
    assert sum_b == 0
    assert sum_wb == 0
    expected = {2: F(260, 5), 3: F(143, 5), 5: F(350, 5),
                7: F(-87, 5), 13: F(-284, 5)}
    assert {p: e for p, e in log_exponents.items() if e} == expected
    numerator = 2**260 * 3**143 * 5**350
    denominator = 7**87 * 13**284
    assert numerator > denominator
    return len(str(numerator)), len(str(denominator))


def atanh_log_interval(r):
    assert 1 <= r <= 2
    y = (r - 1) / (r + 1)
    lower = 2 * sum((y ** (2 * j + 1) / (2 * j + 1) for j in range(TERMS)), F(0))
    degree = 2 * TERMS + 1
    tail = 2 * y**degree / (degree * (1 - y**2))
    return lower, lower + tail


LOG2 = atanh_log_interval(F(2))


def log_interval(q):
    assert q > 0
    k = q.numerator.bit_length() - q.denominator.bit_length()
    power = F(2) ** k
    while q < power:
        k -= 1
        power /= 2
    while q >= 2 * power:
        k += 1
        power *= 2
    rlo, rhi = atanh_log_interval(q / power)
    if k >= 0:
        return k * LOG2[0] + rlo, k * LOG2[1] + rhi
    return k * LOG2[1] + rlo, k * LOG2[0] + rhi


def scale(c, interval):
    lo, hi = interval
    return (c * lo, c * hi) if c >= 0 else (c * hi, c * lo)


def main():
    numerator_digits, denominator_digits = verify_boundary_table()
    K, B = construction()
    lower = upper = F(0)
    probability_sum = F(0)
    jet_sum = F(0)
    for mask in range(1 << N):
        S = {i for i in range(N) if mask & (1 << i)}
        M = [[K[i][j] - (1 if i == j and i not in S else 0)
              for j in range(N)] for i in range(N)]
        event_sign = -1 if (N - len(S)) % 2 else 1
        p = event_sign * determinant(M)
        b = event_sign * determinant_second_jet(M, B)
        assert p > 0
        probability_sum += p
        jet_sum += b
        lo, hi = scale(-b, log_interval(p))
        lower += lo
        upper += hi
    assert probability_sum == 1
    assert jet_sum == 0
    assert F(1, 2) < lower < upper < F(3, 5)
    assert F(5661274627785100258893, 10**22) < lower
    assert upper < F(5661274627785100258923, 10**22)
    print("direct permutation expansion; no external packages")
    print("boundary table rebuilt; sum b_S = sum w_S b_S = 0")
    print("boundary limit factorization and positivity verified")
    print("boundary ratio digit counts:", numerator_digits, denominator_digits)
    print("sum p_S =", probability_sum)
    print("sum b_S =", jet_sum)
    print("proved: 1/2 < D^2 H_K[iB,iB] < 3/5")
    print("display interval:", float(lower), float(upper))


if __name__ == "__main__":
    main()
