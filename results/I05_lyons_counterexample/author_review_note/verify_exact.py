"""Self-contained exact check for the five-point counterexample.

Determinants are expanded directly over permutations. Every asserted inequality
uses fractions.Fraction; floating point appears only in the final display.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations


N = 5
TERMS = 24


def permutation_sign(p):
    inversions = sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))
    return -1 if inversions % 2 else 1


PERMS = [(p, permutation_sign(p)) for p in permutations(range(N))]


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
        term = F(permutation_sign(sigma))
        for i in range(n):
            term *= M[i][sigma[i]]
        ans += term
    return ans


def determinant_second_jet(M, B):
    """Second derivative at zero of det(M+i*t*B), as a rational."""
    ans = F(0)
    for sigma, sign in PERMS:
        for a, b in combinations(range(N), 2):
            term = F(-sign) * B[a][sigma[a]] * B[b][sigma[b]]
            for i in range(N):
                if i != a and i != b:
                    term *= M[i][sigma[i]]
            ans += term
    return 2 * ans


def complex_multiply(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def determinant_at_imaginary_step(M, B, t):
    """Exact determinant of M+i*t*B, returned as (real, imaginary)."""
    ans = (F(0), F(0))
    for sigma, sign in PERMS:
        term = (F(sign), F(0))
        for i in range(N):
            term = complex_multiply(term, (M[i][sigma[i]], t * B[i][sigma[i]]))
        ans = ans[0] + term[0], ans[1] + term[1]
    return ans


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


def construction_at(delta):
    P, B = projection_and_direction()
    K = [[(1 - 2 * delta) * P[i][j] + (delta if i == j else 0)
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


def verify_boundary_data():
    P, B = projection_and_direction()
    Q = [[(1 if i == j else 0) - P[i][j] for j in range(N)] for i in range(N)]
    rows = []
    sum_b = F(0)
    sum_wb = F(0)
    log_exponents = defaultdict(F)
    for mask in range(1 << N):
        S = [i for i in range(N) if mask & (1 << i)]
        C = [i for i in range(N) if i not in S]
        w = abs(len(S) - 2)
        a = (determinant_small(principal(P, S)) if len(S) <= 2
             else determinant_small(principal(Q, C)))
        M = [[P[i][j] - (1 if i == j and i in C else 0)
              for j in range(N)] for i in range(N)]
        event_sign = -1 if len(C) % 2 else 1
        b = event_sign * determinant_second_jet(M, B)
        assert a > 0
        rows.append((mask, len(S), w, a, b))
        sum_b += b
        sum_wb += w * b
        for prime, exponent in rational_prime_exponents(a).items():
            log_exponents[prime] += -b * exponent
    assert sum_b == 0
    assert sum_wb == 0
    expected = {
        2: F(260, 5),
        3: F(143, 5),
        5: F(350, 5),
        7: F(-87, 5),
        13: F(-284, 5),
    }
    assert {p: e for p, e in log_exponents.items() if e} == expected
    numerator = 2**260 * 3**143 * 5**350
    denominator = 7**87 * 13**284
    assert numerator > denominator
    return rows, len(str(numerator)), len(str(denominator))


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


def interval_add(x, y):
    return x[0] + y[0], x[1] + y[1]


def interval_subtract(x, y):
    return x[0] - y[1], x[1] - y[0]


def interval_multiply(x, y):
    products = (x[0] * y[0], x[0] * y[1], x[1] * y[0], x[1] * y[1])
    return min(products), max(products)


def inverse(M):
    n = len(M)
    augmented = [
        list(M[i]) + [F(1 if i == j else 0) for j in range(n)]
        for i in range(n)
    ]
    for column in range(n):
        pivot = next(row for row in range(column, n) if augmented[row][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        value = augmented[column][column]
        augmented[column] = [entry / value for entry in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    augmented[row][j] - factor * augmented[column][j]
                    for j in range(2 * n)
                ]
    return [row[n:] for row in augmented]


def matrix_product(A, B):
    return [
        [sum((A[i][k] * B[k][j] for k in range(len(B))), F(0))
         for j in range(len(B[0]))]
        for i in range(len(A))
    ]


def matrix_trace(A):
    return sum((A[i][i] for i in range(len(A))), F(0))


def matrix_trace_product(A, B):
    return sum(
        (A[i][j] * B[j][i] for i in range(len(A)) for j in range(len(A))),
        F(0),
    )


def real_symmetric_basis():
    basis = []
    for i in range(N):
        A = [[F(0) for _ in range(N)] for _ in range(N)]
        A[i][i] = F(1)
        basis.append(A)
    for i in range(N):
        for j in range(i + 1, N):
            A = [[F(0) for _ in range(N)] for _ in range(N)]
            A[i][j] = A[j][i] = F(1)
            basis.append(A)
    return basis


def real_hessian_intervals(K, basis=None):
    if basis is None:
        basis = real_symmetric_basis()
    size = len(basis)
    hessian = [[(F(0), F(0)) for _ in range(size)] for _ in range(size)]
    for mask in range(1 << N):
        S = {i for i in range(N) if mask & (1 << i)}
        M = [[K[i][j] - (1 if i == j and i not in S else 0)
              for j in range(N)] for i in range(N)]
        event_sign = -1 if (N - len(S)) % 2 else 1
        det_M = determinant(M)
        p = event_sign * det_M
        assert p > 0
        R = inverse(M)
        RA = [matrix_product(R, A) for A in basis]
        traces = [matrix_trace(X) for X in RA]
        log_p = log_interval(p)
        for a in range(size):
            for b in range(a + 1):
                trace_product = matrix_trace_product(RA[a], RA[b])
                dp_a = event_sign * det_M * traces[a]
                dp_b = event_sign * det_M * traces[b]
                d2p = event_sign * det_M * (
                    traces[a] * traces[b] - trace_product
                )
                fisher = -dp_a * dp_b / p
                logarithmic = scale(-d2p, log_p)
                entry = interval_add((fisher, fisher), logarithmic)
                hessian[a][b] = interval_add(hessian[a][b], entry)
                hessian[b][a] = hessian[a][b]
    return hessian


REAL_PRECONDITIONER_INTEGERS = [
    [155038, -170265, -37151, -6572, -43223, -351531, 724595, 46576, -196127, -10491487, 119458, 6854552, -1438924, 3875672, -55116],
    [0, 216285, -7537, -2289, 24716, 349639, -1189199, -116933, 335025, 10457461, 5971466, -1147380, -2147004, 112437, -68825],
    [0, 0, 142813, -29088, -18231, -15759, 463402, -53242, -150430, 2671318, -8146990, -803905, 4206898, 2445879, -4531096],
    [0, 0, 0, 149501, -158435, -24198, -44897, -11954, -73423, -1261, 2040691, 1898657, 7803237, -891984, 9119164],
    [0, 0, 0, 0, 237315, 41850, 46099, 135553, 84955, -2636035, 15376, -6801926, -8424212, -5542001, -4464123],
    [0, 0, 0, 0, 0, 170654, -563863, -38505, 193245, 2593310, 4023787, -2167899, -1562645, 2965712, 1129663],
    [0, 0, 0, 0, 0, 0, 713810, -2426, -148199, -1269388, -8107494, -2491582, -365644, -3118522, 2479487],
    [0, 0, 0, 0, 0, 0, 0, 318525, -62668, -5271321, 4113920, -5546196, 2791303, -3968338, -2349331],
    [0, 0, 0, 0, 0, 0, 0, 0, 232670, -1323464, 3072, -3400475, -4210244, 4751100, 4711411],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 5277273, -4114421, 1287134, -75276, -2832317, -1105152],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4085270, 3800443, 2106227, -5973178, -2446136],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8523092, -5434776, -1442067, -613046],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6752662, 2095140, -2386583],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7525038, -5784705],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6365158],
]


def preconditioned_real_basis():
    original = real_symmetric_basis()
    size = len(original)
    assert all(REAL_PRECONDITIONER_INTEGERS[i][i] for i in range(size))
    transformed = []
    for column in range(size):
        A = [[F(0) for _ in range(N)] for _ in range(N)]
        for row in range(size):
            coefficient = F(REAL_PRECONDITIONER_INTEGERS[row][column], 10**8)
            if coefficient:
                for i in range(N):
                    for j in range(N):
                        A[i][j] += coefficient * original[row][i][j]
        transformed.append(A)
    return transformed


def verify_real_negative_definite(K):
    """Exact Gershgorin certificate after a rational congruence."""
    hessian = real_hessian_intervals(K, preconditioned_real_basis())
    margins = []
    for i in range(len(hessian)):
        diagonal_lower = -hessian[i][i][1]
        radius = sum(
            max(abs(hessian[i][j][0]), abs(hessian[i][j][1]))
            for j in range(len(hessian))
            if j != i
        )
        margins.append(diagonal_lower - radius)
    assert min(margins) > 0
    return min(margins)


def event_data(K, B):
    rows = []
    for mask in range(1 << N):
        S = {i for i in range(N) if mask & (1 << i)}
        M = [[K[i][j] - (1 if i == j and i not in S else 0)
              for j in range(N)] for i in range(N)]
        event_sign = -1 if (N - len(S)) % 2 else 1
        p = event_sign * determinant(M)
        b = event_sign * determinant_second_jet(M, B)
        assert p > 0
        rows.append((p, b))
    assert sum((p for p, _ in rows), F(0)) == 1
    assert sum((b for _, b in rows), F(0)) == 0
    return rows


def hessian_interval(K, B):
    lower = upper = F(0)
    for p, b in event_data(K, B):
        lo, hi = scale(-b, log_interval(p))
        lower += lo
        upper += hi
    return lower, upper


def event_probabilities_at_step(K, B, t):
    probabilities = []
    for mask in range(1 << N):
        S = {i for i in range(N) if mask & (1 << i)}
        M = [[K[i][j] - (1 if i == j and i not in S else 0)
              for j in range(N)] for i in range(N)]
        event_sign = -1 if (N - len(S)) % 2 else 1
        real, imaginary = determinant_at_imaginary_step(M, B, t)
        assert imaginary == 0
        p = event_sign * real
        assert p > 0
        probabilities.append(p)
    assert sum(probabilities, F(0)) == 1
    return probabilities


def entropy_interval(probabilities):
    lower = upper = F(0)
    for p in probabilities:
        lo, hi = scale(-p, log_interval(p))
        lower += lo
        upper += hi
    return lower, upper


def main():
    boundary_rows, numerator_digits, denominator_digits = verify_boundary_data()
    assert len(boundary_rows) == 32

    K, B = construction_at(F(1, 10**6))
    lower, upper = hessian_interval(K, B)
    assert F(1, 2) < lower < upper < F(3, 5)
    assert F(5661274627785100258893, 10**22) < lower
    assert upper < F(5661274627785100258923, 10**22)

    K4, _ = construction_at(F(1, 10**4))
    K5, _ = construction_at(F(1, 10**5))
    assert hessian_interval(K4, B)[1] < 0
    assert hessian_interval(K5, B)[0] > 0

    t = F(1, 10**9)
    center_probabilities = event_probabilities_at_step(K, B, F(0))
    plus_probabilities = event_probabilities_at_step(K, B, t)
    minus_probabilities = event_probabilities_at_step(K, B, -t)
    assert plus_probabilities == minus_probabilities
    center_entropy = entropy_interval(center_probabilities)
    endpoint_entropy = entropy_interval(plus_probabilities)
    gap = (
        endpoint_entropy[0] - center_entropy[1],
        endpoint_entropy[1] - center_entropy[0],
    )
    assert gap[0] > F(1, 10**19)
    assert F(281518, 10**24) < gap[0]
    assert gap[1] < F(281520, 10**24)

    frobenius_norm_squared = sum(B[i][j] ** 2 for i in range(N) for j in range(N))
    assert frobenius_norm_squared == 170
    epsilon = F(1, 10**6)
    assert t * t * frobenius_norm_squared < epsilon * epsilon

    real_block_margin = verify_real_negative_definite(K)
    assert real_block_margin > 0

    print("direct permutation expansion; no external packages")
    print("32 boundary rows rebuilt; sum b_S = sum w_S b_S = 0")
    print("boundary limit factorization and positivity verified")
    print("boundary ratio digit counts:", numerator_digits, denominator_digits)
    print("proved: 1/2 < D^2 H_K[iB,iB] < 3/5")
    print("proved: Hessian is negative at delta=1e-4 and positive at delta=1e-5")
    print("proved: H(K +/- 1e-9*iB) - H(K) > 1e-19")
    print("proved: endpoint laws agree event by event")
    print("proved: endpoint feasibility from ||B||_op <= ||B||_F = sqrt(170)")
    print("proved: the Hessian is negative definite on the real symmetric block")
    print("real-block rational congruence margin:", float(real_block_margin))
    print("display Hessian interval:", float(lower), float(upper))
    print("display chord-gap interval:", float(gap[0]), float(gap[1]))


if __name__ == "__main__":
    main()
