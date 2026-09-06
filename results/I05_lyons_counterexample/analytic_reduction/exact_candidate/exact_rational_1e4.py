#!/usr/bin/env python3
"""Exact certificate for the denominator-10000 six-point chord.

All determinant arithmetic is over pairs of Fraction objects.  Logarithms are
enclosed by a directed atanh series with an explicit rational tail bound.
"""

from __future__ import annotations

import json
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path


N = 6
sys.set_int_max_str_digits(0)
getcontext().prec = 80
ZERO = (F(0), F(0))
ONE = (F(1), F(0))


def cadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def cmul(x, y):
    return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])


def cscale(q, x):
    return (q*x[0], q*x[1])


def parity(p):
    inv = sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p)))
    return -1 if inv % 2 else 1


PERMUTATIONS = [(p, parity(p)) for p in permutations(range(N))]


def determinant(a):
    total = ZERO
    for p, sign in PERMUTATIONS:
        term = ONE
        for i, j in enumerate(p):
            term = cmul(term, a[i][j])
        total = cadd(total, cscale(F(sign), term))
    return total


def determinant_small(a):
    m = len(a)
    if m == 0:
        return ONE
    total = ZERO
    for p in permutations(range(m)):
        term = ONE
        for i, j in enumerate(p):
            term = cmul(term, a[i][j])
        total = cadd(total, cscale(F(parity(p)), term))
    return total


def matmul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return [list(row) for row in zip(*a)]


def matadd(a, b):
    return [[a[i][j]+b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def matvec(a, x):
    return [sum((a[i][j]*x[j] for j in range(len(x))), F(0)) for i in range(len(a))]


def event_probabilities(K, B, sign):
    out = []
    for mask in range(1 << N):
        complement = {i for i in range(N) if not (mask & (1 << i))}
        M = [[(K[i][j]-F(i == j and i in complement), sign*B[i][j])
              for j in range(N)] for i in range(N)]
        value = cscale(F(-1 if len(complement) % 2 else 1), determinant(M))
        assert value[1] == 0
        assert value[0] > 0
        out.append(value[0])
    assert sum(out, F(0)) == 1
    return out


def event_probabilities_by_mobius(K, B, sign):
    """Independent inclusion-exclusion from DPP principal minors."""
    Q = [[(K[i][j], sign*B[i][j]) for j in range(N)] for i in range(N)]
    inclusion = []
    for mask in range(1 << N):
        idx = [i for i in range(N) if mask & (1 << i)]
        minor = [[Q[i][j] for j in idx] for i in idx]
        value = determinant_small(minor)
        assert value[1] == 0
        inclusion.append(value[0])
    out = []
    for mask in range(1 << N):
        complement = ((1 << N)-1) ^ mask
        total = F(0)
        sub = complement
        while True:
            total += F(-1 if sub.bit_count() % 2 else 1) * inclusion[mask | sub]
            if sub == 0:
                break
            sub = (sub-1) & complement
        out.append(total)
    return out


def interval_add(x, y):
    return (x[0]+y[0], x[1]+y[1])


def interval_scale(q, x):
    return (q*x[0], q*x[1]) if q >= 0 else (q*x[1], q*x[0])


def log_interval(x, terms=50):
    assert x > 0
    reduced, exponent = x, 0
    while reduced < 1:
        reduced *= 2
        exponent -= 1
    while reduced >= 2:
        reduced /= 2
        exponent += 1

    def atanh_bounds(y):
        z = (y-1)/(y+1)
        lower = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
        tail = 2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
        return lower, lower+tail

    return interval_add(atanh_bounds(reduced),
                        interval_scale(F(exponent), atanh_bounds(F(2))))


def entropy_interval_grouped(classes, which):
    total = (F(0), F(0))
    for (pc, pe), masks in classes:
        p = pc if which == "center" else pe
        total = interval_add(total, interval_scale(-F(len(masks))*p, log_interval(p)))
    return total


def frac(q):
    return f"{q.numerator}/{q.denominator}" if q.denominator != 1 else str(q.numerator)


def decimal(q, digits=24):
    return format(Decimal(q.numerator)/Decimal(q.denominator), f".{digits}g")


def subset(mask):
    return [i+1 for i in range(N) if mask & (1 << i)]


def main():
    rawK = [
        [5000, 0, 0, 3298, -1672, 3298],
        [0, 5000, 0, -1672, 3298, 3298],
        [0, 0, 5000, 3298, 3298, -1672],
        [3298, -1672, 3298, 5000, 0, 0],
        [-1672, 3298, 3298, 0, 5000, 0],
        [3298, 3298, -1672, 0, 0, 5000],
    ]
    rawB = [
        [0, 11, -11, 0, 0, 0],
        [-11, 0, 11, 0, 0, 0],
        [11, -11, 0, 0, 0, 0],
        [0, 0, 0, 0, -11, 11],
        [0, 0, 0, 11, 0, -11],
        [0, 0, 0, -11, 11, 0],
    ]
    K = [[F(x, 10000) for x in row] for row in rawK]
    B = [[F(x, 10000) for x in row] for row in rawB]
    assert transpose(K) == K
    assert transpose(B) == [[-x for x in row] for row in B]

    # K = 1/2 [[I,C],[C,I]] and B = diag(S,-S).
    C = [[2*K[i][j+3] for j in range(3)] for i in range(3)]
    S = [[B[i][j] for j in range(3)] for i in range(3)]
    assert matadd(matmul(C, S), matmul(S, C)) == [[F(0)]*3 for _ in range(3)]
    assert matmul(K, B) == matmul(B, K)

    eigenpairs = [
        (F(-497, 500), [F(-1, 2), F(-1, 2), F(1)]),
        (F(1231, 1250), [F(1), F(1), F(1)]),
        (F(497, 500), [F(-1), F(1), F(0)]),
    ]
    for lam, v in eigenpairs:
        assert matvec(C, v) == [lam*x for x in v]
    # The three displayed eigenvectors are independent.
    V = [[eigenpairs[j][1][i] for j in range(3)] for i in range(3)]
    detV = (V[0][0]*(V[1][1]*V[2][2]-V[1][2]*V[2][1])
            -V[0][1]*(V[1][0]*V[2][2]-V[1][2]*V[2][0])
            +V[0][2]*(V[1][0]*V[2][1]-V[1][1]*V[2][0]))
    assert detV != 0
    k_eigenvalues = sorted([(1+c)/2 for c, _ in eigenpairs] +
                           [(1-c)/2 for c, _ in eigenpairs])
    assert k_eigenvalues[0] == F(3, 1000)
    assert k_eigenvalues[-1] == F(997, 1000)

    # S=(11/10000)T and ||T||_op=sqrt(3), because T^T T=3I-J.
    a = F(11, 10000)
    T = [[x/a for x in row] for row in S]
    gram = matmul(transpose(T), T)
    assert gram == [[F(2 if i == j else -1) for j in range(3)] for i in range(3)]
    assert 3*a*a < F(3, 1000)**2

    center = event_probabilities(K, B, F(0))
    plus = event_probabilities(K, B, F(1))
    minus = event_probabilities(K, B, F(-1))
    assert center == event_probabilities_by_mobius(K, B, F(0))
    assert plus == event_probabilities_by_mobius(K, B, F(1))
    assert minus == event_probabilities_by_mobius(K, B, F(-1))
    assert plus == minus

    buckets = {}
    for mask, pair in enumerate(zip(center, plus)):
        buckets.setdefault(pair, []).append(mask)
    classes = sorted(buckets.items(), key=lambda item: item[1][0])

    h_center = entropy_interval_grouped(classes, "center")
    h_endpoint = entropy_interval_grouped(classes, "endpoint")
    gap = (h_endpoint[0]-h_center[1], h_endpoint[1]-h_center[0])
    assert gap[0] > F(1, 4_000_000)
    assert gap[0] > F(1, 3_700_000)
    assert gap[1] < F(1, 3_600_000)

    table = []
    for (pc, pe), masks in classes:
        table.append({
            "multiplicity": len(masks),
            "representative_subset": subset(masks[0]),
            "all_masks_decimal": masks,
            "center_probability": frac(pc),
            "endpoint_probability": frac(pe),
            "center_probability_decimal": decimal(pc),
            "endpoint_probability_decimal": decimal(pe),
        })

    result = {
        "status": "NO_HIT_BOUNDED",
        "exact_counterexample_certified": True,
        "certificate_type": "exact rational determinants plus directed rational log intervals",
        "dimension": N,
        "common_denominator": 10000,
        "K_integer_numerators": rawK,
        "B_integer_numerators": rawB,
        "strict_positive_contractions": True,
        "K_eigenvalues": [frac(x) for x in k_eigenvalues],
        "K_spectral_margin": "3/1000",
        "B_operator_norm_squared": "363/100000000",
        "weyl_margin_squared_comparison": "363/100000000 < 9/1000000",
        "all_64_center_probabilities_positive": True,
        "all_64_endpoint_probabilities_positive": True,
        "center_probabilities_sum": frac(sum(center, F(0))),
        "endpoint_probabilities_sum": frac(sum(plus, F(0))),
        "endpoint_laws_equal_eventwise": plus == minus,
        "independent_mobius_crosscheck": True,
        "minimum_center_probability": frac(min(center)),
        "minimum_endpoint_probability": frac(min(plus)),
        "number_of_joint_probability_classes": len(classes),
        "class_table": table,
        "log_series_terms": 50,
        "gap_lower_fraction": frac(gap[0]),
        "gap_upper_fraction": frac(gap[1]),
        "gap_lower_decimal": decimal(gap[0], 50),
        "gap_upper_decimal": decimal(gap[1], 50),
        "gap_interval_width_decimal": decimal(gap[1]-gap[0], 12),
        "proved_simple_bound": "1/3700000 < gap < 1/3600000",
        "proved_requested_bound": "gap > 1/4000000",
    }
    Path("exact_rational_1e4.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    print("denominator-10000 exact certificate")
    print("proved: Q_+ and Q_- are strict positive contractions")
    print("proved: all 64 center and endpoint atoms are positive and normalized")
    print("proved: endpoint laws agree event by event")
    print(f"joint probability classes: {len(classes)} with multiplicities",
          [len(masks) for _, masks in classes])
    print("proved: 1/3700000 < gap < 1/3600000")
    print("directed gap interval:", decimal(gap[0], 30), decimal(gap[1], 30))
    print("minimum center atom:", decimal(min(center), 20))
    print("minimum endpoint atom:", decimal(min(plus), 20))


if __name__ == "__main__":
    main()
