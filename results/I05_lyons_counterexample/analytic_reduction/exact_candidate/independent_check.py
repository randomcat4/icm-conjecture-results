#!/usr/bin/env python3
"""Independent exact audit of the denominator-10000 DPP certificate.

This file deliberately does not import or execute the candidate program.  It uses
Gaussian elimination over Q(i), whereas the candidate uses a Leibniz expansion.
All probability, definiteness, grouping, and logarithm-bound decisions are exact.
"""

from __future__ import annotations

import csv
import json
import sys
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
N = 6
sys.set_int_max_str_digits(0)
Z = (F(0), F(0))
O = (F(1), F(0))


RAW_K = [
    [5000, 0, 0, 3298, -1672, 3298],
    [0, 5000, 0, -1672, 3298, 3298],
    [0, 0, 5000, 3298, 3298, -1672],
    [3298, -1672, 3298, 5000, 0, 0],
    [-1672, 3298, 3298, 0, 5000, 0],
    [3298, 3298, -1672, 0, 0, 5000],
]
RAW_B = [
    [0, 11, -11, 0, 0, 0],
    [-11, 0, 11, 0, 0, 0],
    [11, -11, 0, 0, 0, 0],
    [0, 0, 0, 0, -11, 11],
    [0, 0, 0, 11, 0, -11],
    [0, 0, 0, -11, 11, 0],
]


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def neg(x):
    return (-x[0], -x[1])


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    return (x[0] * y[0] - x[1] * y[1],
            x[0] * y[1] + x[1] * y[0])


def div(x, y):
    d = y[0] * y[0] + y[1] * y[1]
    assert d != 0
    return ((x[0] * y[0] + x[1] * y[1]) / d,
            (x[1] * y[0] - x[0] * y[1]) / d)


def conj(x):
    return (x[0], -x[1])


def det_gaussian(a):
    """Exact determinant over Q(i) using row elimination."""
    m = len(a)
    if m == 0:
        return O
    b = [[tuple(x) for x in row] for row in a]
    ans = O
    for col in range(m):
        pivot = next((r for r in range(col, m) if b[r][col] != Z), None)
        if pivot is None:
            return Z
        if pivot != col:
            b[col], b[pivot] = b[pivot], b[col]
            ans = neg(ans)
        p = b[col][col]
        ans = mul(ans, p)
        for r in range(col + 1, m):
            if b[r][col] == Z:
                continue
            factor = div(b[r][col], p)
            for c in range(col + 1, m):
                b[r][c] = sub(b[r][c], mul(factor, b[col][c]))
            b[r][col] = Z
    return ans


def principal(a, idx):
    return [[a[i][j] for j in idx] for i in idx]


def real_det(a):
    d = det_gaussian(a)
    assert d[1] == 0
    return d[0]


def qmatrix(sign):
    return [[(F(RAW_K[i][j], 10000), F(sign * RAW_B[i][j], 10000))
             for j in range(N)] for i in range(N)]


def identity_minus(a):
    return [[sub((O if i == j else Z), a[i][j]) for j in range(N)]
            for i in range(N)]


def is_hermitian(a):
    return all(a[i][j] == conj(a[j][i]) for i in range(N) for j in range(N))


def sylvester_pivots(a):
    """Leading principal minors; all positive iff Hermitian a is positive definite."""
    assert is_hermitian(a)
    return [real_det(principal(a, list(range(k)))) for k in range(1, N + 1)]


def complete_atoms_direct(q):
    """Use p_A=(-1)^|A^c| det(Q-I_{A^c}), with exact elimination."""
    out = []
    for mask in range(1 << N):
        ac = [i for i in range(N) if not (mask >> i) & 1]
        m = [[sub(q[i][j], (O if i == j and i in ac else Z))
              for j in range(N)] for i in range(N)]
        d = det_gaussian(m)
        assert d[1] == 0
        out.append((-d[0]) if len(ac) % 2 else d[0])
    return out


def complete_atoms_mobius(q):
    """Independent formula from inclusion probabilities det(Q_S)."""
    inc = []
    for mask in range(1 << N):
        idx = [i for i in range(N) if (mask >> i) & 1]
        inc.append(real_det(principal(q, idx)))
    out = []
    full = (1 << N) - 1
    for mask in range(1 << N):
        rest = full ^ mask
        s = rest
        v = F(0)
        while True:
            v += (-1 if s.bit_count() % 2 else 1) * inc[mask | s]
            if s == 0:
                break
            s = (s - 1) & rest
        out.append(v)
    return out


def iadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def iscale(c, x):
    return (c * x[0], c * x[1]) if c >= 0 else (c * x[1], c * x[0])


def log_bounds_positive_rational(x, terms=60):
    """Exact directed enclosure of natural log(x).

    Write x=2^e*r with 1<=r<2 and use
      log(y)=2 sum_{j>=0} z^(2j+1)/(2j+1), z=(y-1)/(y+1).
    Here z is in [0,1/3].  The omitted positive tail is at most
      2*z^(2m+1)/((2m+1)*(1-z^2)).
    Negative e is handled by interval scalar multiplication, which reverses
    endpoints exactly.
    """
    assert x > 0
    r, e = x, 0
    while r < 1:
        r *= 2
        e -= 1
    while r >= 2:
        r /= 2
        e += 1

    def core(y):
        z = (y - 1) / (y + 1)
        assert 0 <= z <= F(1, 3)
        partial = 2 * sum((z ** (2*j + 1) / (2*j + 1)
                           for j in range(terms)), F(0))
        remainder = (2 * z ** (2*terms + 1)
                     / ((2*terms + 1) * (1 - z*z)))
        return (partial, partial + remainder)

    return iadd(core(r), iscale(F(e), core(F(2))))


def entropy_bounds(ps):
    h = (F(0), F(0))
    for p in ps:
        h = iadd(h, iscale(-p, log_bounds_positive_rational(p)))
    return h


def fstr(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def dstr(x, places=30):
    scale = 10 ** places
    n = x.numerator * scale // x.denominator
    sign = "-" if n < 0 else ""
    n = abs(n)
    return f"{sign}{n // scale}.{n % scale:0{places}d}"


def sci(x, digits=30):
    with localcontext() as ctx:
        ctx.prec = digits
        return format(Decimal(x.numerator) / Decimal(x.denominator), f".{digits - 1}E")


def rmatvec(a, v):
    return [sum((a[i][j] * v[j] for j in range(len(v))), F(0))
            for i in range(len(a))]


def rmatmul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def rtranspose(a):
    return [list(row) for row in zip(*a)]


def subset(mask):
    return "{" + ",".join(str(i + 1) for i in range(N) if (mask >> i) & 1) + "}"


def main():
    k_symmetric = all(RAW_K[i][j] == RAW_K[j][i] for i in range(N) for j in range(N))
    b_skew = all(RAW_B[i][j] == -RAW_B[j][i] for i in range(N) for j in range(N))
    b_zero_diag = all(RAW_B[i][i] == 0 for i in range(N))
    assert k_symmetric and b_skew and b_zero_diag

    # Independently audit the candidate's Weyl/norm argument.
    K = [[F(x, 10000) for x in row] for row in RAW_K]
    B = [[F(x, 10000) for x in row] for row in RAW_B]
    C = [[2 * K[i][j + 3] for j in range(3)] for i in range(3)]
    eigendata = [
        (F(-497, 500), [F(-1, 2), F(-1, 2), F(1)]),
        (F(1231, 1250), [F(1), F(1), F(1)]),
        (F(497, 500), [F(-1), F(1), F(0)]),
    ]
    assert all(rmatvec(C, v) == [lam * x for x in v] for lam, v in eigendata)
    # Determinant of the matrix with these vectors as columns is nonzero.
    V = [[eigendata[j][1][i] for j in range(3)] for i in range(3)]
    detv = (V[0][0] * (V[1][1]*V[2][2] - V[1][2]*V[2][1])
            - V[0][1] * (V[1][0]*V[2][2] - V[1][2]*V[2][0])
            + V[0][2] * (V[1][0]*V[2][1] - V[1][1]*V[2][0]))
    assert detv != 0
    k_eigs = sorted([(1 + lam)/2 for lam, _ in eigendata]
                    + [(1 - lam)/2 for lam, _ in eigendata])
    assert k_eigs[0] == F(3, 1000) and k_eigs[-1] == F(997, 1000)
    S = [[B[i][j] for j in range(3)] for i in range(3)]
    T = [[x / F(11, 10000) for x in row] for row in S]
    gram = rmatmul(rtranspose(T), T)
    assert gram == [[F(2 if i == j else -1) for j in range(3)] for i in range(3)]
    bnorm_sq = F(363, 100_000_000)
    assert bnorm_sq < F(9, 1_000_000)
    assert rmatmul(K, B) == rmatmul(B, K)

    qm, q0, qp = qmatrix(-1), qmatrix(0), qmatrix(1)
    assert all(is_hermitian(q) for q in (qm, q0, qp))
    pd = {name: sylvester_pivots(q) for name, q in (("minus", qm), ("center", q0), ("plus", qp))}
    ipd = {name: sylvester_pivots(identity_minus(q))
           for name, q in (("minus", qm), ("center", q0), ("plus", qp))}
    assert all(x > 0 for values in [*pd.values(), *ipd.values()] for x in values)

    direct = {name: complete_atoms_direct(q)
              for name, q in (("minus", qm), ("center", q0), ("plus", qp))}
    mobius = {name: complete_atoms_mobius(q)
              for name, q in (("minus", qm), ("center", q0), ("plus", qp))}
    assert direct == mobius
    assert all(p > 0 for values in direct.values() for p in values)
    assert all(sum(values, F(0)) == 1 for values in direct.values())
    assert direct["minus"] == direct["plus"]

    pairs = {}
    for mask, pair in enumerate(zip(direct["center"], direct["plus"])):
        pairs.setdefault(pair, []).append(mask)
    classes = sorted(pairs.items(), key=lambda x: x[1][0])
    assert len(classes) == 8
    assert sum(len(masks) for _, masks in classes) == 64

    hc = entropy_bounds(direct["center"])
    he = entropy_bounds(direct["plus"])
    gap = (he[0] - hc[1], he[1] - hc[0])
    assert gap[0] > F(1, 3_700_000)
    assert gap[1] < F(1, 3_600_000)

    with (HERE / "all_64_atoms.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["mask", "subset", "p_center", "p_plus", "p_minus",
                    "positive", "endpoints_equal"])
        for mask in range(1 << N):
            pc, pp, pm = direct["center"][mask], direct["plus"][mask], direct["minus"][mask]
            w.writerow([mask, subset(mask), fstr(pc), fstr(pp), fstr(pm),
                        pc > 0 and pp > 0 and pm > 0, pp == pm])

    result = {
        "method": "independent exact Q(i) Gaussian elimination; no candidate import",
        "K_real_symmetric": k_symmetric,
        "B_real_skew_symmetric": b_skew and b_zero_diag,
        "Q_hermitian": True,
        "candidate_weyl_argument": {
            "C_eigenvectors_checked": True,
            "K_eigenvalues": [fstr(x) for x in k_eigs],
            "T_transpose_T_equals_3I_minus_J": True,
            "B_operator_norm_squared": fstr(bnorm_sq),
            "norm_squared_strictly_below_K_margin_squared": True,
            "K_commutes_with_B": True,
            "conclusion": "valid"
        },
        "sylvester_leading_minors_Q": {k: [fstr(x) for x in v] for k, v in pd.items()},
        "sylvester_leading_minors_I_minus_Q": {k: [fstr(x) for x in v] for k, v in ipd.items()},
        "strict_positive_contractions": True,
        "direct_equals_mobius_all_64": True,
        "all_atoms_positive": True,
        "sums": {k: fstr(sum(v, F(0))) for k, v in direct.items()},
        "endpoint_laws_equal_eventwise": True,
        "minimum_atoms": {k: fstr(min(v)) for k, v in direct.items()},
        "number_of_classes": len(classes),
        "classes": [
            {"multiplicity": len(masks), "masks": masks,
             "representative_subset": subset(masks[0]),
             "center_probability": fstr(pair[0]), "endpoint_probability": fstr(pair[1])}
            for pair, masks in classes
        ],
        "log_terms": 60,
        "entropy_center_interval": [fstr(x) for x in hc],
        "entropy_endpoint_interval": [fstr(x) for x in he],
        "gap_interval": [fstr(x) for x in gap],
        "gap_decimal": [dstr(x, 40) for x in gap],
        "gap_interval_width_scientific": sci(gap[1] - gap[0], 30),
        "proved": "1/3700000 < H(Q+)-H(K) < 1/3600000",
        "frozen_primary_threshold": "1/1000000",
        "frozen_primary_threshold_met": gap[0] > F(1, 1_000_000),
        "frozen_stretch_threshold_met": gap[0] > F(1, 10_000),
    }
    (HERE / "independent_results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({
        "classes": [len(m) for _, m in classes],
        "class_representatives": [subset(m[0]) for _, m in classes],
        "min_center": fstr(min(direct["center"])),
        "min_endpoint": fstr(min(direct["plus"])),
        "gap_decimal": result["gap_decimal"],
        "gap_width": result["gap_interval_width_scientific"],
        "frozen_primary_threshold_met": result["frozen_primary_threshold_met"],
        "status": "PASS",
    }, indent=2))


if __name__ == "__main__":
    main()
