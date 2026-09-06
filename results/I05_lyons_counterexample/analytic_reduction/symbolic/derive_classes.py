#!/usr/bin/env python3
"""Symbolically derive the eight probability classes for the three-parameter family.

The script is an audit tool, not a numerical enumerator: every atom is first
formed as a symbolic determinant and then grouped by exact polynomial equality.
Several rational substitutions are checked independently by exact determinants.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
from pathlib import Path

import sympy as sp


N = 6
c, d, s = sp.symbols("c d s", real=True)
I = sp.I


def matrices():
    C = sp.Matrix([
        [c + d, c - 2*d, c + d],
        [c - 2*d, c + d, c + d],
        [c + d, c + d, c - 2*d],
    ]) / 3
    S0 = sp.Matrix([[0, 1, -1], [-1, 0, 1], [1, -1, 0]])
    E = sp.eye(3)
    K = E.row_join(C).col_join(C.row_join(E)) / 2
    B = s * S0.row_join(sp.zeros(3)).col_join(sp.zeros(3).row_join(-S0))
    return C, S0, K, B


def atom(M: sp.Matrix, mask: int) -> sp.Expr:
    """Exact-event probability (-1)^|A^c| det(M-I_{A^c})."""
    complement = [j for j in range(N) if not (mask >> j) & 1]
    X = M.copy()
    for j in complement:
        X[j, j] -= 1
    return sp.factor((-1)**len(complement) * X.det(method="domain-ge"))


def subset(mask: int) -> str:
    xs = [str(j + 1) for j in range(N) if (mask >> j) & 1]
    return "{" + ",".join(xs) + "}"


def class_key(p: sp.Expr, q: sp.Expr):
    return (sp.Poly(sp.expand(p), c, d, s), sp.Poly(sp.expand(q), c, d, s))


def parity(p):
    return -1 if sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p))) % 2 else 1


def det_fraction_complex(a):
    """Independent Leibniz determinant over Python complex pairs of Fraction."""
    total = (Fraction(0), Fraction(0))
    for perm in permutations(range(len(a))):
        x = (Fraction(parity(perm)), Fraction(0))
        for i, j in enumerate(perm):
            y = a[i][j]
            x = (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
        total = (total[0] + x[0], total[1] + x[1])
    return total


def independent_atom(CV, DV, SV, mask, endpoint):
    cv, dv, sv = map(Fraction, (CV, DV, SV))
    CC = [[(cv+dv)/3, (cv-2*dv)/3, (cv+dv)/3],
          [(cv-2*dv)/3, (cv+dv)/3, (cv+dv)/3],
          [(cv+dv)/3, (cv+dv)/3, (cv-2*dv)/3]]
    S0 = [[0, 1, -1], [-1, 0, 1], [1, -1, 0]]
    M = []
    for i in range(N):
        row = []
        for j in range(N):
            re = Fraction(i == j, 2)
            if (i < 3) != (j < 3):
                re = CC[i % 3][j % 3] / 2
            im = Fraction(0)
            if endpoint and i//3 == j//3:
                im = sv * S0[i % 3][j % 3] * (1 if i < 3 else -1)
            if not (mask >> i) & 1 and i == j:
                re -= 1
            row.append((re, im))
        M.append(row)
    val = det_fraction_complex(M)
    sign = -1 if (N - mask.bit_count()) % 2 else 1
    return sign*val[0], sign*val[1]


def main():
    C, S0, K, B = matrices()
    assert sp.simplify(C*S0 + S0*C) == sp.zeros(3)
    assert sp.simplify(K*B - B*K) == sp.zeros(6)

    center = []
    endpoint = []
    for mask in range(1 << N):
        center.append(atom(K, mask))
        endpoint.append(atom(K + I*B, mask))
    assert all(sp.factor(q - q.subs(s, -s)) == 0 for q in endpoint)

    buckets = {}
    for mask, (p, q) in enumerate(zip(center, endpoint)):
        buckets.setdefault(class_key(p, q), []).append(mask)
    assert len(buckets) == 8
    classes = sorted(buckets.items(), key=lambda item: item[1][0])
    assert [len(masks) for _, masks in classes] == [2, 12, 12, 2, 12, 6, 12, 6]
    assert sp.factor(sum(center) - 1) == 0
    assert sp.factor(sum(endpoint) - 1) == 0

    z = s**2
    claimed_P = [
        9*(1-c**2)*(1-d**2)**2,
        3*(1-d**2)*(3+d**2-c**2*(1+3*d**2)),
        3*(1+d**2)*(3-d**2+c**2*(1-3*d**2)),
        9*(1+c**2)*(1+d**2)**2,
        (d**2-2*c*d+3)**2-(c*(3*d**2+1)-2*d)**2,
        (1-d**2)*((3+3*c*d)**2-(c+d)**2),
        9*c**2*d**4-2*c**2*d**2+c**2+8*c*d**3+8*c*d+d**4-2*d**2+9,
        (1+d**2)*(9*c**2*d**2+c**2-16*c*d+d**2+9),
    ]
    claimed_D = [
        -216*z*(1-c**2)*(1+d**2-6*z),
        72*z*(d**2-1-c**2*(3*d**2+1)+6*z*(3*c**2-1)),
        72*z*(3*c**2*d**2+c**2+d**2-1-6*z*(3*c**2+1)),
        -216*z*(1+c**2)*(1+d**2-6*z),
        24*z*(9*c**2*d**2-3*c**2-4*c*d-d**2+3+6*z*(1-9*c**2)),
        24*z*(9*c**2*d**2-3*c**2+8*c*d-d**2+3+6*z*(1-9*c**2)),
        24*z*(-9*c**2*d**2+3*c**2-4*c*d-d**2+3+6*z*(9*c**2+1)),
        24*z*(-9*c**2*d**2+3*c**2+8*c*d-d**2+3+6*z*(9*c**2+1)),
    ]
    for j, (((ppoly, qpoly), _), Pj, Dj) in enumerate(zip(classes, claimed_P, claimed_D), 1):
        assert sp.factor(576*ppoly.as_expr()-Pj) == 0, ("P", j)
        assert sp.factor(576*(qpoly.as_expr()-ppoly.as_expr())-Dj) == 0, ("D", j)
    multiplicities = [len(masks) for _, masks in classes]
    assert sp.factor(sum(m*Pj for m, Pj in zip(multiplicities, claimed_P))-576) == 0
    assert sp.factor(sum(m*Dj for m, Dj in zip(multiplicities, claimed_D))) == 0

    # Coordinate symmetries of the symbolic pair.  A permutation is retained
    # when it fixes K and sends B to either B or -B; the latter is harmless
    # because every determinant is even in s (equivalently Q_+ and Q_- have
    # the same exact-event law).
    coordinate_group = []
    for perm in permutations(range(N)):
        kp = K.extract(perm, perm)
        bp = B.extract(perm, perm)
        if kp == K and (bp == B or bp == -B):
            coordinate_group.append((perm, 1 if bp == B else -1))
    assert len(coordinate_group) == 12
    rho = (1, 2, 0, 5, 3, 4)
    tau = (1, 0, 2, 4, 3, 5)
    sigma = (3, 4, 5, 0, 1, 2)
    generated = {tuple(range(N))}
    changed = True
    while changed:
        changed = False
        for p in list(generated):
            for g in (rho, tau, sigma):
                h = tuple(p[g[i]] for i in range(N))
                if h not in generated:
                    generated.add(h)
                    changed = True
    assert generated == {perm for perm, _ in coordinate_group}

    def permute_mask(mask, perm):
        return sum(1 << i for i in range(N) if (mask >> perm[i]) & 1)

    orbit_partition = []
    unseen = set(range(64))
    while unseen:
        seed = min(unseen)
        orbit = {permute_mask(seed, perm) for perm, _ in coordinate_group}
        orbit |= {63 ^ m for m in list(orbit)}
        orbit_partition.append(sorted(orbit))
        unseen -= orbit
    assert orbit_partition == [masks for _, masks in classes]

    # Independent exact checks at several interior rational points.  These do
    # not prove the formulas; they are intended to catch transcription bugs.
    points = [
        (Fraction(4, 5), Fraction(9, 10), Fraction(1, 100)),
        (Fraction(7, 10), Fraction(4, 5), Fraction(1, 50)),
        (Fraction(1231, 1250), Fraction(497, 500), Fraction(11, 10000)),
    ]
    for cv, dv, sv in points:
        for mask in range(64):
            for ep, expr in [(False, center[mask]), (True, endpoint[mask])]:
                got = independent_atom(cv, dv, sv, mask, ep)
                expected = Fraction(expr.subs({c: sp.Rational(cv.numerator, cv.denominator),
                                               d: sp.Rational(dv.numerator, dv.denominator),
                                               s: sp.Rational(sv.numerator, sv.denominator)}))
                assert got == (expected, Fraction(0)), (cv, dv, sv, mask, ep, got, expected)

    lines = []
    lines.append("# Machine-derived class table\n")
    lines.append("All expressions below are factored exact SymPy polynomials.\n")
    lines.append(f"The symbolic coordinate group has order {len(coordinate_group)}.\n")
    for j, ((ppoly, qpoly), masks) in enumerate(classes, 1):
        p = sp.factor(ppoly.as_expr())
        q = sp.factor(qpoly.as_expr())
        lines += [
            f"## Class {j}\n",
            f"- representative: `{subset(masks[0])}`\n",
            f"- multiplicity: `{len(masks)}`\n",
            f"- masks: `{masks}`\n",
            f"- center: `${sp.latex(p)}$`\n",
            f"- endpoint: `${sp.latex(q)}$`\n",
            f"- displacement: `${sp.latex(sp.factor(q-p))}$`\n",
        ]
    out = Path(__file__).with_name("machine_classes.md")
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"proved symbolically: {len(classes)} classes")
    print("multiplicities:", [len(masks) for _, masks in classes])
    print("representatives:", [subset(masks[0]) for _, masks in classes])
    print("proved: center and endpoint atoms each sum to 1 identically")
    print("proved: human-readable P_j,D_j formulas and scaled normalizations")
    print("proved: coordinate symmetries plus complementation have exactly these orbits")
    print("coordinate symmetries (1-based; sign is B -> sign B):")
    for perm, sign in coordinate_group:
        print(" ", tuple(i+1 for i in perm), sign)
    print("cross-checked at", len(points), "exact rational points")
    print("wrote", out)


if __name__ == "__main__":
    main()
