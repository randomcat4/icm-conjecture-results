#!/usr/bin/env python3
"""Independent exact audit for paper commit 4a58b2b9c6a9a20c18f3f330081b5d3bc8ac3f3e.

The calculation uses only the Python standard library.  The matrices and the
short rational bounds are transcribed from the manuscript.  No author or prior
verifier program is imported or executed.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import re
from fractions import Fraction as F
from pathlib import Path

G = tuple[F, F]  # a + b i
ZERO: G = (F(0), F(0))
ONE: G = (F(1), F(0))


def gadd(x: G, y: G) -> G:
    return (x[0] + y[0], x[1] + y[1])


def gmul(x: G, y: G) -> G:
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def gscale(x: G, q: F) -> G:
    return (q * x[0], q * x[1])


def parity(perm: tuple[int, ...]) -> int:
    inv = sum(perm[i] > perm[j] for i in range(len(perm)) for j in range(i + 1, len(perm)))
    return -1 if inv % 2 else 1


def det_g(a: list[list[G]]) -> G:
    n = len(a)
    if n == 0:
        return ONE
    total = ZERO
    for p in itertools.permutations(range(n)):
        term = ONE
        for i, j in enumerate(p):
            term = gmul(term, a[i][j])
        total = gadd(total, gscale(term, F(parity(p))))
    return total


def padd(a: list[G], b: list[G]) -> list[G]:
    out = [ZERO] * max(len(a), len(b))
    for i in range(len(out)):
        out[i] = gadd(a[i] if i < len(a) else ZERO, b[i] if i < len(b) else ZERO)
    return out


def pmul(a: list[G], b: list[G]) -> list[G]:
    out = [ZERO] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = gadd(out[i + j], gmul(x, y))
    return out


def det_poly(a: list[list[list[G]]]) -> list[G]:
    n = len(a)
    total = [ZERO]
    for p in itertools.permutations(range(n)):
        term = [ONE]
        for i, j in enumerate(p):
            term = pmul(term, a[i][j])
        total = padd(total, [gscale(x, F(parity(p))) for x in term])
    total += [ZERO] * (n + 1 - len(total))
    return total[: n + 1]


def inverse_q(a: list[list[F]]) -> list[list[F]]:
    n = len(a)
    aug = [row[:] + [F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for col in range(n):
        pivot = next(r for r in range(col, n) if aug[r][col])
        aug[col], aug[pivot] = aug[pivot], aug[col]
        q = aug[col][col]
        aug[col] = [x / q for x in aug[col]]
        for r in range(n):
            if r != col and aug[r][col]:
                q = aug[r][col]
                aug[r] = [x - q * y for x, y in zip(aug[r], aug[col])]
    return [row[n:] for row in aug]


def matmul_q(a: list[list[F]], b: list[list[F]]) -> list[list[F]]:
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), F(0)) for j in range(len(b[0]))] for i in range(len(a))]


def submatrix(a, idx: list[int]):
    return [[a[i][j] for j in idx] for i in idx]


def endpoint_matrix(k: list[list[F]], b: list[list[F]], t: F) -> list[list[G]]:
    return [[(k[i][j], t * b[i][j]) for j in range(5)] for i in range(5)]


def identity_minus(a: list[list[G]]) -> list[list[G]]:
    return [[(F(i == j) - a[i][j][0], -a[i][j][1]) for j in range(5)] for i in range(5)]


def leading_minors(a: list[list[G]]) -> list[F]:
    out = []
    for r in range(1, 6):
        d = det_g([row[:r] for row in a[:r]])
        assert d[1] == 0
        out.append(d[0])
    return out


def event_polynomial(k: list[list[F]], b: list[list[F]], mask: str) -> list[F]:
    # Frozen tables use ordinary integer-mask order: the rightmost bit is coordinate 1.
    s = {i for i, bit in enumerate(reversed(mask)) if bit == "1"}
    c = set(range(5)) - s
    matrix: list[list[list[G]]] = []
    for i in range(5):
        row = []
        for j in range(5):
            constant = k[i][j] - F(i == j and i in c)
            row.append([(constant, F(0)), (F(0), b[i][j])])
        matrix.append(row)
    factor = F(-1 if len(c) % 2 else 1)
    poly = [gscale(x, factor) for x in det_poly(matrix)]
    assert all(im == 0 for _, im in poly)
    return [re for re, _ in poly]


def event_direct(q: list[list[G]], mask: str) -> F:
    s = {i for i, bit in enumerate(reversed(mask)) if bit == "1"}
    c = set(range(5)) - s
    m = [[gadd(q[i][j], (F(-(i == j and i in c)), F(0))) for j in range(5)] for i in range(5)]
    d = gscale(det_g(m), F(-1 if len(c) % 2 else 1))
    assert d[1] == 0
    return d[0]


def event_mobius(q: list[list[G]], mask: str) -> F:
    s = {i for i, bit in enumerate(reversed(mask)) if bit == "1"}
    c = [i for i in range(5) if i not in s]
    total = ZERO
    for r in range(len(c) + 1):
        for u in itertools.combinations(c, r):
            idx = sorted(s | set(u))
            total = gadd(total, gscale(det_g(submatrix(q, idx)), F(-1 if r % 2 else 1)))
    assert total[1] == 0
    return total[0]


def log_interval(x: F, terms: int = 80) -> tuple[F, F]:
    assert x > 0
    y, power = x, 0
    while y < 1:
        y *= 2
        power -= 1
    while y >= 2:
        y /= 2
        power += 1
    z = (y - 1) / (y + 1)
    low_y = 2 * sum((z ** (2 * m + 1) / (2 * m + 1) for m in range(terms)), F(0))
    high_y = low_y + 2 * z ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z * z))
    z2 = F(1, 3)
    low_2 = 2 * sum((z2 ** (2 * m + 1) / (2 * m + 1) for m in range(terms)), F(0))
    high_2 = low_2 + 2 * z2 ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z2 * z2))
    if power >= 0:
        return low_y + power * low_2, high_y + power * high_2
    return low_y + power * high_2, high_y + power * low_2


def weighted_log_interval(terms: list[tuple[F, F]]) -> tuple[F, F]:
    lo = hi = F(0)
    for coefficient, argument in terms:
        l, u = log_interval(argument)
        if coefficient >= 0:
            lo += coefficient * l
            hi += coefficient * u
        else:
            lo += coefficient * u
            hi += coefficient * l
    return lo, hi


def parse_fraction(text: str) -> F:
    return F(text)


def read_tsv(path: Path) -> dict[str, dict[str, F]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = {}
        for row in csv.DictReader(handle, delimiter="\t"):
            mask = row.pop("mask")
            parsed = {}
            for key in list(row):
                if key.endswith("_num"):
                    stem = key[:-4]
                    den_key = "q_den" if stem == "q_equals_minus_p2" else f"{stem}_den"
                    parsed[stem] = F(int(row[key]), int(row[den_key]))
            rows[mask] = parsed
        return rows


def ratio_strings(values: list[F]) -> list[str]:
    return [f"{x.numerator}/{x.denominator}" for x in values]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paper_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    paper = args.paper_dir

    kn = [
        [2719, -3449, 1009, -818, -2490],
        [-3449, 5840, -2889, -1513, 1322],
        [1009, -2889, 2148, 2505, 1101],
        [-818, -1513, 2505, 4698, 3958],
        [-2490, 1322, 1101, 3958, 4596],
    ]
    bn = [
        [0, -3, -13, -16, 4],
        [3, 0, -4, 15, 17],
        [13, 4, 0, 6, -15],
        [16, -15, -6, 0, 1],
        [-4, -17, 15, -1, 0],
    ]
    k = [[F(x, 10000) for x in row] for row in kn]
    b = [[F(x, 50) for x in row] for row in bn]
    h = F(1, 100000)
    assert all(k[i][j] == k[j][i] for i in range(5) for j in range(5))
    assert all(b[i][j] == -b[j][i] for i in range(5) for j in range(5))

    q0 = endpoint_matrix(k, b, F(0))
    qp = endpoint_matrix(k, b, h)
    qm = endpoint_matrix(k, b, -h)
    feasibility = {
        "K": leading_minors(q0),
        "I-K": leading_minors(identity_minus(q0)),
        "K+hA": leading_minors(qp),
        "K-hA": leading_minors(qm),
        "I-K-hA": leading_minors(identity_minus(qp)),
        "I-K+hA": leading_minors(identity_minus(qm)),
    }
    assert all(x > 0 for values in feasibility.values() for x in values)
    assert feasibility["K+hA"] == feasibility["K-hA"]
    assert feasibility["I-K-hA"] == feasibility["I-K+hA"]
    expected = {
        "K": [F(2719, 10**4), F(3983359, 10**8), F(24692191, 10**12), F(3908717, 5 * 10**14), F(118205061, 10**20)],
        "I-K": [F(7281, 10**4), F(18393359, 10**8), F(59312197809, 10**12), F(35054548217, 5 * 10**14), F(1709337154939, 10**20)],
        "K+hA": [F(2719, 10**4), F(995839749991, 25 * 10**12), F(123460921839, 5 * 10**15), F(4883537625128308081, 625 * 10**24), F(7353066474454085031, 625 * 10**28)],
        "I-K-hA": [F(7281, 10**4), F(4598339749991, 25 * 10**12), F(296560989039361, 5 * 10**15), F(43818177021217628308081, 625 * 10**24), F(106676909868459401834969, 625 * 10**28)],
    }
    for name, values in expected.items():
        assert feasibility[name] == values

    b4 = [[(b[i][j], F(0)) for j in range(4)] for i in range(4)]
    b5 = [[(b[i][j], F(0)) for j in range(5)] for i in range(5)]
    assert det_g(b4) == (F(58081, 6250000), F(0))
    assert det_g(b5) == ZERO

    center_table = read_tsv(paper / "certificates" / "event_certificate.tsv")
    chord_table = read_tsv(paper / "certificates" / "chord_event_certificate.tsv")
    assert len(center_table) == len(chord_table) == 32
    events = {}
    summed_poly = [F(0)] * 6
    for number in range(32):
        mask = format(number, "05b")
        poly = event_polynomial(k, b, mask)
        assert all(poly[i] == 0 for i in (1, 3, 5))
        summed_poly = [x + y for x, y in zip(summed_poly, poly)]
        p0, p1, p2 = poly[0], poly[1], 2 * poly[2]
        ph_poly = sum((coefficient * h**degree for degree, coefficient in enumerate(poly)), F(0))
        pm_poly = sum((coefficient * (-h) ** degree for degree, coefficient in enumerate(poly)), F(0))
        direct = [event_direct(q, mask) for q in (q0, qp, qm)]
        mobius = [event_mobius(q, mask) for q in (q0, qp, qm)]
        assert direct == mobius
        assert direct == [p0, ph_poly, pm_poly]
        assert ph_poly == pm_poly
        assert p0 > 0 and ph_poly > 0
        expected_center_row = {"p": p0, "p1": p1, "p2": p2, "q_equals_minus_p2": -p2}
        assert center_table[mask] == expected_center_row, (mask, center_table[mask], expected_center_row)
        assert chord_table[mask] == {"p0": p0, "ph": ph_poly}

        c = {i for i, bit in enumerate(reversed(mask)) if bit == "0"}
        m = [[k[i][j] - F(i == j and i in c) for j in range(5)] for i in range(5)]
        invm = inverse_q(m)
        im_b = matmul_q(invm, b)
        trace_square = sum(matmul_q(im_b, im_b)[i][i] for i in range(5))
        assert p2 == p0 * trace_square
        events[mask] = (p0, ph_poly, p2)

    assert summed_poly == [F(1), F(0), F(0), F(0), F(0), F(0)]
    assert sum((v[0] for v in events.values()), F(0)) == 1
    assert sum((v[1] for v in events.values()), F(0)) == 1
    assert sum((v[2] for v in events.values()), F(0)) == 0

    table_text = (paper / "certificates" / "event_table.tex").read_text(encoding="utf-8")
    table_rows = re.findall(r"\\texttt\{([01]{5})\}\s*&\s*(-?\d+)\s*&\s*(-?\d+)\s*&\s*(-?\d+)\\\\", table_text)
    assert len(table_rows) == 32
    for mask, center_num, endpoint_num, second_num in table_rows:
        p0, ph, p2 = events[mask]
        assert (int(center_num), int(endpoint_num), int(second_num)) == (
            p0 * 10**20,
            ph * (625 * 10**28),
            p2 * 312500000000000,
        )

    hessian_terms = [(-p2, p0) for p0, _, p2 in events.values()]
    hess_lo, hess_hi = weighted_log_interval(hessian_terms)
    hess_short_lo = F(5205948204140858303, 10**21)
    hess_short_hi = F(5205948204140858304, 10**21)
    assert hess_short_lo < hess_lo < hess_hi < hess_short_hi
    assert F(1, 200) < hess_short_lo and hess_short_hi < F(1, 190)
    hess_lower_margin = hess_short_lo - F(1, 200)
    hess_upper_margin = F(1, 190) - hess_short_hi
    assert hess_lower_margin == F(205948204140858303, 10**21)
    assert hess_upper_margin == F(16984126895682691, 296875000000000000000)

    chord_terms = []
    for p0, ph, _ in events.values():
        chord_terms.extend([(-ph, ph), (p0, p0)])
    gap_lo, gap_hi = weighted_log_interval(chord_terms)
    gap_short_lo = F(233555253804221762, 10**30)
    gap_short_hi = F(233555253804221763, 10**30)
    assert gap_short_lo < gap_lo < gap_hi < gap_short_hi
    assert gap_short_lo > F(1, 10**13)
    gap_margin = gap_short_lo - F(1, 10**13)
    assert gap_margin == F(66777626902110881, 500000000000000000000000000000)

    summary_json = json.loads((paper / "certificates" / "certificate_summary.json").read_text(encoding="utf-8"))
    chord_json = json.loads((paper / "certificates" / "chord_certificate_summary.json").read_text(encoding="utf-8"))
    assert summary_json["K_numerators"] == kn and summary_json["B_numerators"] == bn
    assert summary_json["leading_principal_minors_K"] == ratio_strings(feasibility["K"])
    assert summary_json["leading_principal_minors_I_minus_K"] == ratio_strings(feasibility["I-K"])
    assert chord_json["leading_principal_minors_K_plus_hA"] == ratio_strings(feasibility["K+hA"])
    assert chord_json["leading_principal_minors_I_minus_K_minus_hA"] == ratio_strings(feasibility["I-K-hA"])

    receipt = {
        "paper_commit": "4a58b2b9c6a9a20c18f3f330081b5d3bc8ac3f3e",
        "method": "stdlib Fraction; Leibniz determinants over Q(i)[t]; independent Mobius cross-check; 80-term directed atanh intervals",
        "rank_A": 4,
        "leading_principal_minors": {name: ratio_strings(values) for name, values in feasibility.items()},
        "event_polynomials": 32,
        "odd_coefficients_all_zero": True,
        "mobius_matches_determinant_at_center_and_both_endpoints": True,
        "center_positive_and_normalized": True,
        "both_endpoints_positive_normalized_and_termwise_equal": True,
        "event_table_tsv_and_tex_rows_matching": 32,
        "sum_p_second": "0/1",
        "log_series_terms": 80,
        "hessian_verified_enclosure": [str(hess_short_lo), str(hess_short_hi)],
        "hessian_bounds": ["1/200", "1/190"],
        "hessian_printed_margins": [str(hess_lower_margin), str(hess_upper_margin)],
        "chord_gap_verified_enclosure": [str(gap_short_lo), str(gap_short_hi)],
        "chord_gap_lower_bound": "1/10000000000000",
        "chord_gap_printed_margin": str(gap_margin),
        "all_assertions_passed": True,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
