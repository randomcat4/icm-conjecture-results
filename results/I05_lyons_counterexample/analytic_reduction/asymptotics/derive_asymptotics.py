#!/usr/bin/env python3
"""Exact symbolic boundary expansion for the split (c,d,s) family.

All determinants and series are computed by SymPy.  The script reconstructs
the 64 complete-event probabilities, groups them by exact polynomial equality,
and expands them under c=1-A*tau, d=1-B*tau, s=R*tau.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


c, d, s = sp.symbols("c d s", real=True)
A, B, R, tau = sp.symbols("A B R tau", positive=True)
LT = sp.symbols("LT", real=True)
I = sp.I
N = 6


def matrices():
    C = sp.Matrix([
        [c + d, c - 2*d, c + d],
        [c - 2*d, c + d, c + d],
        [c + d, c + d, c - 2*d],
    ]) / 3
    S0 = sp.Matrix([[0, 1, -1], [-1, 0, 1], [1, -1, 0]])
    E = sp.eye(3)
    K = E.row_join(C).col_join(C.row_join(E)) / 2
    D = S0.row_join(sp.zeros(3)).col_join(sp.zeros(3).row_join(-S0))
    return C, S0, K, D


def atom(M: sp.Matrix, mask: int) -> sp.Expr:
    complement = [j for j in range(N) if not (mask >> j) & 1]
    X = M.copy()
    for j in complement:
        X[j, j] -= 1
    return sp.factor((-1) ** len(complement) * X.det(method="domain-ge"))


def leading(expr: sp.Expr):
    poly = sp.Poly(sp.expand(expr), tau)
    for (power,), coeff in reversed(poly.terms()):
        if coeff != 0:
            return power, sp.factor(coeff)
    return None, sp.Integer(0)


def entropy_series(expr: sp.Expr, degree: int) -> sp.Expr:
    """Series for -expr*log(expr), with LT standing for log(tau)."""
    order, _ = leading(expr)
    regular = sp.cancel(expr / tau**order)
    log_regular = sp.series(sp.log(regular), tau, 0, degree + 1).removeO()
    ans = -expr * (order * LT + log_regular)
    return sp.series(ans, tau, 0, degree + 1).removeO()


def main():
    C, S0, K, D = matrices()
    assert sp.simplify(C*S0 + S0*C) == sp.zeros(3)
    assert sp.simplify(K*D - D*K) == sp.zeros(6)

    center = [atom(K, mask) for mask in range(1 << N)]
    endpoint = [atom(K + I*s*D, mask) for mask in range(1 << N)]
    groups = {}
    for mask, (p, q) in enumerate(zip(center, endpoint)):
        key = (sp.Poly(sp.expand(p), c, d, s), sp.Poly(sp.expand(q), c, d, s))
        groups.setdefault(key, []).append(mask)
    classes = sorted(groups.items(), key=lambda item: item[1][0])
    assert len(classes) == 8
    assert [len(masks) for _, masks in classes] == [2, 12, 12, 2, 12, 6, 12, 6]
    assert sp.factor(sum(center) - 1) == 0
    assert sp.factor(sum(endpoint) - 1) == 0

    scaling = {c: 1-A*tau, d: 1-B*tau, s: R*tau}
    data = []
    scaled_rows = []
    for j, ((ppoly, qpoly), masks) in enumerate(classes, 1):
        p = sp.factor(ppoly.as_expr())
        q = sp.factor(qpoly.as_expr())
        ps = sp.factor(p.subs(scaling))
        qs = sp.factor(q.subs(scaling))
        dp = sp.factor(q-p)
        dps = sp.factor(dp.subs(scaling))
        p_order, p_lead = leading(ps)
        q_order, q_lead = leading(qs)
        dp_order, dp_lead = leading(dps)
        data.append({
            "class": j,
            "representative_mask": masks[0],
            "multiplicity": len(masks),
            "p": str(p),
            "q": str(q),
            "delta_p": str(dp),
            "p_order": p_order,
            "p_lead": str(p_lead),
            "q_order": q_order,
            "q_lead": str(q_lead),
            "delta_order": dp_order,
            "delta_lead": str(dp_lead),
            "p_series_4": str(sp.series(ps, tau, 0, 5).removeO()),
            "q_series_4": str(sp.series(qs, tau, 0, 5).removeO()),
        })
        scaled_rows.append((len(masks), ps, qs))

    entropy_delta_series = sp.Integer(0)
    for multiplicity, ps, qs in scaled_rows:
        entropy_delta_series += multiplicity * (
            entropy_series(qs, 3) - entropy_series(ps, 3)
        )
    entropy_delta_series = sp.collect(sp.expand(entropy_delta_series), [tau, LT])

    x, r = sp.symbols("x r", positive=True)
    X_general = B*(2*A+B)
    Y_general = X_general - 12*R**2
    general_boundary_coefficient = sp.Rational(1, 2) * (
        X_general*sp.log(X_general) - Y_general*sp.log(Y_general)
    ) + R**2 * (
        sp.log(sp.Rational(3**6, 2**4)
               / ((A+2*B)**6 * (3*A+2*B)**2 * B**4)) - 6
    )
    X = 2*x + 1
    Y = X - 12*r**2
    boundary_coefficient = sp.Rational(1, 2) * (X*sp.log(X) - Y*sp.log(Y))
    boundary_coefficient += r**2 * (
        sp.log(sp.Rational(3**6, 2**4) / ((x+2)**6 * (3*x+2)**2)) - 6
    )
    small_r_quadratic = sp.log(
        sp.Rational(3**6, 2**4) * (2*x+1)**6
        / ((x+2)**6 * (3*x+2)**2)
    )
    small_r_series = sp.series(boundary_coefficient, r, 0, 7).removeO()
    assert sp.simplify(sp.expand_log(
        small_r_series.coeff(r, 2) - small_r_quadratic, force=True
    )) == 0
    assert sp.simplify(small_r_series.coeff(r, 4) + 36/(2*x+1)) == 0
    expected_g2_derivative = (
        12*(-x**2 + 2*x + 2) / ((2*x+1)*(x+2)*(3*x+2))
    )
    assert sp.simplify(sp.diff(small_r_quadratic, x)
                       - expected_g2_derivative) == 0

    tau2 = sp.expand(entropy_delta_series).coeff(tau, 2)
    tau3 = sp.expand(entropy_delta_series).coeff(tau, 3)
    tau3_log = sp.expand(tau3).coeff(LT)
    tau3_constant = sp.expand(tau3).subs(LT, 0)
    tau2_difference = sp.simplify(sp.expand_log(
        tau2 - general_boundary_coefficient, force=True
    ))
    assert tau2_difference == 0
    assert sp.simplify(tau3_log - 6*R**2*(A+2*B)) == 0

    # Numerical check at the currently observed finite maximizer, normalized
    # with tau=beta and B=1.
    cv = sp.Rational("0.9849709285")
    dv = sp.Rational("0.9939746109")
    sv = sp.Rational("0.0011147170")
    beta_v = 1-dv
    x_v = (1-cv)/beta_v
    r_v = sv/beta_v
    current_subs = {A: x_v, B: 1, R: r_v}
    tau2_direct_v = sp.N(tau2.subs(current_subs), 40)
    F_v = sp.N(boundary_coefficient.subs({x: x_v, r: r_v}), 40)
    J_v = sp.N(tau3_constant.subs(current_subs), 40)
    log_v = sp.N(tau3_log.subs(current_subs), 40)
    approximation_2 = sp.N(beta_v**2 * F_v, 40)
    approximation_3 = sp.N(
        beta_v**2 * F_v + beta_v**3 * (log_v*sp.log(beta_v) + J_v), 40
    )
    exact_current = sp.N(sum(
        multiplicity * (-q.subs({c: cv, d: dv, s: sv})
                        * sp.log(q.subs({c: cv, d: dv, s: sv}))
                        + p.subs({c: cv, d: dv, s: sv})
                        * sp.log(p.subs({c: cv, d: dv, s: sv})))
        for multiplicity, p, q in [
            (len(masks), ppoly.as_expr(), qpoly.as_expr())
            for (ppoly, qpoly), masks in classes
        ]
    ), 40)

    out = Path(__file__).with_name("asymptotic_data.json")
    payload = {
        "classes": data,
        "delta_series_through_tau3": str(entropy_delta_series),
        "boundary_coefficient_general": str(general_boundary_coefficient),
        "boundary_coefficient_B_equals_1": str(boundary_coefficient),
        "small_r_quadratic_coefficient": str(small_r_quadratic),
        "small_r_series_through_r6": str(small_r_series),
        "tau3_log_coefficient": str(sp.factor(tau3_log)),
        "current_check": {
            "alpha_over_beta": str(x_v),
            "s_over_beta": str(r_v),
            "F": str(F_v),
            "tau2_direct": str(tau2_direct_v),
            "tau2_F": str(approximation_2),
            "tau3_constant_J": str(J_v),
            "tau3_log_coefficient": str(log_v),
            "through_tau3": str(approximation_3),
            "exact_delta": str(exact_current),
        },
    }
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print("classes:", len(data))
    for row in data:
        print(
            row["class"], row["multiplicity"],
            "p ~ tau^", row["p_order"], "*", row["p_lead"],
            "; q ~ tau^", row["q_order"], "*", row["q_lead"],
            "; q-p ~ tau^", row["delta_order"], "*", row["delta_lead"],
        )
    print("boundary coefficient for B=1:")
    print(boundary_coefficient)
    print("small-r quadratic coefficient:")
    print(small_r_quadratic)
    print("tau^3 log(tau) coefficient:", sp.factor(tau3_log))
    print("current normalized ratios x,r:", sp.N(x_v, 16), sp.N(r_v, 16))
    print("current tau^2 F:", approximation_2)
    print("current direct tau^2 coefficient:", tau2_direct_v)
    print("proved: symbolic tau^2 coefficient equals the displayed boundary formula")
    print("full tau^3 coefficient is stored in asymptotic_data.json")
    print("current expansion through tau^3:", approximation_3)
    print("current exact Delta:", exact_current)
    print("wrote", out)


if __name__ == "__main__":
    main()
