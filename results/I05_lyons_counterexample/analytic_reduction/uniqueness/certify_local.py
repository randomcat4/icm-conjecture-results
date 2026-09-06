#!/usr/bin/env python3
"""Rigorous local stationary-point certificate using Arb ball arithmetic.

Run with
    uv run --offline --with sympy --with python-flint python certify_local.py

The script implements a Krawczyk inclusion for grad(Delta) and a Sylvester
criterion for negative definiteness of Hess(Delta) on the whole input box.
All centers and preconditioner entries are exact dyadic rationals.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp
from flint import arb, ctx, fmpz

HERE = Path(__file__).resolve().parent
ctx.prec = 256

c, d, s = sp.symbols("c d s", real=True)
VARS = (c, d, s)
rows = json.loads((HERE / "formulae.json").read_text(encoding="utf-8"))
delta = sp.Integer(0)
probabilities = []
weighted_polynomials = []
for row in rows:
    loc = {"c": c, "d": d, "s": s}
    p = sp.sympify(row["p"], locals=loc)
    q = sp.sympify(row["q"], locals=loc)
    m = int(row["multiplicity"])
    probabilities.extend((p, q))
    weighted_polynomials.extend(((m, p), (-m, q)))
    delta += m * (-q * sp.log(q) + p * sp.log(p))

ARB_MODULE = {"log": lambda x: x.log()}
def compile_expr(expr):
    return sp.lambdify(VARS, expr, modules=[ARB_MODULE], cse=True)

prob_f = [compile_expr(x) for x in probabilities]
terms = []
for weight, poly in weighted_polynomials:
    r0 = compile_expr(poly)
    r1 = [compile_expr(sp.diff(poly, VARS[i])) for i in range(3)]
    r2 = [[compile_expr(sp.diff(poly, VARS[i], VARS[j]))
           for j in range(3)] for i in range(3)]
    r3 = [[[compile_expr(sp.diff(poly, VARS[i], VARS[j], VARS[k]))
            for k in range(3)] for j in range(3)] for i in range(3)]
    terms.append((weight, r0, r1, r2, r3))

def eval_delta(*x):
    ans = arb(0)
    for w, r0, _, _, _ in terms:
        r = r0(*x)
        ans += w*r*r.log()
    return ans

def eval_gradient(*x):
    ans = [arb(0) for _ in range(3)]
    for w, r0, r1, _, _ in terms:
        r = r0(*x); factor = w*(r.log()+1)
        for i in range(3): ans[i] += factor*r1[i](*x)
    return ans

def eval_hessian(*x):
    ans = [[arb(0) for _ in range(3)] for _ in range(3)]
    for w, r0, r1, r2, _ in terms:
        r = r0(*x); lr = r.log()+1
        dr = [r1[i](*x) for i in range(3)]
        for i in range(3):
          for j in range(3):
            ans[i][j] += w*(lr*r2[i][j](*x) + dr[i]*dr[j]/r)
    return ans

def eval_third(*x, probabilities_override=None):
    ans = [[[arb(0) for _ in range(3)] for _ in range(3)] for _ in range(3)]
    for term_index, (w, r0, r1, r2, r3) in enumerate(terms):
        r = (r0(*x) if probabilities_override is None
             else probabilities_override[term_index])
        lr = r.log()+1
        dr = [r1[i](*x) for i in range(3)]
        d2 = [[r2[i][j](*x) for j in range(3)] for i in range(3)]
        for i in range(3):
          for j in range(3):
            for k in range(3):
              ans[i][j][k] += w*(lr*r3[i][j][k](*x)
                + (d2[i][j]*dr[k]+d2[i][k]*dr[j]+d2[j][k]*dr[i])/r
                - dr[i]*dr[j]*dr[k]/(r*r))
    return ans

DEN_BITS = 160
ROOT_NUM = [
    1439536244858795825703648138477772648244438200086,
    1452695378559335846296495333266843919245350418322,
    1629166850848927691821404480252899240504455452,
]
INV_HESS_NUM = [
    [-434874911733729747703024153429984373737076435497182,
     -124469171711489335261734429677849203969046201753677,
       13545991794332354879423005061808997597528705106382],
    [-124469171711489335261734429677849203969046201753677,
      -68755059027591314344113359510915350893849740028392,
        5817103300398370867166561395481288721898524257312],
    [  13545991794332354879423005061808997597528705106382,
        5817103300398370867166561395481288721898524257312,
       -1334570903704311347201177471950935055938544194648],
]

def dyadic(num: int, bits: int = DEN_BITS) -> arb:
    return arb(fmpz(num)) / (arb(2) ** bits)

mid = [dyadic(n) for n in ROOT_NUM]
# A single rational radius makes the box readily reproducible.
radius = arb(1) / (arb(2) ** 30)
box = [arb(x, radius) for x in mid]
Y = [[dyadic(INV_HESS_NUM[i][j]) for j in range(3)] for i in range(3)]

def mmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), arb(0))
             for j in range(len(B[0]))] for i in range(len(A))]

def mvec(A, x):
    return [sum((A[i][j] * x[j] for j in range(len(x))), arb(0))
            for i in range(len(A))]

def det3(A):
    return (A[0][0] * (A[1][1]*A[2][2] - A[1][2]*A[2][1])
          - A[0][1] * (A[1][0]*A[2][2] - A[1][2]*A[2][0])
          + A[0][2] * (A[1][0]*A[2][1] - A[1][1]*A[2][0]))

def positive(x):
    return x.lower() > 0

def negative(x):
    return x.upper() < 0

det_Y = det3(Y)
assert positive(det_Y) or negative(det_Y), "preconditioner is singular"

# Domain checks ensure that every logarithm evaluated below has positive input.
prob_box = [f(*box) for f in prob_f]
assert all(positive(x) for x in prob_box), "a probability interval meets zero"
sqrt3 = arb(3).sqrt()
assert positive(box[0]) and positive(box[1]) and positive(box[2])
assert negative(box[0] - box[1])
assert positive(arb(1) - box[1] - 2*sqrt3*box[2])

F0 = eval_gradient(*mid)
JX = eval_hessian(*box)

# Krawczyk operator K(x0,X) = x0-YF(x0)+(I-YJ(X))(X-x0).
YF = mvec(Y, F0)
YJ = mmul(Y, JX)
M = [[(arb(1) if i == j else arb(0)) - YJ[i][j]
      for j in range(3)] for i in range(3)]
D = [arb(0, radius) for _ in range(3)]
MD = mvec(M, D)
K = [mid[i] - YF[i] + MD[i] for i in range(3)]
krawczyk_ok = all(box[i].contains_interior(K[i]) for i in range(3))
assert krawczyk_ok, "Krawczyk image is not in the interior of the box"

# Hessian is negative definite everywhere in the box.  For A=-H, verify the
# three leading principal minors are positive (Sylvester criterion).
A = [[-JX[i][j] for j in range(3)] for i in range(3)]
minor1 = A[0][0]
minor2 = A[0][0]*A[1][1] - A[0][1]*A[1][0]
minor3 = det3(A)
assert positive(minor1) and positive(minor2) and positive(minor3)

gap_box = eval_delta(*box)
assert positive(gap_box)

def line(label, x):
    print(f"{label}: {x.str(30) if isinstance(x, arb) else x}")

print("PROVED_LOCAL")
print(f"precision_bits: {ctx.prec}")
print(f"box_radius: 2^-30 = {radius}")
for i, name in enumerate(("c", "d", "s")):
    line(f"box_{name}", box[i])
    line(f"box_{name}_lower", box[i].lower())
    line(f"box_{name}_upper", box[i].upper())
    line(f"krawczyk_{name}", K[i])
line("min_probability_enclosure_lower", min((x.lower() for x in prob_box), key=float))
line("feasibility_margin", arb(1) - box[1] - 2*sqrt3*box[2])
line("preconditioner_determinant", det_Y)
line("minus_hessian_minor_1", minor1)
line("minus_hessian_minor_2", minor2)
line("minus_hessian_minor_3", minor3)
line("gap_on_box", gap_box)
line("gap_on_box_lower", gap_box.lower())
line("gap_on_box_upper", gap_box.upper())
print("proved: exactly one stationary point lies in the box")
print("proved: Hessian is negative definite throughout the box")
print("proved: that stationary point is the unique maximizer on the box")
