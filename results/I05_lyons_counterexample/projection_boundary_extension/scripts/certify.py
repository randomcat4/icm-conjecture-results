#!/usr/bin/env python3
"""Exact, standard-library certificate for the frozen real/complex comparison.

Run: python certify.py [--output certificate.json]
No floating-point arithmetic, numerical eigensolver, imported old certificate,
or network access is used. Logs use 80 rational atanh terms; every interval
operation rounds outwards. This is an author self-check, not independent review.
"""
from __future__ import annotations
import argparse
import json
from fractions import Fraction as F
from itertools import combinations, permutations
from math import factorial
from pathlib import Path

K_INT = [[2719,-3449,1009,-818,-2490],[-3449,5840,-2889,-1513,1322],
         [1009,-2889,2148,2505,1101],[-818,-1513,2505,4698,3958],
         [-2490,1322,1101,3958,4596]]
B_INT = [[0,-3,-13,-16,4],[3,0,-4,15,17],[13,4,0,6,-15],
         [16,-15,-6,0,1],[-4,-17,15,-1,0]]
BASE = '9c27d3f25cb68ae8fe1d6596ffbb49d023d45892'
NLOG = 80
GRID = 10**40
Matrix = list[list[F]]
Interval = tuple[F, F]
ZERO: Interval = (F(0), F(0))


def need(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def floor(x: F) -> int:
    return x.numerator // x.denominator


def ceil(x: F) -> int:
    return -floor(-x)


def outward(a: F, b: F) -> Interval:
    need(a <= b, 'reversed interval')
    return F(floor(a*GRID), GRID), F(ceil(b*GRID), GRID)


def add(a: Interval, b: Interval) -> Interval:
    return outward(a[0]+b[0], a[1]+b[1])


def neg(a: Interval) -> Interval:
    return -a[1], -a[0]


def mul(a: Interval, b: Interval) -> Interval:
    v = [x*y for x in a for y in b]
    return outward(min(v), max(v))


def scale(a: Interval, b: F | int) -> Interval:
    return mul(a, (F(b), F(b)))


def div(a: Interval, b: Interval) -> Interval:
    need(b[0] > 0 or b[1] < 0, 'division interval contains zero')
    return mul(a, (1/b[1], 1/b[0]))


def atanh_log(y: F) -> Interval:
    need(1 <= y <= 2, 'log range reduction failed')
    z = (y-1)/(y+1)
    zz = z*z
    power = z
    total = F(0)
    for k in range(NLOG):
        total += 2*power/(2*k+1)
        power *= zz
    remainder = 2*power/((2*NLOG+1)*(1-zz))
    return outward(total, total+remainder)


def log_interval(x: F, log2: Interval) -> Interval:
    need(x > 0, 'log argument is not positive')
    y, k = x, 0
    while y < 1:
        y *= 2
        k -= 1
    while y >= 2:
        y /= 2
        k += 1
    return add(atanh_log(y), scale(log2, k))


def eye(n: int) -> Matrix:
    return [[F(i == j) for j in range(n)] for i in range(n)]


def matmul(a: Matrix, b: Matrix) -> Matrix:
    n, m, p = len(a), len(b), len(b[0])
    return [[sum((a[i][k]*b[k][j] for k in range(m)), F(0))
             for j in range(p)] for i in range(n)]


def trace(a: Matrix) -> F:
    return sum((a[i][i] for i in range(len(a))), F(0))


def trace_product(a: Matrix, b: Matrix) -> F:
    return sum((a[i][j]*b[j][i] for i in range(len(a))
                for j in range(len(a))), F(0))


def inverse_det(a: Matrix) -> tuple[Matrix, F]:
    n = len(a)
    x = [r[:] + e for r, e in zip(a, eye(n))]
    determinant = F(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if x[i][j]), None)
        need(pivot is not None, 'singular matrix')
        if pivot != j:
            x[j], x[pivot] = x[pivot], x[j]
            determinant = -determinant
        d = x[j][j]
        determinant *= d
        x[j] = [v/d for v in x[j]]
        for i in range(n):
            if i != j:
                c = x[i][j]
                x[i] = [u-c*v for u, v in zip(x[i], x[j])]
    return [row[n:] for row in x], determinant


def event_matrix(k: Matrix, mask: int) -> Matrix:
    a = [row[:] for row in k]
    for i in range(len(k)):
        a[i][i] -= 1-((mask >> i) & 1)
    return a


def determinant_polynomial(m: Matrix, a: Matrix) -> list[F]:
    """Independent Leibniz polynomial det(m+x*a), including every coefficient."""
    n = len(m)
    out = [F(0)]*(n+1)
    for perm in permutations(range(n)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(n) for j in range(i+1,n))
        poly = [F(sign)]
        for i, j in enumerate(perm):
            nxt = [F(0)]*(len(poly)+1)
            for d, c in enumerate(poly):
                nxt[d] += c*m[i][j]
                nxt[d+1] += c*a[i][j]
            poly = nxt
        out = [u+v for u,v in zip(out, poly)]
    return out


def leibniz_jet2(m: Matrix, pairs: list[tuple[int, int]]) -> tuple[list[F], Matrix]:
    """All real first and mixed second determinant derivatives, without inverses."""
    n, size = len(m), len(pairs)
    lookup = {pair: a for a, pair in enumerate(pairs)}
    first = [F(0)]*size
    second = [[F(0) for _ in range(size)] for _ in range(size)]
    for perm in permutations(range(n)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(n) for j in range(i+1,n))
        variables = [lookup[tuple(sorted((i,perm[i])))] for i in range(n)]
        for i in range(n):
            value = F(sign)
            for row in range(n):
                if row != i:
                    value *= m[row][perm[row]]
            first[variables[i]] += value
        for i,j in combinations(range(n),2):
            value = F(sign)
            for row in range(n):
                if row != i and row != j:
                    value *= m[row][perm[row]]
            a,b = variables[i],variables[j]
            second[a][b] += value
            second[b][a] += value
    return first, second


def imaginary_cycle_second(m: Matrix, b: Matrix) -> dict[int, F]:
    """Second derivatives grouped by the cycle containing the two B factors.

    Key 0 means different cycles; it cancels exactly for a skew direction.
    All other cycles in a permutation are retained, including their K factors.
    """
    n = len(m)
    out = {length: F(0) for length in range(n+1)}
    for perm in permutations(range(n)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(n) for j in range(i+1,n))
        owner, lengths = {}, {}
        for i in range(n):
            if i not in owner:
                u, length = i, 0
                while u not in owner:
                    owner[u] = i
                    length += 1
                    u = perm[u]
                lengths[i] = length
        for i,j in combinations(range(n),2):
            value = -2*sign*b[i][perm[i]]*b[j][perm[j]]
            for row in range(n):
                if row != i and row != j:
                    value *= m[row][perm[row]]
            length = lengths[owner[i]] if owner[i] == owner[j] else 0
            out[length] += value
    return out


def derivative_ratios(t: Matrix) -> list[F]:
    powers, current = [], eye(len(t))
    for _ in range(4):
        current = matmul(current, t)
        powers.append(trace(current))
    s1,s2,s3,s4 = powers
    return [s1, s1*s1-s2, s1**3-3*s1*s2+2*s3,
            s1**4-6*s1*s1*s2+3*s2*s2+8*s1*s3-6*s4]


def entropy_234(events: list[dict], direction: Matrix, imaginary: bool) -> list[Interval]:
    hs = [ZERO, ZERO, ZERO]
    for event in events:
        p, r, lp = event['p'], event['R'], event['logp']
        ratios = derivative_ratios(matmul(r, direction))
        if imaginary:
            need(ratios[0] == ratios[2] == 0, 'odd skew derivative')
            ratios = [F(0), -ratios[1], F(0), ratios[3]]
        a,b,c,d = [p*v for v in ratios]
        constants = [-a*a/p, -3*a*b/p+a**3/p**2,
                     -4*a*c/p-3*b*b/p+6*a*a*b/p**2-2*a**4/p**3]
        for i, derivative in enumerate([b,c,d]):
            hs[i] = add(hs[i], add(scale(lp,-derivative),
                                    (constants[i], constants[i])))
    return hs


def interval_ldlt(a: list[list[Interval]]) -> list[Interval]:
    n = len(a)
    l = [[ZERO for _ in range(n)] for _ in range(n)]
    d: list[Interval] = []
    for j in range(n):
        pivot = a[j][j]
        for k in range(j):
            pivot = add(pivot, neg(mul(mul(l[j][k], l[j][k]), d[k])))
        need(pivot[0] > 0, f'nonpositive LDL pivot {j}: {pivot}')
        d.append(pivot)
        l[j][j] = (F(1),F(1))
        for i in range(j+1,n):
            v = a[i][j]
            for k in range(j):
                v = add(v, neg(mul(mul(l[i][k],l[j][k]),d[k])))
            l[i][j] = div(v, pivot)
    return d


def decimal_bound(x: F, places: int, upper: bool) -> str:
    g = 10**places
    q = ceil(x*g) if upper else floor(x*g)
    sign = '-' if q < 0 else ''
    q = abs(q)
    return f'{sign}{q//g}.{q%g:0{places}d}'


def display(a: Interval, places: int = 18) -> list[str]:
    return [decimal_bound(a[0],places,False), decimal_bound(a[1],places,True)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('certificate.json'))
    args = parser.parse_args()
    k = [[F(x,10000) for x in row] for row in K_INT]
    b = [[F(x,50) for x in row] for row in B_INT]
    n = len(k)
    need(all(k[i][j] == k[j][i] and b[i][j] == -b[j][i]
             for i in range(n) for j in range(n)), 'symmetry failure')
    sylvester = []
    for a in [k, [[F(i == j)-k[i][j] for j in range(n)] for i in range(n)]]:
        minors = [inverse_det([row[:j] for row in a[:j]])[1] for j in range(1,n+1)]
        need(all(v > 0 for v in minors), 'not an interior contraction')
        sylvester.append([str(v) for v in minors])
    log2 = atanh_log(F(2))
    events: list[dict] = []
    summed_poly = [F(0)]*(n+1)
    h4_log, h4_fisher = ZERO, F(0)
    cycle_h2 = {length: ZERO for length in range(n+1)}
    for mask in range(1 << n):
        m = event_matrix(k,mask)
        r,detm = inverse_det(m)
        sign = (-1)**(n-mask.bit_count())
        p = sign*detm
        need(p > 0, 'nonpositive event')
        coeff = [sign*v for v in determinant_polynomial(m,b)]
        need(all(coeff[j] == 0 for j in range(1,n+1,2)), 'odd polynomial term')
        rb = matmul(r,b)
        rb2 = matmul(rb,rb)
        r2,r4 = trace(rb2),trace_product(rb2,rb2)
        d2,d4 = p*r2,p*(3*r2*r2-6*r4)
        need(coeff[0] == p and -2*coeff[2] == d2 and 24*coeff[4] == d4,
             'Leibniz/trace derivative mismatch')
        lp = log_interval(p,log2)
        cycle_parts = imaginary_cycle_second(m,b)
        need(cycle_parts[0] == cycle_parts[1] == 0, 'odd marked cycle did not cancel')
        need(sign*sum(cycle_parts.values()) == d2, 'cycle decomposition mismatch')
        for length, value in cycle_parts.items():
            cycle_h2[length] = add(cycle_h2[length],scale(lp,-sign*value))
        events.append(dict(p=p,R=r,logp=lp,p2=d2,p4=d4))
        summed_poly = [u+v for u,v in zip(summed_poly,coeff)]
        h4_log = add(h4_log,scale(lp,-d4))
        h4_fisher -= 3*d2*d2/p
    need(summed_poly == [F(1)]+[F(0)]*n, 'polynomial normalization')
    # Normalization differentiates to sum_S p_S R_S = 0 entrywise.
    need(all(sum(e['p']*e['R'][i][j] for e in events) == 0
             for i in range(n) for j in range(n)), 'matrix normalization')
    hi = entropy_234(events,b,True)
    need(hi[0][0] > F(1,200) and hi[0][1] < F(1,190), 'old H2 not reproduced')
    need(hi[1][0] <= 0 <= hi[1][1], 'H3 parity mismatch')
    need(hi[2][0] > -64154777 and hi[2][1] < -64154776, 'H4 bracket')
    # Full real Hessian in the unnormalised symmetric coordinate basis.
    pairs = [(i,j) for i in range(n) for j in range(i,n)]
    basis: list[Matrix] = []
    for i,j in pairs:
        e = [[F(0) for _ in range(n)] for _ in range(n)]
        e[i][j] = e[j][i] = F(1)
        basis.append(e)
    size = len(basis)
    h = [[ZERO for _ in range(size)] for _ in range(size)]
    q = [[ZERO for _ in range(size)] for _ in range(size)]
    fisher = [[F(0) for _ in range(size)] for _ in range(size)]
    for mask,event in enumerate(events):
        p,r,lp = event['p'],event['R'],event['logp']
        re = [matmul(r,e) for e in basis]
        v = [trace(x) for x in re]
        direct_first,direct_second = leibniz_jet2(event_matrix(k,mask),pairs)
        sign = (-1)**(n-mask.bit_count())
        need(all(sign*direct_first[i] == p*v[i] for i in range(size)), 'full Jacobian check')
        for i in range(size):
            for j in range(i+1):
                f = p*v[i]*v[j]
                second = p*(v[i]*v[j]-trace_product(re[i],re[j]))
                need(sign*direct_second[i][j] == second, 'full mixed derivative check')
                fisher[i][j] += f
                q[i][j] = add(q[i][j],scale(lp,-second))
    for i in range(size):
        for j in range(i+1):
            h[i][j] = add(q[i][j],(-fisher[i][j],-fisher[i][j]))
            h[j][i] = h[i][j]
    weights = [1 if i == j else 2 for i,j in pairs]
    shifted = [[add(neg(h[i][j]), (-F(16*weights[i]),-F(16*weights[i])))
                if i == j else neg(h[i][j]) for j in range(size)] for i in range(size)]
    pivots = interval_ldlt(shifted)
    # A real direction can have a positive determinant/log part.
    first,second = pairs.index((0,0)),pairs.index((1,1))
    real_q = add(add(q[first][first],q[second][second]),scale(q[second][first],2))
    real_f = fisher[first][first]+fisher[second][second]+2*fisher[second][first]
    need(real_q[0] > 24 and real_q[1] < 25, 'positive real log part not certified')
    # Rational rank-two decomposition, pivot (0,3); no spectral decomposition.
    u,v = 0,3
    b1 = [[(b[i][u]*b[j][v]-b[i][v]*b[j][u])/b[u][v]
           for j in range(n)] for i in range(n)]
    b2 = [[b[i][j]-b1[i][j] for j in range(n)] for i in range(n)]
    for a in [b1,b2]:
        # Vanishing four-index Pluecker relations plus a nonzero entry gives rank 2.
        need(any(a[i][j] != 0 for i in range(n) for j in range(n)), 'zero summand')
        need(all(a[i][j]*a[l][m]-a[i][l]*a[j][m]+a[i][m]*a[j][l] == 0
                 for i,j,l,m in combinations(range(n),4)), 'summand not rank two')
    h1,h2 = entropy_234(events,b1,True)[0],entropy_234(events,b2,True)[0]
    cross = scale(add(add(hi[0],neg(h1)),neg(h2)),F(1,2))
    need(h1[1] < 0 and h2[1] < 0 and cross[0] > 0, 'interference signs')
    # Realification H4, evaluated on all 32 x 32 product events without a 10! expansion.
    lift4 = F(0)
    rbs = [matmul(e['R'],b) for e in events]
    for i,e in enumerate(events):
        for j,f in enumerate(events):
            tau = trace_product(rbs[i],rbs[j])
            lift4 -= 12*e['p']*f['p']*tau*tau
    need(-120599633 < lift4 < -120599632, 'realification H4 bracket')
    # A non-real triangle product obstructs every same-coordinate real DPP kernel.
    linear = b[0][1]*k[1][2]*k[2][0]+k[0][1]*b[1][2]*k[2][0]+k[0][1]*k[1][2]*b[2][0]
    cubic = -b[0][1]*b[1][2]*b[2][0]
    step = F(1,100000)
    need(linear*step+cubic*step**3 != 0, 'triangle obstruction disappeared')
    # A deterministic real-direction derivative check, independent polynomial versus trace.
    real_a = [[b[min(i,j)][max(i,j)] for j in range(n)] for i in range(n)]
    for mask,event in enumerate(events):
        poly = determinant_polynomial(event_matrix(k,mask),real_a)
        ratios = derivative_ratios(matmul(event['R'],real_a))
        sign = (-1)**(n-mask.bit_count())
        need(all(sign*factorial(d)*poly[d] == event['p']*ratios[d-1]
                 for d in range(1,5)), 'real derivative polynomial mismatch')
    real234 = entropy_234(events,real_a,False)
    printed_h = [[tuple(F(x) for x in display(h[i][j],8)) for j in range(size)] for i in range(size)]
    printed_shifted = [[add(neg(printed_h[i][j]),(-F(16*weights[i]),-F(16*weights[i])))
                        if i == j else neg(printed_h[i][j]) for j in range(size)] for i in range(size)]
    printed_pivots = interval_ldlt(printed_shifted)
    result = {
        'status': 'EXACT_SELF_CHECK_PASS; INDEPENDENT_REVIEW_PENDING',
        'base_commit': BASE, 'K_integer':K_INT,'K_denominator':10000,
        'B_integer':B_INT,'B_denominator':50,
        'method': {'arithmetic':'Fraction only','atanh_terms':NLOG,
                   'outward_grid_denominator':str(GRID),'real_frobenius_margin':16},
        'sylvester_K_and_I_minus_K':sylvester,
        'normalization_and_Leibniz_trace_checks':True,
        'imaginary_H2_H3_H4': [display(x) for x in hi],
        'imaginary_H4_log_part':display(h4_log),
        'imaginary_H2_by_marked_cycle_length':{str(k):display(v) for k,v in cycle_h2.items()},
        'imaginary_H4_fisher_part':display((h4_fisher,h4_fisher)),
        'real_basis_upper_triangular_zero_based':pairs,
        'minus_real_Hessian_minus_16W_LDL_pivots':[display(x,12) for x in pivots],
        'real_Hessian_lower_triangle_intervals_8dp':
            [[display(h[i][j],8) for j in range(i+1)] for i in range(size)],
        'printed_matrix_LDL_pivots': [display(x,8) for x in printed_pivots],
        'real_E00_plus_E11':{'log_part':display(real_q),'fisher_penalty':display((real_f,real_f))},
        'rank_two_decomposition_pivot_zero_based':[u,v],
        'rank_two_H2': [display(h1),display(h2)], 'rank_two_cross_Hessian':display(cross),
        'realification_H2_H3':['0','0'],'realification_H4':display((lift4,lift4)),
        'triangle_012_imaginary_product':{'coefficient_t':str(linear),'coefficient_t3':str(cubic)},
        'real_symmetric_upper_B_test_H2_H3_H4':[display(x) for x in real234],
        'scope':'Fixed K real Hessian only; not global real concavity or real transference.'
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in [
        'status','imaginary_H2_H3_H4','minus_real_Hessian_minus_16W_LDL_pivots',
        'rank_two_H2','rank_two_cross_Hessian','imaginary_H2_by_marked_cycle_length','realification_H4',
        'triangle_012_imaginary_product']},ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
