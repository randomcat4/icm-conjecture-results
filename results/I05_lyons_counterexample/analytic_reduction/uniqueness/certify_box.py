#!/usr/bin/env python3
"""Certify strict concavity on a macroscopic box around the stationary point.

Run with
    uv run --offline --with sympy --with python-flint python certify_box.py

This imports certify_local.py, so its Krawczyk certificate is checked first.
The larger box is then covered by exact dyadic subboxes.  On every subbox we
verify positivity of all probability formulae and Sylvester's criterion for
minus the Hessian.
"""
from __future__ import annotations

import certify_local as L
from flint import arb

RADII = [arb(1)/(arb(2)**13),
         arb(1)/(arb(2)**14),
         arb(1)/(arb(2)**16)]

worst_p = None
worst_m1 = None
worst_m2 = None
worst_m3 = None
checked = 0
split_count = 0
max_axis_level = 0
stack = [(L.mid, RADII, (0,0,0))]
while stack:
      centers, radii, levels = stack.pop()
      max_axis_level = max(max_axis_level, *levels)
      X = [arb(centers[z], radii[z]) for z in range(3)]
      probs = []
      for _, r0, r1, _, _ in L.terms:
        p0 = r0(*centers)
        perr = sum((arb(r1[z](*X)).abs_upper()*radii[z]
                    for z in range(3)), arb(0))
        probs.append(arb(p0, perr))
      ok = all(L.positive(z) for z in probs)
      if ok:
        # Centered mean-value enclosure of the Hessian.  Direct natural
        # interval evaluation loses too much to cancellation near the upper
        # corner of this box.
        H0 = L.eval_hessian(*centers)
        T = L.eval_third(*X, probabilities_override=probs)
        J = [[None for b in range(3)] for a in range(3)]
        for a in range(3):
          for b in range(3):
            err = sum((T[a][b][z].abs_upper()*radii[z]
                       for z in range(3)), arb(0))
            J[a][b] = arb(H0[a][b], err)
        A = [[-J[a][b] for b in range(3)] for a in range(3)]
        m1 = A[0][0]
        m2 = A[0][0]*A[1][1] - A[0][1]*A[1][0]
        m3 = L.det3(A)
        ok = L.positive(m1) and L.positive(m2) and L.positive(m3)
        if checked == 0 and split_count == 0:
          print("initial_parent_minors",m1.str(20),m2.str(20),m3.str(20),flush=True)
      if not ok:
        axis = min(range(3), key=lambda z: levels[z])
        if levels[axis] >= 12:
          print("FAILED_CELL", levels)
          for name,z,r in zip(("c","d","s"),centers,radii):
            print(name,z.str(25),r.str(10))
          print("probability_intervals_positive", all(L.positive(z) for z in probs))
          if not all(L.positive(z) for z in probs):
            for idx,z in enumerate(probs):
              if not L.positive(z): print("bad_probability",idx,z.str(30))
          else:
            print("minors",m1.str(30),m2.str(30),m3.str(30))
          raise AssertionError("adaptive cover exceeded depth limit")
        new_radii = list(radii)
        new_radii[axis] = radii[axis] / 2
        new_levels = list(levels)
        new_levels[axis] += 1
        for sign in (-1, 1):
          new_centers = list(centers)
          new_centers[axis] = centers[axis] + sign*new_radii[axis]
          stack.append((new_centers, new_radii, tuple(new_levels)))
        split_count += 1
        continue
      pmin = min((z.lower() for z in probs), key=float)
      worst_p = pmin if worst_p is None or float(pmin) < float(worst_p) else worst_p
      worst_m1 = m1.lower() if worst_m1 is None or float(m1.lower()) < float(worst_m1) else worst_m1
      worst_m2 = m2.lower() if worst_m2 is None or float(m2.lower()) < float(worst_m2) else worst_m2
      worst_m3 = m3.lower() if worst_m3 is None or float(m3.lower()) < float(worst_m3) else worst_m3
      checked += 1

parent = [arb(L.mid[i], RADII[i]) for i in range(3)]
sqrt3 = arb(3).sqrt()
assert L.positive(parent[0]) and L.positive(parent[1]) and L.positive(parent[2])
assert L.negative(parent[0]-parent[1])
assert L.positive(arb(1)-parent[1]-2*sqrt3*parent[2])

print("PROVED_BOX")
print(f"terminal_subboxes_checked: {checked}")
print(f"adaptive_splits: {split_count}")
print(f"maximum_axis_subdivision_level: {max_axis_level}")
for name, X in zip(("c","d","s"), parent):
    print(f"parent_box_{name}: {X.str(30)}")
print(f"worst_probability_lower_bound: {worst_p.str(20)}")
print(f"worst_minus_hessian_minor_1_lower_bound: {worst_m1.str(20)}")
print(f"worst_minus_hessian_minor_2_lower_bound: {worst_m2.str(20)}")
print(f"worst_minus_hessian_minor_3_lower_bound: {worst_m3.str(20)}")
print("proved: Hessian is negative definite throughout the parent box")
print("proved: the Krawczyk root is the unique stationary point in the parent box")
print("proved: the Krawczyk root is the unique global maximizer on the parent box")
print("not proved: exclusion of larger values outside the parent box")
