# Analytic reduction of the symmetric six-point family

This directory records the structural explanation found after the original
counterexample and the shorter author-review note.

The family is

\[
C(c,d)=\frac13
\begin{pmatrix}
c+d&c-2d&c+d\\
c-2d&c+d&c+d\\
c+d&c+d&c-2d
\end{pmatrix},
\quad
K=\frac12\begin{pmatrix}I&C\\ C&I\end{pmatrix},
\quad
B=s\,\mathrm{diag}(S_0,-S_0),
\]

where

\[
S_0=\begin{pmatrix}0&1&-1\\-1&0&1\\1&-1&0\end{pmatrix}.
\]

Its useful features are

\[
\mathrm{spec}(C)=\{c,d,-d\},\qquad CS_0=-S_0C.
\]

The \(c\)-mode lies in \(\ker S_0\), while the \(d,-d\) modes carry the
imaginary motion.  Separating their boundary margins explains why the
two-parameter center improves the earlier locked choice \(c=d\).

## Proved results

- A coordinate group of order 12, together with complementation, partitions
  all 64 complete events into eight orbits of sizes
  \(2,12,12,2,12,6,12,6\).
- All eight center and endpoint probabilities are explicit integer
  polynomials after a common scaling.  This reduces stationarity to three
  seven-term logarithmic equations.
- With
  \(x=(1-c)/(1-d)\) and \(r=s/(1-d)\), the boundary expansion has leading
  term \((1-d)^2G(x,r)\).  Its small-\(r\) quadratic coefficient has a
  unique maximum at \(x=1+\sqrt3\).
- Arb ball arithmetic isolates one stationary point near
  \[
  (c,d,s)=(0.9849706685842502,\,
  0.9939745132358184,\,
  0.001114721194445069)
  \]
  and proves it is the unique maximizer on the explicit box in
  `uniqueness/box_definition.md`.
- The accompanying denominator-\(10^4\) rational point gives a rigorously
  certified entropy gap
  \[
  \frac1{3\,700\,000}<\Delta<\frac1{3\,600\,000}.
  \]

The interval result proves uniqueness only on the displayed box.  It does
not prove global uniqueness throughout the full feasible chamber, and the
family is not claimed to classify every six-point counterexample.

## Reading order

1. `analysis_zh.md` gives the Chinese synthesis.
2. `symbolic/class_formulas.md` and
   `symbolic/stationary_equations.md` contain the exact reduction.
3. `asymptotics/asymptotic_derivation.md` explains the passive and moving
   modes and proves the boundary expansion.
4. `uniqueness/proof_or_gap.md` explains the Krawczyk and strict-concavity
   certificates.
5. `exact_candidate/verification.md` records the independent exact check of
   the larger rational candidate.

## Reproduction

From this directory:

```console
uv run --offline --with sympy python symbolic/derive_classes.py
uv run --offline --with sympy python asymptotics/derive_asymptotics.py
uv run --offline --with sympy --with python-flint python uniqueness/certify_box.py
uv run --offline python exact_candidate/exact_rational_1e4.py
uv run --offline python exact_candidate/independent_check.py
```

The first two programs reconstruct the formulas from all 64 determinants.
The third performs the interval proof.  The final two are independent
standard-library checks of the rational candidate.
