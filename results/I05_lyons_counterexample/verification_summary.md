# Public verification summary

Status: **independently verified exact counterexample**.

## Statement checked

For the displayed rational (5\times5) matrix (K), rational real antisymmetric matrix (B), Hermitian direction (A=iB), and (h=1/100000), the matrices (K-hA) and (K+hA) are strict positive contractions and

\[
H(K)<\frac{H(K-hA)+H(K+hA)}2.
\]

This refutes the universal assertion in the original Lyons Conjecture 2.6 over finite-dimensional complex Hilbert spaces.

## Exact mathematical checks

Two independent mathematical reconstructions returned `CORRECT`. A later fresh manuscript reconstruction also returned `CORRECT` for the frozen paper and independently rebuilt the decisive calculations. The public checker reproduces the latter reconstruction using only `fractions.Fraction` and other Python standard-library modules.

The exact checks establish:

- (K^T=K), (B^T=-B), so (A=iB) is purely imaginary Hermitian; a nonzero (4\times4) principal determinant and the zero (5\times5) antisymmetric determinant give `rank(A) = 4`.
- Every leading principal minor needed by Sylvester's criterion is a positive rational for (K), (I-K), (K\pm hA), and (I-K\mp hA).
- All 32 event polynomials are reconstructed exactly. Their odd coefficients vanish, determinant and Möbius formulas agree at the midpoint and endpoints, all 96 evaluated probabilities are positive, each 32-event distribution sums to one, and the two endpoints agree term by term.
- The 32 center probabilities, 32 endpoint probabilities, and 32 second derivatives agree with the published TSV and TeX certificates.
- An 80-term positive `atanh` expansion with an explicit rational remainder gives

  \[
  \frac{5205948204140858303}{10^{21}}
  <D^2H(K)[A,A]<
  \frac{5205948204140858304}{10^{21}},
  \]

  which lies strictly between (1/200) and (1/190).
- The same directed rational interval method gives

  \[
  \frac{233555253804221762}{10^{30}}
  <H(K+hA)-H(K)<
  \frac{233555253804221763}{10^{30}},
  \]

  so the gap is strictly greater than (10^{-13}). Endpoint equality completes the midpoint-concavity violation.

## Manuscript checks

The frozen manuscript was reviewed against the exact data and the primary source statement. Its definitions, quantifiers, event formula, feasibility proof, entropy differentiation, interval bounds, finite-chord conclusion, source locations, and stated limitations were accepted. An independent no-shell-escape four-step build produced six A4 pages without unresolved citations, unresolved references, or box warnings.

## Boundaries

The verification does not cover a real-symmetric-only conjecture, minimality of dimension five, absolute novelty or priority, formal proof-assistant verification, authorship, or submission. It also does not rerun the numerical search that originally found the rational example; discovery history is unnecessary for checking the exact certificate.
