# Larger entropy gap search: round 1

## Verdict

The frozen primary target \(\Delta>10^{-6}\) was not reached. The round is
therefore NO_HIT_BOUNDED, not a nonexistence result.

The round did produce a new exact six-point counterexample with
\[
\frac1{3{,}700{,}000}<\Delta<\frac1{3{,}600{,}000}.
\]
Its gap is approximately \(2.7433770261572459\times10^{-7}\), about 10.9
times the previously certified symmetric six-point gap. It has common matrix
denominator \(10^4\), strict spectral margin, 64 positive atoms, eventwise
equal endpoint laws, and the same eight probability classes.

## Exact promoted construction

\[
K=\frac1{10000}
\begin{pmatrix}
5000&0&0&3298&-1672&3298\\
0&5000&0&-1672&3298&3298\\
0&0&5000&3298&3298&-1672\\
3298&-1672&3298&5000&0&0\\
-1672&3298&3298&0&5000&0\\
3298&3298&-1672&0&0&5000
\end{pmatrix},
\]
\[
B=\frac1{10000}
\begin{pmatrix}
0&11&-11&0&0&0\\
-11&0&11&0&0&0\\
11&-11&0&0&0&0\\
0&0&0&0&-11&11\\
0&0&0&11&0&-11\\
0&0&0&-11&11&0
\end{pmatrix},
\qquad Q_\pm=K\pm iB.
\]

Writing \(K=\frac12\left(\begin{smallmatrix}I&C\\C&I\end{smallmatrix}\right)\),
the eigenvalues of \(C\) are
\(-497/500,1231/1250,497/500\). The exact checker proves
\(CS=-SC\), hence \(KB=BK\), and
\(\lVert B\rVert=11\sqrt3/10000<3/1000\). Thus both endpoints are
strict positive contractions.

The candidate checker and a fresh independent implementation agree on all 64
events, normalization, endpoint equality, the eight classes, and the directed
logarithm interval. The independent verifier marked the full frozen target
CRITICAL_GAPS only because the verified gap remains below \(10^{-6}\); it
found no defect in the exact secondary result.

## Route outcomes

- **A: direct finite chord.** 85,898 objective-level evaluations over
  \(n=6,7,8\). Best floating gap \(2.7472382897\times10^{-7}\); rounding to
  denominator \(10^4\) gave the exact construction above.
- **B: algebraic family.** The first involution-family work unit used at least
  76,117 finite-chord evaluations. A second 154,552-evaluation search of the
  relaxed three-parameter family
  \[
  C(a,b)=\begin{pmatrix}a&b&a\\b&a&a\\a&a&b\end{pmatrix},
  \qquad B=s\,\mathrm{diag}(S_0,-S_0)
  \]
  found a numerical optimum \(2.7472385457\times10^{-7}\). Six positive
  runs converged to this value within \(4.45\times10^{-16}\). This is
  numerical saturation evidence for the family, not a proved upper bound.
- **C: full imaginary Hessian.** A 750-center scan found 64 positive-Hessian
  centers; random centers gave small finite gains. Two local refinement rounds
  around the known center used 1,026 further center evaluations and reached
  \(1.6780097825\times10^{-7}\), confirmed at 90 digits.
- **D: unrestricted endpoints.** General Hermitian optimization, dimensional
  lifts, 96,000 random near-projection chords, and 145,200 permutation-orbit
  probes did not exceed the same six-point basin. Its best endpoints were
  related by a coordinate permutation and had gap
  \(2.7472319647\times10^{-7}\).

## Structural conclusion

The old symmetric point is not locally optimal, but its natural enlargement
has a stable ceiling near \(2.74724\times10^{-7}\). Direct optimization, the
relaxed anticommuting family, and unrestricted endpoints independently
converged to the same scale. More restarts in this basin are not justified.
A further order-of-magnitude improvement will require a different mechanism:
multiple active imaginary modes with rare-event control, a center with
positive curvature away from the projection boundary, or a noncommuting
finite-chord construction.
