# Larger-gap exact counterexamples

This directory contains two exact replacements for the very small finite chord in the five-page author-review note.  In both constructions the midpoint `K` is real symmetric, the direction is purely imaginary Hermitian, and the two endpoint laws agree event by event.

## Five labelled points

The file [`certify_optimized_five.py`](certify_optimized_five.py) contains rational matrices `K` and `B`, both with denominator `10^6`, and sets

\[
Q_+=K+iB,\qquad Q_-=K-iB.
\]

Exact Sylvester certificates prove that both endpoints are strict positive contractions.  Directed rational logarithm bounds prove

\[
\frac1{130\,000\,000}
<
H(Q_+)-H(K)
<
\frac1{120\,000\,000}.
\]

The gap is approximately `7.965330164295168e-9`, about `2.83e10` times the gap displayed in the five-page review note.  The smallest endpoint event probability is approximately `4.58e-9`.

This option retains the five-point dimension, but its matrices are best regarded as a compact exact certificate rather than a conceptual construction.

## Six labelled points with a symmetric construction

Put

\[
R=\frac13
\begin{pmatrix}
2&-1&2\\
-1&2&2\\
2&2&-1
\end{pmatrix},\qquad
P=\frac12
\begin{pmatrix}
I_3&R\\
R&I_3
\end{pmatrix},
\]

and

\[
S=\begin{pmatrix}
0&1&-1\\
-1&0&1\\
1&-1&0
\end{pmatrix},\qquad
B=\begin{pmatrix}
S&0\\
0&-S
\end{pmatrix}.
\]

The matrix `R` is a symmetric orthogonal reflection.  Hence `P` is a rank-three orthogonal projection.  Direct multiplication gives `B^T=-B` and `PB=BP`.  Define

\[
\varepsilon=\frac1{475},\qquad
K=\varepsilon I_6+(1-2\varepsilon)P,
\qquad
t=\frac1{2300},
\qquad
Q_\pm=K\pm itB.
\]

The file [`certify_symmetric_six.py`](certify_symmetric_six.py) proves

\[
\frac1{40\,000\,000}
<
H(Q_+)-H(K)
<
\frac1{39\,000\,000}.
\]

The gap is approximately `2.517182345129681e-8`, about `8.94e10` times the review-note gap.  The 64 complete events collapse to only eight distinct center/endpoint probability pairs.  This option is one dimension larger but has the shorter structural description and the larger gap.

## Reproduction

Both scripts use only Python's standard library and exact `Fraction` arithmetic:

```console
python certify_optimized_five.py
python certify_symmetric_six.py
```

Each script reconstructs every complete-event probability, proves strict positivity and normalization, checks eventwise endpoint equality, and evaluates the entropy difference with directed rational logarithm intervals.  Independent verification is summarized in [`verification_summary.md`](verification_summary.md).

Block-diagonal direct sums add entropy gaps.  Thus the six-point example can be amplified to any prescribed absolute gap, at the cost of increasing the dimension.  The compact six-point construction above is the preferred nontrivial example.
