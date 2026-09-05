# Verification summary

Review date: 2026-09-05. Frozen mathematical source:
`8e2488d238f0fd79dab4adf51f1f9a863a00a112`.

Two fresh read-only reviews were performed separately.

## Analytic proof review

The reviewer checked the full quantifiers and proof chain for Theorem A, Corollary B,
Theorem C, and Theorem D in [proof_zh.md](proof_zh.md), including:

- the event vanishing order
  \(w_S=r+|S|-2\,\mathrm{rank}(F_S)\);
- the uniform boundary expansion and its operator-norm remainder;
- separation of the block-preserving and transverse logarithmic coefficients;
- the full-spark formulas for \(\Lambda_P(D)\) and \(T_P(E)\);
- positive definiteness of the real Fisher forms under coordinate-graph connectivity;
- the cross-component zero space in the disconnected case;
- the fourth-order polarization and strictness argument for independent blocks.

Verdict: **CORRECT**. No missing hypothesis, weakened quantifier, or unsupported local-to-global
step was found. The review does not promote the local real theorem to global real concavity.

## Exact certificate review

The second reviewer extracted the two scripts from the frozen source into an independent
temporary directory and ran both with a bundled Python interpreter. They use only Python's
standard library when kept together in the same directory.

The projection-boundary checker returned `PASS` with \(N=32\). Its rebuilt JSON was
byte-for-byte identical to the frozen certificate:

```text
SHA-256 537E9621C41D7667889AE43F50255B6B6F3E0877B070D00EAF046FE23EFC55FD
```

The review checked that the estimates cover every \(0<\varepsilon\le2^{-32}\), rather than
only the endpoint. The strict sufficient inequalities retain positive margins for the safe
probability interval, the real internal and transverse directions, and the imaginary error.
It also independently recomputed

\[
\ell(P,I-2P,\widehat B)
=0.03403629510449026624798215159924063818\ldots
\]

inside the directed rational interval stored in the certificate, and confirmed
\(1/30<\ell<1/29\).

Verdict: **CORRECT** for the explicit all-parameter family. The finite non-full-spark and mixed
fourth-order program checks are diagnostic cross-checks; the corresponding arbitrary-dimensional
claims rest on the analytic proof.

## Notation and matroid boundary

The proof defines \(\Lambda_P(D)\) for real block-preserving directions. It is not the imaginary
finite limit \(\ell(P,N,B)\), and the frozen files do not define a quantity literally named
\(\Lambda(P,B)\).

The rank function determines \(w_S\) and the zero/nonzero basis pattern. It does not determine
the Plücker magnitudes in \(a_S\), the probability jets \(b_S\), or the direction \(B\).
Consequently the ordinary underlying matroid is insufficient to recover
\(\ell(P,N,B)\), and the proof does not claim otherwise.

As a separate exact check, the reviewer considered rank-one projections \(P=uu^T\) on
\(\mathbb R^3\), with \(Bv=u\times v\). Direct Möbius inversion gives

\[
\ell(P,I-2P,B)
=2\sum_i\left[x_i\log x_i-(1-x_i)\log(1-x_i)\right],
\qquad x_i=u_i^2.
\]

The choices \(u=(1,1,1)/\sqrt3\) and \(u=(1,1,2)/\sqrt6\) have the same oriented
matroid \(U_{1,3}\) and the same \(\lVert B\rVert_F^2=2\), but give
\(\log(9/16)\) and
\(\log(144\,5^{2/3}/625)\), respectively. This is a strict representation-dependence
counterexample, rather than merely an absence of a matroidal formula in the proof.

These reviews do not certify absolute novelty, priority, formal proof-assistant verification,
authorship, submission, or the unresolved global real-symmetric restriction.
