# Projection-boundary extension / 投影边界推广

This public research note develops a second mechanism behind the five-dimensional complex
counterexample. It works near projection kernels and separates real symmetric directions from
purely imaginary Hermitian directions. The full Chinese proof is in
[proof_zh.md](proof_zh.md); [mechanism_background_zh.md](mechanism_background_zh.md) derives
the interior second-, third-, and fourth-order formulas used by the extension.
The frozen mathematics and the exact certificate received separate fresh read-only checks;
their scope and limitations are recorded in [verification_summary.md](verification_summary.md).

The main statements are:

1. Near every connected real projection, every sufficiently small uniformly inward
   perturbation has strictly negative entropy Hessian in all nonzero real symmetric directions.
2. For every real projection, the isotropic regularization
   \(K_\varepsilon=\varepsilon I+(1-2\varepsilon)P\) is concave in all real directions for
   sufficiently small \(\varepsilon>0\), modulo cross-component zero directions.
3. One explicit rational rank-two projection and commuting skew direction satisfy, throughout
   the whole certified interval \(0<\varepsilon\le2^{-32}\),
   \[
   D^2H(K_\varepsilon)[C,C]<-\lVert C\rVert_F^2
   \]
   for every nonzero real symmetric \(C\), while
   \[
   D^2H(K_\varepsilon)[i\widehat B,i\widehat B]>\frac1{40}.
   \]
4. At a direct sum of two interior real kernels, the mixed fourth derivative in cross-block
   directions is a negative square expectation. This blocks one natural real-counterexample
   construction.

These results do not prove global concavity on the real-symmetric kernel domain. They also do
not give a real counterexample.

## What is actually known about \(\Lambda\)

Two different quantities must be kept separate.

For a real symmetric direction \(C=D+E\), split it relative to
\(\mathrm{ran}(P)\oplus\ker(P)\), with \(D\) block diagonal and \(E\) block off-diagonal.
The logarithmic coefficient in the projection-boundary Hessian expansion contains

\[
\Lambda_P(D)=\sum_{S\subseteq[n]}w_S\,D^2p_S(P)[D,D],
\qquad
w_S=r+|S|-2\,\mathrm{rank}(F_S).
\]

The exponents \(w_S\) are determined by the rank function of the represented matroid. For a
full-spark projection, the proof reduces the whole sum to

\[
\Lambda_P(D)=4\,\mathrm{tr}(U)\,\mathrm{tr}(V),
\]

where \(D=\mathrm{diag}(U,V)\) in the range/kernel splitting. For a non-full-spark projection,
the remaining expression is

\[
\Lambda_P(D)
=4\,\mathrm{tr}(U)\,\mathrm{tr}(V)
+2\!\!\sum_{\substack{|S|=r\\ \mathrm{rank}(F_S)=r-1}}
D^2p_S(P)[D,D].
\]

The extra terms come precisely from zero bases of defect two. The proof locates them and checks
one exact non-generic example, but it does not turn them into a representation-free invariant
of an abstract matroid.

The purely imaginary counterexample uses a different finite limit:

\[
\ell(P,N,B)=-\sum_S b_S(P,B)\log a_S(N),
\qquad
b_S(P,B)=
\left.\frac{d^2}{dt^2}p_S(P+itB)\right|_{t=0}.
\]

The cancellation \(\sum_S w_Sb_S=0\) removes the logarithmic divergence. For the explicit
rational family, the certificate proves

\[
\frac1{30}<\ell(P,I-2P,\widehat B)<\frac1{29},
\qquad
\ell=0.034036295104490266\ldots.
\]

Thus the extension computes the explicit example rigorously and gives a general finite formula.
It does **not** show that \(\ell\), sometimes informally denoted \(\Lambda(P,B)\), depends only on
the underlying matroid. The leading coefficients \(a_S\) contain Gram/Plücker magnitudes, and
\(b_S\) depends on the chosen representation and direction \(B\). The support matroid alone
does not contain this metric data.

There is a strict three-coordinate witness. Let \(P=uu^T\) have rank one, let
\(Bv=u\times v\), and write \(x_i=u_i^2\). Then

\[
\ell(P,I-2P,B)
=2\sum_{i=1}^3
\left[x_i\log x_i-(1-x_i)\log(1-x_i)\right].
\]

Both \(u=(1,1,1)/\sqrt3\) and \(u=(1,1,2)/\sqrt6\) represent the same oriented
uniform matroid \(U_{1,3}\), and both directions have \(\lVert B\rVert_F^2=2\). Nevertheless,

\[
\ell=\log\frac9{16}
\quad\text{and}\quad
\ell=\log\left(\frac{144\,5^{2/3}}{625}\right),
\]

respectively. The two values differ, proving that even the isotropic specialization cannot be
recovered from the ordinary or oriented matroid plus the normalization of \(B\).

## Exact reproduction

Only Python's standard library is required. From this directory, run:

```console
python scripts/certify_projection_boundary.py --output certificates/projection_boundary_certificate.json
```

The program uses rational arithmetic and directed rational logarithm intervals. It checks the
projection and commuting skew direction, all event vanishing orders and leading coefficients,
all 225 entries of the logarithmic Hessian form, an independent event-jet calculation at
\(\varepsilon=1/4\), the positive imaginary limit, the complete interval
\(0<\varepsilon\le2^{-32}\), one non-full-spark stratum, and the mixed fourth-order identity.

The earlier mechanism certificate can be rebuilt with:

```console
python scripts/certify.py --output certificates/mechanism_certificate.json
```

The checked outputs are in [certificates/](certificates/). Both scripts retain their original
frozen status strings so that the generated files can be compared directly with the research
artifacts.

## Verification and provenance

The mathematical source was frozen in private commit
`8e2488d238f0fd79dab4adf51f1f9a863a00a112`. The public files preserve its theorem statements,
explicit matrices, and exact algorithms; only repository-local links, publication metadata, and
unsupported Markdown macros were adapted. The exact certificate was rerun from a clean public
path before this package was pushed.

The bounded source review and these calculations do not certify absolute novelty, priority,
authorship, or a resolution of the real-symmetric restriction.
