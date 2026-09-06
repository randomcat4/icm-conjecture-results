# Proof and remaining gap

## Status: `PROVED_BOX`

There is a rigorously isolated stationary point in the explicit box `X_0`
defined in `box_definition.md`.  It is the unique global maximizer of
`Delta` on the larger explicit box `X`.  This is a computer-assisted proof
using Arb ball arithmetic, not a conclusion from optimizer convergence.

## Derivative representation

Write the sixteen center/endpoint probability polynomials as `r_l` and give
them signed multiplicity `w_l`: `+m_j` for `p_j` and `-m_j` for `q_j`.  Then

\[
\Delta=\sum_l w_l r_l\log r_l.
\]

For coordinates `x_i in {c,d,s}`, the verifier evaluates

\[
\partial_i\Delta=
 \sum_l w_l(1+\log r_l)\partial_i r_l,
\]

\[
\partial_{ij}\Delta=
 \sum_l w_l\left((1+\log r_l)\partial_{ij}r_l+
 \frac{\partial_i r_l\partial_j r_l}{r_l}\right),
\]

and the exact third-derivative identity

\[
\begin{aligned}
\partial_{ijk}(r\log r)
={}&(1+\log r)r_{ijk}\\
&+\frac{r_{ij}r_k+r_{ik}r_j+r_{jk}r_i}{r}
-\frac{r_i r_j r_k}{r^2}.
\end{aligned}
\]

All derivatives of `r_l` are derivatives of explicit rational
polynomials.  SymPy differentiates only these polynomials; Arb evaluates all
transcendental and interval operations at 256-bit precision.

## Local existence and uniqueness

Let `F=grad Delta`.  The verifier uses the exact dyadic point `x_hat` and an
exact dyadic approximation `Y` to the inverse Hessian.  It evaluates the
Krawczyk operator

\[
\mathcal K(\widehat x,X_0)=
\widehat x-YF(\widehat x)
+(I-YF'(X_0))(X_0-\widehat x).
\]

Arb proves coordinatewise strict inclusion

\[
\mathcal K(\widehat x,X_0)\subset \mathrm{int}(X_0).
\]

The sixteen probability intervals are strictly positive on `X_0`, so `F`
is continuously differentiable there.  The Krawczyk theorem therefore gives
existence and uniqueness of a zero of `F` in `X_0`.

For `A=-Hess Delta`, Arb also proves that all three leading principal minors
are strictly positive throughout `X_0`.  Sylvester's criterion proves that
the Hessian is negative definite there.

## Unique maximum on the larger box

Natural interval evaluation of the Hessian is needlessly wide because the
eight entropy contributions strongly cancel.  The larger-box verifier uses
a centered mean-value enclosure instead.  On each subbox `Z` with center
`z` and coordinate radii `rho_k`, it encloses

\[
H_{ij}(Z)\subset H_{ij}(z)+
[-1,1]\sum_k\rho_k\sup_Z|\partial_kH_{ij}|.
\]

The third derivatives in the right side are evaluated from the displayed
identity.  An adaptive exact-dyadic cover terminates with 128 subboxes.  On
every terminal subbox:

1. every probability interval is strictly positive;
2. all three leading principal minors of `-Hess Delta` have positive lower
   endpoints.

The smallest certified lower endpoints over the cover are approximately

\[
6.2962\times10^{-3},\qquad
1.1386\times10^{-4},\qquad
1.2796\times10^{-5}.
\]

Thus `Delta` is strictly concave on the convex box `X`.  The stationary point
already proved to lie in `X_0 subset X` is consequently the unique stationary
point and unique global maximizer on `X`.

## What is not proved

This run does **not** prove that the point is the global maximizer on the full
feasible chamber.  Repeated numerical convergence is not used for that
claim.

The exact missing statement is an upper bound on the complement:

\[
\sup\{\Delta(c,d,s):(c,d,s)\text{ feasible},\ (c,d,s)\notin X\}
<\Delta(x_*).
\]

A full proof must treat boundary collars separately because some `p_j` or
`q_j` tend to zero and derivatives of `r log r` become singular.  A viable
completion would combine:

- interval branch-and-bound on a compact interior where all probabilities
  have a uniform positive lower bound;
- analytic `x log x` bounds in collars near `s=0`, `c=0,1`, `d=0,1`, and
  `d+2 sqrt(3)s=1`;
- a certified lower bound for `Delta(x_*)` against which every exterior box
  is compared.

The last two bullets are the present bottleneck.  A naive continuity bound
is much larger than the target height, about `2.74724e-7`, and therefore
cannot exclude the complement.
