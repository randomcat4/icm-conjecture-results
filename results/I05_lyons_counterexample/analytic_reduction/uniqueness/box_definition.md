# Certified boxes

## Function and chamber

The function is

\[
\Delta(c,d,s)=\sum_{j=1}^{8}m_j
   \{p_j(c,d)\log p_j(c,d)-q_j(c,d,s)\log q_j(c,d,s)\},
\]

with multiplicities

\[
(m_1,\ldots,m_8)=(2,12,12,2,12,6,12,6).
\]

The exact polynomials `p_j,q_j` are in `formulae.json`.  The domain is the
open chamber

\[
0<c<d<1,\qquad s>0,\qquad d+2\sqrt3s<1.
\]

## Exact dyadic center

All interval endpoints used by the certificates are dyadic.  Let
\(D=2^{160}\) and

\[
\widehat x=D^{-1}\begin{pmatrix}
1439536244858795825703648138477772648244438200086\\
1452695378559335846296495333266843919245350418322\\
1629166850848927691821404480252899240504455452
\end{pmatrix}.
\]

Numerically,

\[
\widehat x\approx
(0.98497066858425021191,\ 0.99397451323581841410,\
0.0011147211944450686890).
\]

## Krawczyk box

The local box is

\[
X_0=\widehat x+[-2^{-30},2^{-30}]^3.
\]

`certify_local.py` proves that `grad Delta` has exactly one zero in `X_0`
and that the Hessian is negative definite throughout `X_0`.

## Concavity box

The larger box is

\[
X=\widehat x+
[-2^{-13},2^{-13}]\times
[-2^{-14},2^{-14}]\times
[-2^{-16},2^{-16}].
\]

Its approximate coordinate ranges are

\[
\begin{aligned}
0.98484859827&<c<0.98509273890,\\
0.99391347808&<d<0.99403554839,\\
0.001099462405&<s<0.001129979984.
\end{aligned}
\]

`certify_box.py` proves strict concavity on all of `X`.  Together with the
Krawczyk zero in `X_0 subset X`, this proves that the zero is the unique
stationary point and unique global maximizer **on `X`**.

No claim is made about the complement of `X` in the feasible chamber.
