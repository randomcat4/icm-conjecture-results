# Boundary asymptotics for the split six-point family

**Status: PROVED_ASYMPTOTIC.**  This note proves the exact eight-class
reduction and the boundary expansion stated below.  It does not assert global
optimality or uniqueness of the observed finite stationary point.

Throughout, `log` is the natural logarithm and

\[
h(u)=-u\log u,\qquad \Delta(c,d,s)=H(Q_+)-H(K).
\]

## 1. The passive and moving modes

Put \(e=(1,1,1)^T\).  Direct multiplication gives

\[
C e=c e,\qquad S_0e=0,\qquad C S_0=-S_0C.
\]

On \(e^\perp\), the eigenvalues of \(C\) are \(d,-d\), while
\(S_0^2=-3I\).  Thus the \(c\)-mode belongs to \(\ker S_0\) and is not moved
by the imaginary perturbation.  The two \(d\)-modes are the moving modes.
Since \(K\) commutes with \(iB(s)\), their eigenvalues can be read
simultaneously.  The passive pair is

\[
\frac{1-c}{2},\quad \frac{1+c}{2},
\]

and the four moving eigenvalues are

\[
\frac{1\pm d}{2}\pm\sqrt3s.
\]

Consequently the endpoint feasibility condition at the upper boundary is
exactly

\[
1-d>2\sqrt3s.
\]

This spectral decomposition already explains why the two regularizations
should be separated: \(1-c\) pads a mode which is invisible to the motion,
whereas \(1-d\) is the spectral margin actually spent by the motion.

## 2. Exact eight-class formula

Write \(z=s^2\).  For each class let \(p_j\) be its center probability and
\(\eta_j=q_j-p_j\) its endpoint displacement.  Exact symbolic determinants of
all 64 complete events give the following eight polynomial classes.  The
multiplicities sum to 64, and both \(\sum m_jp_j\) and
\(\sum m_j(p_j+\eta_j)\) are identically one.

| \(j\) | representative | \(m_j\) | \(p_j\) | \(\eta_j\) |
|---:|:---:|---:|:---|:---|
|1|\(\varnothing\)|2|\(-\frac{(c-1)(c+1)(d-1)^2(d+1)^2}{64}\)|\(-\frac{3z(c-1)(c+1)(-d^2+6z-1)}8\)|
|2|\(\{1\}\)|12|\(\frac{(d-1)(d+1)(3c^2d^2+c^2-d^2-3)}{192}\)|\(\frac{z(-3c^2d^2+18c^2z-c^2+d^2-6z-1)}8\)|
|3|\(\{1,2\}\)|12|\(-\frac{(d^2+1)(3c^2d^2-c^2+d^2-3)}{192}\)|\(-\frac{z(-3c^2d^2+18c^2z-c^2-d^2+6z+1)}8\)|
|4|\(\{1,2,3\}\)|2|\(\frac{(c^2+1)(d^2+1)^2}{64}\)|\(\frac{3z(c^2+1)(-d^2+6z-1)}8\)|
|5|\(\{1,4\}\)|12|\(-\frac{(3cd^2-2cd+c+d^2-2d+3)(3cd^2+2cd+c-d^2-2d-3)}{576}\)|\(-\frac{z(-9c^2d^2+54c^2z+3c^2+4cd+d^2-6z-3)}{24}\)|
|6|\(\{2,4\}\)|6|\(-\frac{(d-1)(d+1)(3cd-c-d+3)(3cd+c+d+3)}{576}\)|\(-\frac{z(-9c^2d^2+54c^2z+3c^2-8cd+d^2-6z-3)}{24}\)|
|7|\(\{1,2,4\}\)|12|\(\frac{9c^2d^4-2c^2d^2+c^2+8cd^3+8cd+d^4-2d^2+9}{576}\)|\(\frac{z(-9c^2d^2+54c^2z+3c^2-4cd-d^2+6z+3)}{24}\)|
|8|\(\{1,3,4\}\)|6|\(\frac{(d^2+1)(9c^2d^2+c^2-16cd+d^2+9)}{576}\)|\(\frac{z(-9c^2d^2+54c^2z+3c^2+8cd-d^2+6z+3)}{24}\)|

The accompanying script reconstructs these formulas from the six-by-six
determinants; the table is not entered as input.

## 3. Boundary scaling of the atoms

Set

\[
1-c=\alpha=a\tau,\qquad 1-d=\beta=b\tau,
\qquad s=\rho\tau,
\]

where \(a>b>2\sqrt3\rho>0\) are fixed and \(\tau\downarrow0\).  The leading
terms of the exact class formulas are:

| \(j\) | \(m_j\) | center \(p_j\) | endpoint \(q_j\) | \(q_j-p_j\) |
|---:|---:|:---|:---|:---|
|1|2|\(\frac{ab^2}{8}\tau^3\)|\(\frac{a(b^2-12\rho^2)}8\tau^3\)|\(-\frac32a\rho^2\tau^3\)|
|2|12|\(\frac{b(2a+b)}{24}\tau^2\)|\(\frac{2ab+b^2-12\rho^2}{24}\tau^2\)|\(-\frac12\rho^2\tau^2\)|
|3|12|\(\frac{a+2b}{24}\tau\)|same to order \(\tau\)|\(\frac12\rho^2\tau^2\)|
|4|2|\(\frac18\)|same to order one|\(-\frac32\rho^2\tau^2\)|
|5|12|\(\frac{3a+2b}{72}\tau\)|same to order \(\tau\)|\(\frac16\rho^2\tau^2\)|
|6|6|\(\frac b9\tau\)|same to order \(\tau\)|\(\frac23\rho^2\tau^2\)|
|7|12|\(\frac1{18}\)|same to order one|\(-\frac13\rho^2\tau^2\)|
|8|6|\(\frac1{72}\)|same to order one|\(\frac16\rho^2\tau^2\)|

The feasibility inequality makes every displayed endpoint coefficient
positive.  Class 2 is the important rare class: its center mass and its
displacement are both of order \(\tau^2\), so a Taylor expansion merely to
second order in \(s\) is not uniform on this boundary scale.

## 4. Entropy expansion

Define

\[
X=b(2a+b),\qquad Y=X-12\rho^2.
\]

On every compact subset of \(a>b>2\sqrt3\rho>0\), the exact class formulas
give the uniform expansion

\[
\boxed{
\Delta(1-a\tau,1-b\tau,\rho\tau)
=\tau^2F(a,b,\rho)
+6\rho^2(a+2b)\tau^3\log\tau
+\tau^3J(a,b,\rho)
+O(\tau^4|\log\tau|),
}
\]

where \(J\) is an explicit analytic function on this chamber and

\[
\boxed{
F(a,b,\rho)
=\frac12\{X\log X-Y\log Y\}
+\rho^2\left[
\log\frac{3^6}{2^4(a+2b)^6(3a+2b)^2b^4}-6
\right].
}
\]

The symbolic script writes the complete, unsimplified exact formula for
\(J\) to `asymptotic_data.json`.

Here is the proof of the error order.  Each atom is a polynomial in \(\tau\).
After its leading power \(\tau^k\) is factored out, the remaining factor is
positive and analytic uniformly on compact subsets of the chamber.  Hence

\[
-p(\tau)\log p(\tau)
=-p(\tau)\{k\log\tau+\log[p(\tau)/\tau^k]\}
\]

has a convergent power-log expansion.  Summing the eight finite classes gives
the displayed remainder.  The apparent \(\tau^2\log\tau\) terms cancel:
class 2 contributes \(+12\rho^2\tau^2\log\tau\), while classes 3, 5 and 6
together contribute \(-12\rho^2\tau^2\log\tau\).  Collecting the remaining
order-\(\tau^2\) terms gives \(F\).  Direct symbolic expansion gives the next
logarithmic coefficient \(6\rho^2(a+2b)\).

## 5. Why splitting \(\alpha\) and \(\beta\) helps

Normalize by taking \(\tau=\beta\), and put

\[
x=\frac\alpha\beta>1,qquad r=\frac s\beta<\frac1{2\sqrt3}.
\]

Homogeneity gives \(F(a,b,\rho)=b^2G(x,r)\), where

\[
\begin{aligned}
G(x,r)
={}&\frac12\Big[(2x+1)\log(2x+1)
 -(2x+1-12r^2)\log(2x+1-12r^2)\Big]\\
&+r^2\left[
\log\frac{3^6}{2^4(x+2)^6(3x+2)^2}-6
\right].
\end{aligned}
\]

For small \(r\),

\[
\boxed{
G(x,r)=g_2(x)r^2-\frac{36}{2x+1}r^4+O(r^6),
}
\]

with

\[
g_2(x)=
\log\frac{3^6(2x+1)^6}{2^4(x+2)^6(3x+2)^2}.
\]

The positive quadratic term is the second-order entropy gain.  The negative
quartic term is exactly the first rare-event/Fisher penalty from class 2:

\[
-m_2\frac{\eta_2^2}{2p_2}
\sim -\frac{36}{2x+1}r^4\beta^2.
\]

Thus increasing the passive regularization \(\alpha\), while keeping the
moving margin \(\beta\) fixed, both improves the quadratic coefficient and
raises the class-2 baseline against the quartic penalty.

The old isotropic constraint is \(x=1\).  There

\[
g_2(1)=\log(729/400)=0.6002091849\ldots .
\]

Moreover

\[
g_2'(x)=
\frac{12(-x^2+2x+2)}{(2x+1)(x+2)(3x+2)},
\]

so the unique maximum for \(x>1\) occurs at

\[
x_*=1+\sqrt3=2.7320508075\ldots,
\qquad g_2(x_*)=1.0464962875\ldots .
\]

The split therefore increases the infinitesimal quadratic gain by about 74%,
while the magnitude of the quartic coefficient falls from \(12\) at \(x=1\)
to \(36/(3+2\sqrt3)=5.569\ldots\) at \(x=x_*\).

This proves a concrete positive-gap scale relation: choose any fixed \(x\)
near \(1+\sqrt3\), then choose a sufficiently small fixed \(r>0\).  Since
\(g_2(x)>0\), one has \(G(x,r)>0\), and therefore

\[
\alpha=x\beta,qquad s=r\beta,qquad
\Delta=G(x,r)\beta^2+O(\beta^3|\log\beta|)>0
\]

for all sufficiently small \(\beta\).  The improvement changes the leading
constant, not the exponent \(\beta^2\).

For fixed \(x\), the leading coefficient also has the exact derivative

\[
\frac{\partial G}{\partial r}
=2r\log\frac{3^6(2x+1-12r^2)^6}
{2^4(x+2)^6(3x+2)^2}.
\]

At the current \(x\), its interior zero is near \(r=0.2821\), close to the
spectral limit \(1/(2\sqrt3)\).  The finite maximizer lies farther inside
because the next boundary term is large and negative.

## 6. Check at the current finite maximizer

For

\[
(c,d,s)=(0.9849709285,0.9939746109,0.0011147170),
\]

we have

\[
\beta=0.0060253891,\qquad
x=2.4942906177,\qquad r=0.1850033220.
\]

The exact eight-class formula gives

\[
\Delta=2.747238545728012\ldots\times10^{-7}.
\]

At these ratios,

\[
G(x,r)=0.02846834752\ldots,
\qquad \beta^2G=1.0335522902\ldots\times10^{-6}.
\]

The leading term alone overestimates the finite gap because \(\beta\) is not
yet in the uniform leading regime.  The next coefficients are

\[
6r^2(x+2)=0.9229357236\ldots,
\qquad J(x,1,r)=1.2000506116\ldots .
\]

Numerically,

\[
\begin{aligned}
\beta^3[6r^2(x+2)\log\beta]
  &=-1.0320442707\ldots\times10^{-6},\\
\beta^3J(x,1,r)
  &= 2.6251544196\ldots\times10^{-7}.
\end{aligned}
\]

Thus the expansion through order \(\beta^3\) gives

\[
2.6402346150\ldots\times10^{-7},
\]

within about 4% of the exact value.  This is also why the observed finite
ratio \(r=0.1850\) is much smaller than the leading-order optimizer: the
negative \(\beta^3\log\beta\) correction nearly cancels the leading gain at
this value of \(\beta\).

## 7. Reproduction

From this directory run

```text
uv run --with sympy python -B derive_asymptotics.py
```

The script reconstructs all 64 determinants, proves the eight polynomial
classes by exact equality, checks both normalization identities, derives the
power-log series, verifies symbolically that its \(\tau^2\) coefficient equals
the displayed \(F\), and reproduces the finite-point comparison.
