# Exact stationary equations

**Status: PROVED as an exact reduction.** These equations are necessary for
every interior stationary point. They do not prove uniqueness or global
maximality of the observed solution.

Use $P_j,Q_j=P_j+D_j$ and $m_j$ from `class_formulas.md`, and put $z=s^2$.
Because $s>0$ in the chamber, stationarity in $s$ is equivalent to
stationarity in $z$.

Let $h(x)=-x\log x$. Normalization gives

\[
\Delta(c,d,z)
=\sum_{j=1}^8m_j\left[h\left(\frac{Q_j}{576}\right)
                         -h\left(\frac{P_j}{576}\right)\right]
=\frac1{576}\sum_{j=1}^8m_j(P_j\log P_j-Q_j\log Q_j).
\]

The terms involving $\log576$ cancel because
$\sum m_jP_j=\sum m_jQ_j=576$.

For $\theta=c,d$, differentiated normalization gives

\[
576\,\partial_\theta\Delta
=\sum_{j=1}^8m_j
 \{P_{j,\theta}\log P_j-Q_{j,\theta}\log Q_j\}.
\]

Since the center probabilities do not move with $z$,

\[
576\,\partial_z\Delta
=-\sum_{j=1}^8m_jQ_{j,z}\log Q_j.
\]

Use class 8 as a reference to remove one redundant term and every dimensional
logarithm. The exact stationary equations are

\[
\boxed{
F_c:=\sum_{j=1}^{7}m_j
\left[P_{j,c}\log\frac{P_j}{P_8}
-Q_{j,c}\log\frac{Q_j}{Q_8}\right]=0,}
\]

\[
\boxed{
F_d:=\sum_{j=1}^{7}m_j
\left[P_{j,d}\log\frac{P_j}{P_8}
-Q_{j,d}\log\frac{Q_j}{Q_8}\right]=0,}
\]

\[
\boxed{
F_z:=\sum_{j=1}^{7}m_jQ_{j,z}
\log\frac{Q_j}{Q_8}=0.}
\]

These are logarithms of rational functions with integer-polynomial
numerators and denominators. There is no hidden determinant or 64-event sum.

For explicit coefficients, write

\[
D_j=z(A_j+zB_j).
\]

Then

\[
\begin{array}{c|l|l}
j&A_j&B_j\\ \hline
1&-216(1-c^2)(1+d^2)&1296(1-c^2)\\
2&72\{d^2-1-c^2(3d^2+1)\}&432(3c^2-1)\\
3&72(3c^2d^2+c^2+d^2-1)&-432(3c^2+1)\\
4&-216(1+c^2)(1+d^2)&1296(1+c^2)\\
5&24(9c^2d^2-3c^2-4cd-d^2+3)&144(1-9c^2)\\
6&24(9c^2d^2-3c^2+8cd-d^2+3)&144(1-9c^2)\\
7&24(-9c^2d^2+3c^2-4cd-d^2+3)&144(9c^2+1)\\
8&24(-9c^2d^2+3c^2+8cd-d^2+3)&144(9c^2+1).
\end{array}
\]

Thus all coefficients in the boxed equations are explicit:

\[
Q_{j,z}=A_j+2zB_j,
\]

\[
Q_{j,c}=P_{j,c}+zA_{j,c}+z^2B_{j,c},
\qquad
Q_{j,d}=P_{j,d}+zA_{j,d}+z^2B_{j,d}.
\]

Together with the eight displayed $P_j,A_j,B_j$, this is a complete explicit
system of three logarithmic equations in $c,d,z$. The factored form retains
the boundary scales $1-c$, $1-d$, and $z$ and is easier to audit than a full
coefficient expansion.

The equations apply where all $P_j,Q_j>0$, in particular in the stated
strict-contraction chamber. If an atom vanishes on the boundary, one-sided
limits are required; the interior equations alone do not treat that case.
