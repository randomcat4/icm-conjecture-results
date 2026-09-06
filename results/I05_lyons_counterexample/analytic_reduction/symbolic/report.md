# Symbolic route report

## Verdict

**PROVED for the assigned symbolic target.** The 64 exact-event determinants
reduce identically to eight polynomial classes with multiplicities

\[
2,12,12,2,12,6,12,6.
\]

Closed formulas for all center and endpoint probabilities are in
`class_formulas.md`. The three exact interior stationary equations are in
`stationary_equations.md`.

This route does **not** prove local or global uniqueness of the observed
stationary point. It reduces that question to a three-variable system of
seven-term logarithmic sums suitable for a separate interval analysis.

## Main structural gain

With $z=s^2$ and probabilities scaled by $576$, every endpoint class has the
form

\[
Q_j=P_j+zA_j+z^2B_j.
\]

Thus the apparent 64-determinant optimization is exactly a three-variable
entropy built from eight integer polynomials. Normalization removes the
constant terms after differentiation, and the stationary equations become
logarithms of ratios $P_j/P_8$ and $Q_j/Q_8$.

## Exhaustion proof

The order-12 coordinate group preserving $K$ and preserving $B$ up to sign,
together with particle-hole complementation, has exactly eight subset orbits.
They are classified by subset size, block occupancy, and whether cross-block
entries are of $\alpha=(c+d)/3$ or $\beta=(c-2d)/3$ type. Their sizes sum to
64. The script checks exact symbolic determinant equality on each orbit and
distinct generic polynomial pairs between the eight orbits.

## Reproducibility

Run `derive_classes.py` with Python and SymPy. The script:

1. constructs $C,K,B$ symbolically;
2. verifies $CS_0=-S_0C$ and $KB=BK$;
3. expands all 64 center and endpoint signed determinants exactly;
4. groups them by exact polynomial equality;
5. verifies the multiplicities, orbit partition, evenness in $s$, the
   human-readable formulas, and both normalization identities;
6. cross-checks all 128 values at three rational points using an independent
   Leibniz determinant over `Fraction` pairs.

The rational-point checks guard against implementation errors but are not used
to infer polynomial identities.

## Numerical orientation only

Solving the three derived equations at high precision from the recorded
candidate refines the stationary point to

$$(c,d,s)=(0.98497066858425021191,
0.99397451323581841410,
0.0011147211944450686890).$$

with gap approximately
$2.74723854723614806595\times10^{-7}$.  The numerical Hessian eigenvalues are
approximately $-1.83589488,-0.04786910,-0.00308570$.  This is evidence for a
strict local maximum, not an interval proof.

## Remaining analytic obstacle

The logarithms keep the stationary system transcendental. A rigorous
uniqueness result will require interval bounds for the three log-ratio
equations and their Jacobian on a stated box, plus exclusion of the rest of the
chamber. Repeated numerical convergence cannot supply that step.
