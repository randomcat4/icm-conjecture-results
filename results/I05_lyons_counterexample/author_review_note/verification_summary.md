# Verification summary

The five-page author-review note was checked in fresh, read-only review contexts after the mathematical files were frozen.  The added real-block certificate received a separate focused review after the final revision.

## Frozen files

- `main.tex`: `EBDFCD0CB12D521DA2771C37415F9CE07830E60F838421065D4F1AF7A7710523`
- `verify_exact.py`: `FBEA3920226D096A9490DE89EFFD6A8BA39C5AC68EDCA9DE251E445E5520148C`
- `main.pdf`: `6FEE4EC5C23BCF8CBDECBB2B0E5A3C9C85056F6C172595BF5FC2BFFCCC68B961`
- `author_review_email.md`: `E3933E78B466EA5177CAD0874704DBDCDCF27F0E677FD3D14FFC6DB01D62610E`

## Mathematical verdict

**CORRECT.** The reviewer independently checked:

- strict feasibility of the explicit endpoints \(Q_\pm=K\pm10^{-9}iB\);
- exact eventwise equality of the endpoint laws;
- the entropy Hessian identity and the disappearance of its Fisher term;
- the fact that \(PB=BP\), rather than the coefficient \(-1/2\), removes the boundary logarithmic divergence;
- the \(U_{2,5}\) boundary orders and leading coefficients;
- the exact prime-factor form and positivity of the boundary limit;
- every direction choice in the rational logarithm intervals;
- the strict finite entropy-gap enclosure.

The independently inspected directed interval is

\[
2.8151854998206740972\times10^{-19}
<
H(Q_+)-H(K)
<
2.8151976274156489613\times10^{-19}.
\]

Its width is about \(1.21\times10^{-24}\), so the positive sign is not hidden by the logarithm-enclosure error.

## Real-symmetric block verdict

**CORRECT.** A separate reviewer checked all 15 real symmetric directions and all 225 mixed second derivatives.  The rational preconditioner is upper triangular with nonzero diagonal, hence invertible.  After the associated congruence, the interval enclosure for the negative Hessian is strictly diagonally dominant, with an exact positive Gershgorin margin whose decimal value is approximately \(0.9999634098734166\).  Therefore the Hessian at the displayed midpoint is negative definite on the entire real symmetric block.  This is a local statement at that kernel and does not settle the real-symmetric restriction in general.

## Exposition verdict

**CORRECT.** The second reviewer found that the note and email use normal mathematical language, contain no repository link or Gu 2020 reference, distinguish Conjecture 2.6 as stated from the unresolved real-symmetric restriction, and do not overstate the corrected boundary mechanism. The forbidden `operatorname` macro is absent.

These are internal independent checks. They are not external peer review, a priority determination, or an endorsement by Russell Lyons.
