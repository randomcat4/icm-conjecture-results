# Verification summary

The five-page author-review note was checked in two fresh, read-only review contexts after the mathematical files were frozen.

## Frozen files

- `main.tex`: `73A86A4E620B9E3AB7E937D0EBD2218B758A428E1CED3C2276A81D62EC7C973B`
- `verify_exact.py`: `BDD005F16D9F5E1D84E3ADC4491757E52EFCB3E751DEF261B049380C3687DDD8`
- `main.pdf`: `D7EF292A5B77550C264998BD9CD6C4B8C6CA28A7B0CE62DF5F674EDFF40D1CD0`

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

## Exposition verdict

**CORRECT.** The second reviewer found that the note and email use normal mathematical language, contain no repository link or Gu 2020 reference, distinguish Conjecture 2.6 as stated from the unresolved real-symmetric restriction, and do not overstate the corrected boundary mechanism. The forbidden `operatorname` macro is absent.

These are internal independent checks. They are not external peer review, a priority determination, or an endorsement by Russell Lyons.
