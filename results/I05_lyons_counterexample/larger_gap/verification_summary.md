# Verification summary

Both candidates were checked in fresh read-only contexts using a second calculation based on exact principal minors, inclusion-exclusion, rational Gaussian elimination, and independent directed logarithm bounds.

## Five-point candidate

**CORRECT.** All four Hermitian matrices `Q_+`, `Q_-`, `I-Q_+`, and `I-Q_-` have positive Sylvester minors.  All 96 center and endpoint event probabilities are strictly positive and exactly normalized, and the two endpoint vectors agree entry by entry.  The independently recomputed interval proves

\[
\frac1{130\,000\,000}<H(Q_+)-H(K)<\frac1{120\,000\,000}.
\]

Its interval width is approximately `7.22e-46`.

## Symmetric six-point candidate

**CORRECT.** The independent calculation verified the identities defining `P` and `B`, strict endpoint feasibility, all 64 positive normalized event probabilities, endpoint equality, and

\[
\frac1{40\,000\,000}<H(Q_+)-H(K)<\frac1{39\,000\,000}.
\]

Its interval width is below `1.453e-49`.

## Frozen checker hashes

- `certify_optimized_five.py`: `0D7D113D882891122502A00EA0874E12F87F772B7712744E6928F9CB10D8FCC1`
- `certify_symmetric_six.py`: `2BBD28941ED1B4FA275F2F05271F15806F90BF3F087767BFCD856EDDF2AA2FFA`

These are exact finite certificates.  The reviews do not claim dimension minimality, a largest possible gap in fixed dimension, or a result for the separate all-real restriction.
