# Author review note, version 2

This directory contains the six-point rewrite intended for mathematical review.

- `main.pdf`: compiled review note.
- `main.tex`: LaTeX source.
- `main_zh.pdf`: Chinese version of the review note.
- `main_zh.tex`: Chinese LaTeX source.
- `verify_exact.py`: standard-library checker for endpoint feasibility, all 64 complete-event probabilities, endpoint-law equality, and the exact entropy-gap enclosure.
- `references.bib`: the two original sources for the conjecture.
- `verification_summary.md`: compilation, visual, symbol, and exact-arithmetic checks.

The proof uses the symmetric six-point construction.  It introduces only one auxiliary scalar, `u = epsilon(1-epsilon)`, to shorten the eight-row event table.  The information-geometric mechanism is explained after the direct proof rather than used as a prerequisite.
