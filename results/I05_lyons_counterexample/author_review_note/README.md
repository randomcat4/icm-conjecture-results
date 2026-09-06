# Author review package

This directory contains the short, human-readable version intended for initial review by Russell Lyons.

- [`main.pdf`](main.pdf): the five-page review note.
- [`main.tex`](main.tex): LaTeX source.
- [`author_review_email.md`](author_review_email.md): a brief email draft with no repository link; it has not been sent.
- [`verify_exact.py`](verify_exact.py): a self-contained standard-library checker for the boundary limit, imaginary-direction Hessian, negative definiteness of the entire real symmetric Hessian block, endpoint feasibility, eventwise endpoint equality, and finite midpoint entropy gap.
- [`boundary_table.csv`](boundary_table.csv): an optional readable copy of the 32 boundary terms; the checker does not read it.
- [`verification_summary.md`](verification_summary.md): fresh read-only review verdicts and frozen file hashes.

The note gives the information-geometric mechanism, the four-vector construction, the exact positive boundary limit, and an explicit midpoint violation at \(\varepsilon=10^{-6}\). It is designed to be sent as a PDF together with the single checker file, without exposing the repository.
