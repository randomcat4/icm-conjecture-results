# I05: exact counterexample to Lyons's Conjecture 2.6

This package gives an exact counterexample to the original finite-dimensional complex-Hilbert-space form of Lyons's Conjecture 2.6. For the explicit matrices in the paper, (K\pm hA) are strict positive contractions on five fixed labelled coordinates, their midpoint is (K), and

\[
H(K)<\frac{H(K-hA)+H(K+hA)}2.
\]

Here (K) is real symmetric, (A=iB) is purely imaginary Hermitian of rank four, and (h=1/100000). All 32 atom probabilities are exact positive rationals, each distribution sums to one, and the two endpoint distributions agree atom by atom. Directed rational logarithm bounds certify both the strict finite-chord entropy violation and a strictly positive Hessian direction.

For a detailed Chinese account of the conjecture's 2003 origin, its relation to the Lyons--Steif stationary conjecture, its 2014 ICM restatement, later partial results, author profiles, and the exact scope of this project, see [history_zh.md](history_zh.md). A compact index of primary sources with direct PDF links is in [sources/README.md](sources/README.md).

The [projection-boundary extension](projection_boundary_extension/README.md) supplies a second,
local-curvature proof mechanism together with broader structural results: real concavity in inward
neighbourhoods of connected projections, an exact one-parameter real/complex separation family,
and a mixed fourth-order obstruction for independent real blocks.

For initial review by the conjecture's author, the
[author-review package](author_review_note/README.md) gives a shorter five-page proof organized
around the information-geometric mechanism, together with its own dependency-free exact checker
and a brief email draft.

The [larger-gap package](larger_gap/README.md) gives two further exact finite chords: an optimized
five-point example with gap greater than `1/130000000`, and a symmetric six-point example with
gap greater than `1/40000000`.  The latter has a short block construction and only eight distinct
center/endpoint probability pairs.

## Reproduce the exact check

Only Python's standard library is required. From this directory, run:

```console
python scripts/verify_exact.py . --output reproduction/exact_check.json
```

The command reconstructs the matrices and all event polynomials with rational arithmetic. It checks symmetry and skew-symmetry, hence the Hermitian property of (A=iB); proves rank four; verifies strict feasibility of (K) and both endpoints; checks all 32 center events and both 32-event endpoints; checks exact normalization and endpoint equality; and verifies the directed logarithm enclosures for the entropy gap and Hessian direction. A successful checked run is recorded in [`reproduction/exact_check.json`](reproduction/exact_check.json).

To rebuild the paper with shell escape disabled:

```console
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
```

The checked build outcome is recorded in [`reproduction/build_receipt.md`](reproduction/build_receipt.md).

## Original source and scope

Russell Lyons stated Conjecture 2.6 in *Determinantal Probability: Basic Properties and Conjectures*, Proceedings of the International Congress of Mathematicians 2014, volume IV, pages 137--161. The coordinate DPP definition begins on printed page 137; the real-or-complex Hilbert-space setting appears on page 138; positive contractions are defined on page 140; and the full subset-distribution entropy and Conjecture 2.6 appear on page 141. See the [published PDF](https://rdlyons.pages.iu.edu/pdf/icm-pub.pdf) and [arXiv:1406.2707](https://arxiv.org/abs/1406.2707).

The paper source, bibliography, PDF, and certificate tables are byte-for-byte copies of the frozen manuscript at commit `4a58b2b9c6a9a20c18f3f330081b5d3bc8ac3f3e`. The public checker is a byte-for-byte copy of the final fresh independent paper verifier at commit `832cbf1d75c6529da25813f0021f8ba5d822a7f3`. The underlying exact construction was frozen at `9c27d3f25cb68ae8fe1d6596ffbb49d023d45892` and separately reconstructed by two independent mathematical checks frozen at `e4288b2886093f0b4e07675bdf83f5c0d8c7c141` and `eb6f5aa54e8c462a48cbae06069852fbc95f4d8a`. The concise [verification summary](verification_summary.md) records what those checks establish without publishing review prompts, internal reports, or unrelated research data.

The frozen manuscript still contains `Author placeholder`, `Internal draft`, and a paragraph describing the state before fresh paper review. Those strings are retained because the verified mathematical manuscript is copied without alteration. The later fresh review returned `CORRECT`; this public README records that outcome separately. Authorship and any submission decision remain to be decided.

## Limits

This package does not resolve the restriction to real symmetric kernels, prove that dimension five is minimal, establish absolute novelty or priority, provide a formal proof-assistant development, decide authorship, or imply submission to any journal, proceedings, preprint server, or repository release service. The bounded source review used during verification found no earlier result in its search, but that is not an exhaustive priority determination.

No license or DOI is asserted by this package.

## Files

- [history_zh.md](history_zh.md): detailed Chinese history, author context, and a precise account of this project's contribution.
- [sources/](sources/): primary-source index, direct external PDF links, and a reusable BibTeX bibliography.

- [`main.pdf`](main.pdf), [`main.tex`](main.tex), and [`references.bib`](references.bib): frozen six-page manuscript and source.
- [`certificates/`](certificates/): exact JSON, TSV, and TeX certificate data for all 32 events.
- [`scripts/verify_exact.py`](scripts/verify_exact.py): standard-library exact checker.
- [`verification_summary.md`](verification_summary.md): public-safe verification coverage and limitations.
- [`reproduction/`](reproduction/): outputs from the public clean run and paper build.
- [`projection_boundary_extension/`](projection_boundary_extension/): the generalized
  projection-boundary proof, exact certificates, and a careful account of the quantities
  \(\Lambda_P(D)\) and \(\ell(P,N,B)\).
- [`author_review_note/`](author_review_note/): five-page human-readable note, exact checker,
  boundary data, and an unsent email draft.
- [`larger_gap/`](larger_gap/): independently verified five-point and symmetric six-point
  counterexamples with substantially larger finite entropy gaps.
