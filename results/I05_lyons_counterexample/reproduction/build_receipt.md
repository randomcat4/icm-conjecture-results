# Clean paper build receipt

Build date: 2026-09-05.

The build used a clean directory containing only `main.tex`, `references.bib`, and `certificates/event_table.tex`. Shell escape was disabled. The prescribed recipe completed as follows:

| Step | Program | Exit status |
|---:|---|---:|
| 1 | `pdflatex -no-shell-escape` | 0 |
| 2 | `bibtex` | 0 |
| 3 | `pdflatex -no-shell-escape` | 0 |
| 4 | `pdflatex -no-shell-escape` | 0 |

Tool versions:

- MiKTeX-pdfTeX 4.27 (MiKTeX 26.5), pdfTeX engine 1.40.29.
- MiKTeX-BibTeX 4.2 (MiKTeX 26.5), BibTeX engine 0.99e.

The rebuilt artifact has 6 A4 pages and is 272,766 bytes. Its extracted-text SHA-256 is `EFF07BD2BA1C8267DDC114824B4A8A263AA9896DE7484F7E34756425317A453F`, identical to the extracted-text hash of the frozen `main.pdf` in this package.

The final log contains zero undefined citations, zero undefined references, zero overfull boxes, zero underfull boxes, and zero LaTeX warnings. The `epstopdf` package reports one package warning that shell escape is not enabled; this is expected and confirms the requested build setting. MiKTeX also reports its local update-status notice.
