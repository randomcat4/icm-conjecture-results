# Verification summary for version 2

The English and Chinese notes each compile to four A4 pages. Their final
LaTeX logs contain no errors, undefined references, overfull boxes, underfull
boxes, or LaTeX warnings. All four pages of both PDFs were rendered and
inspected after the final compilation.

The standard-library checker proves:

- both endpoints are strict positive contractions;
- the endpoint laws agree for all 64 complete events;
- the 64 center/endpoint pairs agree with the eight rows and representatives
  printed in the note;
- the event probabilities are positive and sum to one;
- the directed rational entropy interval lies strictly between
  1/40000000 and 1/39000000;
- the midpoint entropy lies strictly between
  2.877726131217012 and 2.877726131217013.

The sources contain no custom LaTeX macros and do not use the forbidden
`operatorname` macro. The only auxiliary mathematical abbreviation in the
main construction is `u = epsilon(1-epsilon)` in the event table.

SHA-256:

- main.pdf: F5FA161A9995FDD4F408F68951494974B08E33DA731E930C0907481F9F9881CC
- main.tex: B4A8DAEBE55D31F875CFB7979D8C14BFDD7B726CA155551275C30378312F946C
- main_zh.pdf: EE5FA0031D25ABA15631664A8C93680F181E1E2E5CA0DFAD995B4B174D306EE4
- main_zh.tex: 5EF0EDE341F78DEC09606E9915D5AA695DCA675C653ED8DC64715CDEC569C4E4
- verify_exact.py: 060B6881D9E0DECD0BA6CE67D42128670D2963F3386F202B9A603D6DA259782F
