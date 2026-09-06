# Uniqueness-route report

## Verdict

`PROVED_BOX`

The observed numerical maximizer has been upgraded to a rigorous local
object.  Arb interval arithmetic proves:

- exactly one stationary point in the explicit Krawczyk box `X_0`;
- negative-definite Hessian throughout `X_0`;
- negative-definite Hessian throughout a much larger explicit box `X`;
- hence a unique stationary point and unique global maximizer on `X`.

The proof does not exclude a larger value elsewhere in the full feasible
chamber.

## Reproduction

From this directory run:

```powershell
uv run --offline --with sympy --with python-flint python certify_local.py
uv run --offline --with sympy --with python-flint python certify_box.py
```

The first command should finish with `PROVED_LOCAL`.  The second reruns the
local certificate, checks the adaptive 128-subbox cover, and finishes with
`PROVED_BOX`.

`formulae.json` contains the eight exact center/endpoint probability classes.
The sibling program `../main/derive.py` independently derives those formulas
from the six-by-six determinants and verifies their multiplicities and total
mass.

## Files

- `certify_local.py`: Krawczyk and local Hessian certificate.
- `certify_box.py`: strict-concavity cover of the larger box.
- `formulae.json`: exact probability polynomials and orbit data.
- `local_certificate.txt`, `box_certificate.txt`: captured successful runs.
- `box_definition.md`: exact dyadic definitions of both boxes.
- `proof_or_gap.md`: proof explanation and precise global gap.
- `explore_stationary.py`: non-rigorous high-precision exploration, retained
  only for provenance; it is not used by either certificate.

No floating-point comparison is used in any proof decision.  Floating-point
conversion in `certify_box.py` is used only to select which already certified
lower endpoint to print as the worst one.
