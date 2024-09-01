# TODO

General

- Warn on `0` tags
- Why does `Node` call `exit()` on failure in `getEigenvector`?
- Some methods are preemptively marked `virtual`; I've changed
  some of these to use a macro, `VIRTUAL`, which is defined empty.

Organizing:

- `material/damage/` <- `SRC/damage`
- `transform` <- coordTransformation
- `container` <- tagged/
- `analysis/convergence` <- convergenceTest/
- `builder`   <- modelbuilder/
- `streams`   <- handler/

- `analysis/algebra/` (<- SRC/system_of_eqn)
  - `linear/`
  - `eigen/`

- `algebra`   <- `matrix/`
  - `inverse/`
  - `rotation/`
  - `euclid/`
  - `Vector.h`
  - `Matrix.h`


Commands:

- `converge` command
  # analysis.cpp (config)
  strategy test <convergence>
  strategy iter <algorithm>
  strategy incr <integrator>

  # solve.cpp
  solve <step> <n>
  solve  iter  <n>  ; # Iterate
  solve  tang       ; # Invert
  solve  test       ; # Converge
  solve  incr  <dt> ; # Increment

  # status.cpp
  stat norm ; # testNorm
  stat iter ; # testIter


