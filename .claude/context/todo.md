# TODO

- [ ] Follow-on: GPU batched evaluator (`strains_stresses_from_B_gpu`, spec §5b). Needs
      `[weakdeps] CUDA` + `ext/GXBeamCSCUDAExt.jl` and a `julia` compat bump from 1.7 to 1.9.
      Not attempted yet — no CUDA-capable machine available to verify it.
- [ ] Regularized `cc.F` solve (spec §6) if the near-singular solve (`det(cc.F) ~ 1e-9` on
      real turbine cross-sections) ever throws `SingularException` in practice. Marked with
      a `# TODO(near-singular cc.F)` comment at the one call site in
      `_assemble_operators` (src/clt.jl).
- [ ] Pre-existing, unrelated test failure: `test/meshtools.jl:89`, `"Mesh Translation"`
      testset, `find_direction` returns NaN on one case. Confirmed present on a clean HEAD
      checkout before this session's changes — not something introduced by the operator work.
- [ ] Decide whether to revive the legacy commented-out `test/clt.jl` / `test/fem.jl` suites
      (currently excluded from `runtests.jl`) — left alone this session per explicit scope
      decision.
