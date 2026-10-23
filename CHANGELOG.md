# Changelog

## [0.3.0] - 2026-10-21

### Breaking: grids that `f64` cannot resolve are rejected

`DiscretizedGrid` stores grid points as `f64` (`lower_bound + g * step`), so a
resolution is only usable while neighbouring points stay distinguishable in
`f64`. That limit depends on the bound magnitudes, not on the resolution alone:
`[0, 1)` supports `R = 50`, symmetric `O(1)` bounds `R = 51`, and
`[1e10, 1e10 + 1)` only `R = 17`. Beyond it the coordinate round trip silently
returned a neighbouring index (off by one at `R = 52`, off by hundreds at
`R = 62`) — see tensor4all/quanticsgrids-rs#4.

- `DiscretizedGridBuilder::build()` now returns the new
  `QuanticsGridError::GridNotRepresentable` for such grids instead of returning
  a grid whose conversions are wrong. The check is conservative — it requires
  the step to be at least four ulps of the bounds and checks that the last two
  grid points stay inside the domain and map back to their own index — so a grid
  whose arithmetic happens to round nowhere is rejected too once its step falls
  below that limit. **Migration:** reduce the resolution; the error names the
  dimension, step and bounds.
- The grid step is now computed and validated once in `build()` from the
  declared bounds instead of being re-derived from a rounded internal upper
  bound. This also keeps `include_endpoint(true)` coordinates on their endpoint:
  the last coordinate is exactly `upper_bound`, where it used to miss it by an
  ulp (and could land above it, so that `origcoord_to_grididx` rejected the
  grid's own output).
- `origcoord_to_grididx` documents the rounding it implements:
  `round((c - lower_bound) / grid_step)`, clamped to `0..=N-1`, ties to the
  larger index, and no wrap-around in the half-open case.
- `origcoord_to_grididx`/`grididx_to_origcoord`/`grid_step`/`grid_max` and
  `grid_origcoords` are unchanged otherwise.

`QuanticsGridError` gained a variant, which breaks exhaustive `match`es.

### Tests

Index-conversion tests at `R = 58`/`R = 62` moved to `InherentDiscreteGrid`,
where coordinates are not involved; `DiscretizedGrid` now has a round-trip test
over every accepted resolution and bound regime.

## [0.2.0] - 2026-08-12

### Breaking: 0-indexed grid indices and quantics digits

Migrates the public API from the 1-indexed (Julia) convention to 0-indexed,
matching the rest of the Rust workspace. This removes the only 1-indexed
surface in the tensor4all-rs stack (tensor4all/tensor4all-rs#584).

**Porting note for QuanticsGrids.jl / QuanticsTCI.jl scripts: subtract 1 from
grid indices and quantics digits at the call boundary.**

- Grid indices are now `0..N` (first grid point is `[0, 0, ...]`) instead of
  `1..=N`.
- Quantics digits are now `0`/`1` bit values (and `0..base-1` digits for
  non-binary bases) instead of `1`/`2`.
- Index-table bit numbers (`QuanticsIndex`) are 0-indexed (0 = most
  significant bit) instead of 1-indexed.
- Grid-index and quantics parameters/results are typed `usize` instead of
  `i64`, making the convention visible in the type. Original coordinates
  (origin/step/origcoord) remain `i64` (InherentDiscreteGrid) or `f64`
  (DiscretizedGrid) and are unchanged.
- `max_grididx()` now returns the maximum 0-based index (`base^R - 1`)
  instead of the point count (`base^R`).
- Conversion formulas: `origcoord = origin + grididx * step` and
  `grididx = (coord - origin) / step` (no `±1`).
