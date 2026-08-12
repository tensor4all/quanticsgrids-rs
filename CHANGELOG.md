# Changelog

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
