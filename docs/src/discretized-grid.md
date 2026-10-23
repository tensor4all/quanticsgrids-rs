# Discretized Grid

`DiscretizedGrid` discretizes a continuous `d`-dimensional domain and supports
conversions between:

- quantics indices
- grid indices
- original coordinates

Unlike the Julia package, the Rust API returns an error when original
coordinates fall outside the configured domain.

## Creating a one-dimensional grid

We can discretize the interval `[0, 1)` with `R` bits as follows:

```rust
{{#include ../examples/discretized-grid-1d.rs}}
```

## Supported resolutions

Coordinates are `f64`, so grid points must stay distinguishable: the rounding
error of a conversion is about one ulp of the bounds, and a grid whose step falls
below that limit cannot map coordinates to indices exactly. The limit therefore
depends on the bounds, not on `R` alone:

| Bounds | Highest supported `R` |
|---|---|
| `[0, 1)` | 50 |
| symmetric `O(1)`, e.g. `[-1.5, 1.5)` | 51 |
| `[0.3, 10.3)` | 50 |
| `[1e10, 1e10 + 1)` | 17 |
| `[1e12, 1e12 + 1e6)` | 30 |

`build()` returns `QuanticsGridError::GridNotRepresentable` above these, instead
of producing a grid whose conversions silently return a neighbouring index. The
check is conservative (step of at least four ulps of the bounds, plus the last
two grid points checked), so it can also reject a grid that would happen to be
exact. Within the accepted range,
`origcoord_to_grididx(&grididx_to_origcoord(&[g])?)` returns `g` exactly.

## Mapping coordinates to indices

A coordinate `c` maps to the index `round((c - lower_bound) / grid_step)`, clamped
to `0..=N-1`; a normalized index exactly halfway between two indices rounds to
the larger one. The rounding applies to that computed ratio, so the result is not
promised to be the geometrically nearest emitted coordinate in every case.

The grid does not wrap around: in the half-open case a coordinate in the last step
below the upper bound keeps the last grid point. Coordinates outside the domain
(`[lower, upper)`, or `[lower, upper]` with `include_endpoint(true)`) return
`QuanticsGridError::CoordinateOutOfBounds`.

The Rust builder also supports including the upper endpoint:

```rust
{{#include ../examples/discretized-grid-1d-endpoint.rs}}
```

## Creating a `d`-dimensional grid

A `d`-dimensional grid is created in the same way by supplying one resolution
per dimension. You can choose either fused or interleaved unfolding.

### Fused representation

```rust
{{#include ../examples/discretized-grid-3d-fused.rs}}
```

Incrementing the least significant fused index increments the `x` coordinate
first.

### Interleaved representation

```rust
{{#include ../examples/discretized-grid-3d-interleaved.rs}}
```
