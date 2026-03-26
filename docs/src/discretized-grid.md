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
