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
# extern crate quanticsgrids;
use quanticsgrids::DiscretizedGrid;

fn main() -> quanticsgrids::Result<()> {
    let r = 4;
    let grid = DiscretizedGrid::builder(&[r])
        .with_bounds(0.0, 1.0)
        .build()?;

    let quantics = vec![1; r];
    let grididx = vec![1];

    assert_eq!(grid.quantics_to_grididx(&quantics)?, grididx);
    assert_eq!(grid.grididx_to_quantics(&grididx)?, quantics);
    assert!((grid.grididx_to_origcoord(&grididx)?[0] - 0.0).abs() < 1e-12);
    assert_eq!(grid.origcoord_to_grididx(&[0.0])?, grididx);

    let quantics = vec![2; r];
    let grididx = vec![2_i64.pow(r as u32)];
    let x = 1.0 - 1.0 / 2.0_f64.powi(r as i32);

    assert_eq!(grid.quantics_to_grididx(&quantics)?, grididx);
    assert_eq!(grid.grididx_to_quantics(&grididx)?, quantics);
    assert!((grid.quantics_to_origcoord(&quantics)?[0] - x).abs() < 1e-12);
    assert!((grid.grididx_to_origcoord(&grididx)?[0] - x).abs() < 1e-12);
    assert_eq!(grid.origcoord_to_grididx(&[x])?, grididx);
    assert_eq!(grid.origcoord_to_quantics(&[x])?, quantics);

    Ok(())
}
```

The Rust builder also supports including the upper endpoint:

```rust
# extern crate quanticsgrids;
use quanticsgrids::DiscretizedGrid;

fn main() -> quanticsgrids::Result<()> {
    let r = 4;
    let grid = DiscretizedGrid::builder(&[r])
        .with_bounds(0.0, 1.0)
        .include_endpoint(true)
        .build()?;

    assert!((grid.grididx_to_origcoord(&[1])?[0] - 0.0).abs() < 1e-12);
    assert!((grid.grididx_to_origcoord(&[2_i64.pow(r as u32)])?[0] - 1.0).abs() < 1e-12);

    Ok(())
}
```

## Creating a `d`-dimensional grid

A `d`-dimensional grid is created in the same way by supplying one resolution
per dimension. You can choose either fused or interleaved unfolding.

### Fused representation

```rust
# extern crate quanticsgrids;
use quanticsgrids::{DiscretizedGrid, UnfoldingScheme};

fn main() -> quanticsgrids::Result<()> {
    let r = 4;
    let grid = DiscretizedGrid::builder(&[r, r, r])
        .with_lower_bound(&[0.0, 0.0, 0.0])
        .with_upper_bound(&[1.0, 1.0, 1.0])
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()?;

    let quantics = vec![1; r];
    let grididx = vec![1, 1, 1];

    assert_eq!(grid.quantics_to_grididx(&quantics)?, grididx);
    assert_eq!(grid.grididx_to_quantics(&grididx)?, quantics);
    assert_eq!(grid.origcoord_to_grididx(&[0.0, 0.0, 0.0])?, grididx);

    let quantics = vec![1, 1, 1, 2];
    let grididx = vec![2, 1, 1];

    assert_eq!(grid.quantics_to_grididx(&quantics)?, grididx);
    assert_eq!(grid.grididx_to_quantics(&grididx)?, quantics);

    let coord = grid.quantics_to_origcoord(&quantics)?;
    assert!((coord[0] - 1.0 / 16.0).abs() < 1e-12);
    assert!((coord[1] - 0.0).abs() < 1e-12);
    assert!((coord[2] - 0.0).abs() < 1e-12);

    Ok(())
}
```

Incrementing the least significant fused index increments the `x` coordinate
first.

### Interleaved representation

```rust
# extern crate quanticsgrids;
use quanticsgrids::{DiscretizedGrid, UnfoldingScheme};

fn main() -> quanticsgrids::Result<()> {
    let r = 4;
    let grid = DiscretizedGrid::builder(&[r, r, r])
        .with_lower_bound(&[0.0, 0.0, 0.0])
        .with_upper_bound(&[1.0, 1.0, 1.0])
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()?;

    let quantics = vec![1; 3 * r];
    let grididx = vec![1, 1, 1];

    assert_eq!(grid.quantics_to_grididx(&quantics)?, grididx);
    assert_eq!(grid.grididx_to_quantics(&grididx)?, quantics);
    assert_eq!(grid.origcoord_to_grididx(&[0.0, 0.0, 0.0])?, grididx);

    let quantics = vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1];
    let grididx = vec![2, 1, 1];

    assert_eq!(grid.quantics_to_grididx(&quantics)?, grididx);
    assert_eq!(grid.grididx_to_quantics(&grididx)?, quantics);

    let coord = grid.quantics_to_origcoord(&quantics)?;
    assert!((coord[0] - 1.0 / 16.0).abs() < 1e-12);
    assert!((coord[1] - 0.0).abs() < 1e-12);
    assert!((coord[2] - 0.0).abs() < 1e-12);

    Ok(())
}
```
