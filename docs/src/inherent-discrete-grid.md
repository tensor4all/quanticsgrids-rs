# Inherent Discrete Grid

`InherentDiscreteGrid` is the low-level grid type for discrete target spaces. It
has an interface parallel to `DiscretizedGrid`, but original coordinates are
integer-valued and are defined by an origin and a step size.

```rust
{{#include ../examples/inherent-discrete-grid.rs}}
```

As with `DiscretizedGrid`, the Rust API returns an error when original
coordinates fall outside the representable range.
