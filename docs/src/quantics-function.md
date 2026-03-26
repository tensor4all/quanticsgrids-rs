# Quantics Function

When working with algorithms that expect a function of quantics indices, you can
wrap a coordinate-space function with `quantics_function`.

```rust
{{#include ../examples/quantics-function.rs}}
```

This is the Rust equivalent of first defining a function in the original
coordinate space and then obtaining a quantics-index wrapper around it.
