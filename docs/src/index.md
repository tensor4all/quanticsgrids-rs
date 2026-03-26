# quanticsgrids user guide

This library provides utilities for working with functions on quantics grids and
for converting between quantics indices, grid indices, and original
coordinates.

`quanticsgrids` provides two main layers:

1. Low-level discrete grids with integer coordinates
2. High-level discretized grids over continuous domains

Most users will want the high-level interface built around
`DiscretizedGrid`. The following chapters mirror the Julia guide while using the
Rust API and Rust idioms.
