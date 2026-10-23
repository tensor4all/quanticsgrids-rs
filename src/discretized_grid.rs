//! Discretized grid implementation with continuous coordinate support

use crate::error::{QuanticsGridError, Result};
use crate::inherent_discrete_grid::InherentDiscreteGridBuilder;
use crate::{IndexTable, InherentDiscreteGrid, UnfoldingScheme};

/// A discretized grid with continuous domain and floating-point coordinates.
///
/// This structure wraps an [`InherentDiscreteGrid`] and adds support for
/// continuous coordinate systems with specified lower and upper bounds.
///
/// # Resolution supported by `f64`
///
/// Grid points are `f64` values `lower_bound + g * grid_step`, so the grid is
/// only usable while neighbouring points stay distinguishable and the extreme
/// points still convert to their own index. That limit depends on the bounds, not
/// on the resolution alone: `[0, 1)` supports `R = 50`, symmetric `O(1)` bounds
/// `R = 51`, `[1e12, 1e12 + 1e6)` `R = 30` and `[1e10, 1e10 + 1)` only `R = 17`.
/// [`DiscretizedGrid::builder`] applies that check, which rejects grids whose
/// coordinate round trip is inexact instead of returning a neighbouring index
/// later; within the accepted range
/// `origcoord_to_grididx(&grididx_to_origcoord(&[g])?) == [g]` holds exactly.
///
/// The check is a conservative supported-resolution condition, not a proof of
/// distinguishability: it requires the step to be at least four ulps of the
/// bounds and checks the last two grid points, so a grid whose arithmetic happens
/// to round nowhere (for example `[0, 1)` with a power-of-two step) is still
/// rejected once its step falls below that limit.
///
/// # Example
/// ```
/// use quanticsgrids::{DiscretizedGrid, UnfoldingScheme};
///
/// // Create a 2D grid discretizing [0, 1) x [0, 1)
/// let grid = DiscretizedGrid::builder(&[3, 2])
///     .with_variable_names(&["x", "y"])
///     .build()
///     .unwrap();
///
/// // Convert between coordinate systems
/// let quantics = grid.origcoord_to_quantics(&[0.5, 0.25]).unwrap();
/// let coords = grid.quantics_to_origcoord(&quantics).unwrap();
/// ```
#[derive(Debug, Clone)]
pub struct DiscretizedGrid {
    /// The underlying discrete grid
    discrete_grid: InherentDiscreteGrid,
    /// Lower bounds for each dimension
    lower_bound: Vec<f64>,
    /// Public upper bounds for each dimension
    upper_bound: Vec<f64>,
    /// Grid step for each dimension
    step: Vec<f64>,
    /// Whether each dimension includes the upper endpoint
    include_endpoint: Vec<bool>,
}

impl DiscretizedGrid {
    /// Create a new builder for a DiscretizedGrid.
    ///
    /// # Example
    /// ```
    /// use quanticsgrids::DiscretizedGrid;
    ///
    /// let grid = DiscretizedGrid::builder(&[3, 2])
    ///     .with_variable_names(&["x", "y"])
    ///     .build()
    ///     .unwrap();
    /// ```
    pub fn builder(rs: &[usize]) -> DiscretizedGridBuilder {
        DiscretizedGridBuilder::new(rs)
    }

    /// Create a grid from an explicit index table.
    pub fn from_index_table(
        variable_names: &[&str],
        index_table: IndexTable,
    ) -> DiscretizedGridBuilder {
        DiscretizedGridBuilder::from_index_table(variable_names, index_table)
    }

    /// Number of dimensions
    pub fn ndims(&self) -> usize {
        self.discrete_grid.ndims()
    }

    /// Number of tensor sites (cores)
    pub fn len(&self) -> usize {
        self.discrete_grid.len()
    }

    /// Returns true if the grid has no tensor sites
    pub fn is_empty(&self) -> bool {
        self.discrete_grid.is_empty()
    }

    /// Resolution (number of bits) per dimension
    pub fn rs(&self) -> &[usize] {
        self.discrete_grid.rs()
    }

    /// Variable names
    pub fn variable_names(&self) -> &[String] {
        self.discrete_grid.variable_names()
    }

    /// Per-dimension bases
    pub fn bases(&self) -> &[usize] {
        self.discrete_grid.bases()
    }

    /// Uniform numeric base (panics if bases differ between dimensions)
    pub fn base(&self) -> usize {
        self.discrete_grid.base()
    }

    /// Index table
    pub fn index_table(&self) -> &IndexTable {
        self.discrete_grid.index_table()
    }

    /// Lower bounds for each dimension
    pub fn lower_bound(&self) -> &[f64] {
        &self.lower_bound
    }

    /// Upper bounds for each dimension
    pub fn upper_bound(&self) -> &[f64] {
        &self.upper_bound
    }

    /// Local dimension of a tensor site
    pub fn site_dim(&self, site: usize) -> Result<usize> {
        self.discrete_grid.site_dim(site)
    }

    /// Local dimensions of all tensor sites
    pub fn local_dimensions(&self) -> Vec<usize> {
        self.discrete_grid.local_dimensions()
    }

    /// Grid step size in each dimension
    pub fn grid_step(&self) -> Vec<f64> {
        self.step.clone()
    }

    /// Minimum grid coordinates (same as lower_bound)
    pub fn grid_min(&self) -> &[f64] {
        &self.lower_bound
    }

    /// Maximum representable grid coordinates in each dimension
    pub fn grid_max(&self) -> Vec<f64> {
        (0..self.ndims())
            .map(|dim| self.origcoord_at(dim, self.grid_point_count(dim) - 1))
            .collect()
    }

    /// Get original coordinates for a dimension by index
    pub fn grid_origcoords(&self, dim: usize) -> Result<Vec<f64>> {
        if dim >= self.ndims() {
            return Err(QuanticsGridError::SiteIndexOutOfBounds {
                site: dim,
                max: self.ndims(),
            });
        }

        Ok((0..self.grid_point_count(dim))
            .map(|g| self.origcoord_at(dim, g))
            .collect())
    }

    /// Get original coordinates for a dimension by variable name
    pub fn grid_origcoords_by_name(&self, name: &str) -> Result<Vec<f64>> {
        let dim = self
            .discrete_grid
            .variable_names()
            .iter()
            .position(|n| n == name)
            .ok_or_else(|| QuanticsGridError::UnknownVariable {
                variable: name.to_string(),
                valid: self.discrete_grid.variable_names().to_vec(),
            })?;
        self.grid_origcoords(dim)
    }

    // ========================================================================
    // Core conversion functions
    // ========================================================================

    /// Convert quantics indices to grid indices (0-indexed).
    pub fn quantics_to_grididx(&self, quantics: &[usize]) -> Result<Vec<usize>> {
        self.discrete_grid.quantics_to_grididx(quantics)
    }

    /// Convert grid indices to quantics indices (0-indexed).
    pub fn grididx_to_quantics(&self, grididx: &[usize]) -> Result<Vec<usize>> {
        self.discrete_grid.grididx_to_quantics(grididx)
    }

    /// Convert grid indices (0-indexed) to original coordinates.
    pub fn grididx_to_origcoord(&self, grididx: &[usize]) -> Result<Vec<f64>> {
        let grididx = self.expand_grididx(grididx)?;
        self.validate_grididx(&grididx)?;

        Ok((0..self.ndims())
            .map(|dim| self.origcoord_at(dim, grididx[dim]))
            .collect())
    }

    /// Convert original coordinates to grid indices (0-indexed).
    ///
    /// Each coordinate maps to the index `round((c - lower_bound) / grid_step)`,
    /// clamped to `0..=N-1`. A normalized index exactly halfway between two
    /// indices rounds to the larger one. The rounding applies to that computed
    /// ratio, which is not bit-for-bit the midpoint between two emitted
    /// coordinates, so the result is not promised to be the geometrically
    /// nearest grid point in every case.
    ///
    /// The grid is not periodic: in the half-open case a coordinate in the last
    /// step below `upper_bound` maps to `N-1`, the last grid point, rather than
    /// wrapping to `0`.
    ///
    /// # Errors
    ///
    /// Returns [`QuanticsGridError::CoordinateOutOfBounds`] when a coordinate
    /// lies outside the configured domain (`[lower, upper)`, or
    /// `[lower, upper]` with `include_endpoint(true)`), and
    /// [`QuanticsGridError::DimensionMismatch`] when `coord` has a length other
    /// than one or the number of dimensions.
    pub fn origcoord_to_grididx(&self, coord: &[f64]) -> Result<Vec<usize>> {
        let coord = self.expand_coord(coord)?;
        self.validate_origcoord(&coord)?;

        Ok((0..self.ndims())
            .map(|dim| self.grididx_at(dim, coord[dim]))
            .collect())
    }

    /// Convert original coordinates to quantics indices (0-indexed).
    pub fn origcoord_to_quantics(&self, coord: &[f64]) -> Result<Vec<usize>> {
        let grididx = self.origcoord_to_grididx(coord)?;
        self.grididx_to_quantics(&grididx)
    }

    /// Convert quantics indices to original coordinates.
    pub fn quantics_to_origcoord(&self, quantics: &[usize]) -> Result<Vec<f64>> {
        let grididx = self.quantics_to_grididx(quantics)?;
        self.grididx_to_origcoord(&grididx)
    }

    // ========================================================================
    // Private helper methods
    // ========================================================================

    /// Number of grid points in one dimension.
    fn grid_point_count(&self, dim: usize) -> usize {
        self.discrete_grid.bases()[dim].pow(self.discrete_grid.rs()[dim] as u32)
    }

    /// Coordinate of grid index `g` (0-indexed) in one dimension.
    fn origcoord_at(&self, dim: usize, g: usize) -> f64 {
        // The included endpoint is the upper bound by definition, while rounding
        // of `lower_bound + (N-1) * step` can miss it by an ulp in either
        // direction.
        if self.include_endpoint[dim] && g == self.grid_point_count(dim) - 1 {
            return self.upper_bound[dim];
        }
        self.lower_bound[dim] + g as f64 * self.step[dim]
    }

    /// Grid index nearest to the coordinate `c` (0-indexed), clamped to the grid.
    ///
    /// This is the inverse of [`Self::origcoord_at`] as far as `f64` allows; see
    /// [`Self::origcoord_to_grididx`] for the rounding rule.
    fn grididx_at(&self, dim: usize, c: f64) -> usize {
        let idx = ((c - self.lower_bound[dim]) / self.step[dim]).round() as i64;
        idx.clamp(0, self.grid_point_count(dim) as i64 - 1) as usize
    }

    fn validate_grididx(&self, grididx: &[usize]) -> Result<()> {
        let max_grididx = self.discrete_grid.max_grididx();
        for (dim, (&val, &max)) in grididx.iter().zip(max_grididx.iter()).enumerate() {
            if val > max {
                return Err(QuanticsGridError::GridIndexOutOfBounds {
                    dim,
                    value: val,
                    max,
                });
            }
        }
        Ok(())
    }

    fn validate_origcoord(&self, coord: &[f64]) -> Result<()> {
        for (dim, (((&c, &lo), &hi), &include_endpoint)) in coord
            .iter()
            .zip(self.lower_bound.iter())
            .zip(self.upper_bound.iter())
            .zip(self.include_endpoint.iter())
            .enumerate()
        {
            let out_of_bounds = if include_endpoint {
                c < lo || c > hi
            } else {
                c < lo || c >= hi
            };

            if out_of_bounds {
                return Err(QuanticsGridError::CoordinateOutOfBounds {
                    dim,
                    value: c,
                    lower: lo,
                    upper: hi,
                });
            }
        }
        Ok(())
    }

    /// Expand a scalar or lower-dimensional input to full dimensions
    fn expand_grididx(&self, grididx: &[usize]) -> Result<Vec<usize>> {
        let ndims = self.ndims();
        if grididx.len() == 1 && ndims > 1 {
            Ok(vec![grididx[0]; ndims])
        } else if grididx.len() == ndims {
            Ok(grididx.to_vec())
        } else {
            Err(QuanticsGridError::DimensionMismatch {
                expected: ndims,
                actual: grididx.len(),
            })
        }
    }

    fn expand_coord(&self, coord: &[f64]) -> Result<Vec<f64>> {
        let ndims = self.ndims();
        if coord.len() == 1 && ndims > 1 {
            Ok(vec![coord[0]; ndims])
        } else if coord.len() == ndims {
            Ok(coord.to_vec())
        } else {
            Err(QuanticsGridError::DimensionMismatch {
                expected: ndims,
                actual: coord.len(),
            })
        }
    }
}

impl std::fmt::Display for DiscretizedGrid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let ndims = self.ndims();
        let rs = self.discrete_grid.rs();
        let bases = self.discrete_grid.bases();

        // Calculate total points
        let total_points: u64 = rs
            .iter()
            .zip(bases.iter())
            .map(|(&r, &b)| (b as u64).pow(r as u32))
            .product();

        if ndims <= 1 {
            write!(
                f,
                "DiscretizedGrid{{{}}} with {} grid points",
                ndims, total_points
            )?;
        } else {
            let sizes: Vec<String> = rs
                .iter()
                .zip(bases.iter())
                .map(|(&r, &b)| format!("{}", (b as u64).pow(r as u32)))
                .collect();
            write!(
                f,
                "DiscretizedGrid{{{}}} with {} = {} grid points",
                ndims,
                sizes.join(" x "),
                total_points
            )?;
        }

        // Variable names
        let var_names = self.discrete_grid.variable_names();
        if var_names
            .iter()
            .any(|n| !n.chars().all(|c| c.is_ascii_digit()))
        {
            write!(f, "\n  Variables: ({})", var_names.join(", "))?;
        }

        // Resolutions
        if ndims == 1 {
            write!(f, "\n  Resolution: {} bits", rs[0])?;
        } else {
            let res_str: Vec<String> = var_names
                .iter()
                .zip(rs.iter())
                .map(|(n, &r)| format!("{}: {}", n, r))
                .collect();
            write!(f, "\n  Resolutions: ({})", res_str.join(", "))?;
        }

        // Domain
        let step = self.grid_step();
        if ndims == 1 {
            let right_bracket = if self.include_endpoint[0] { ']' } else { ')' };
            write!(
                f,
                "\n  Domain: [{}, {}{}",
                self.lower_bound[0], self.upper_bound[0], right_bracket
            )?;
            write!(f, "\n  Grid spacing: {}", step[0])?;
        } else {
            let bounds_str: Vec<String> = self
                .lower_bound
                .iter()
                .zip(self.upper_bound.iter())
                .zip(self.include_endpoint.iter())
                .map(|((&lo, &hi), &include_endpoint)| {
                    let right_bracket = if include_endpoint { ']' } else { ')' };
                    format!("[{}, {}{}", lo, hi, right_bracket)
                })
                .collect();
            write!(f, "\n  Domain: {}", bounds_str.join(" x "))?;

            let step_str: Vec<String> = var_names
                .iter()
                .zip(step.iter())
                .map(|(n, &s)| format!("d{} = {}", n, s))
                .collect();
            write!(f, "\n  Grid spacing: ({})", step_str.join(", "))?;
        }

        // Base (if not all binary)
        let all_same_base = bases.iter().all(|&b| b == bases[0]);
        if all_same_base && bases[0] != 2 {
            write!(f, "\n  Base: {}", bases[0])?;
        } else if !all_same_base {
            let base_str: Vec<String> = bases.iter().map(|b| b.to_string()).collect();
            write!(f, "\n  Bases: ({})", base_str.join(", "))?;
        }

        // Tensor structure
        let num_sites = self.len();
        let sitedims = self.local_dimensions();
        write!(f, "\n  Tensor train: {} sites", num_sites)?;
        if !sitedims.is_empty() {
            if sitedims.iter().all(|&d| d == sitedims[0]) {
                write!(f, " (uniform dimension {})", sitedims[0])?;
            } else {
                let dims_str: Vec<String> = sitedims.iter().map(|d| d.to_string()).collect();
                write!(f, " (dimensions: {})", dims_str.join("-"))?;
            }
        }

        Ok(())
    }
}

/// Builder for DiscretizedGrid
#[derive(Debug, Clone)]
pub struct DiscretizedGridBuilder {
    inner: InherentDiscreteGridBuilder,
    lower_bound: Option<Vec<f64>>,
    upper_bound: Option<Vec<f64>>,
    include_endpoint: Option<Vec<bool>>,
}

impl DiscretizedGridBuilder {
    /// Create a new builder with given resolutions
    pub fn new(rs: &[usize]) -> Self {
        Self {
            inner: InherentDiscreteGridBuilder::new(rs),
            lower_bound: None,
            upper_bound: None,
            include_endpoint: None,
        }
    }

    /// Create a builder from an explicit index table
    pub fn from_index_table(variable_names: &[&str], index_table: IndexTable) -> Self {
        Self {
            inner: InherentDiscreteGridBuilder::from_index_table(variable_names, index_table),
            lower_bound: None,
            upper_bound: None,
            include_endpoint: None,
        }
    }

    /// Set the lower bounds for each dimension
    pub fn with_lower_bound(mut self, lower_bound: &[f64]) -> Self {
        self.lower_bound = Some(lower_bound.to_vec());
        self
    }

    /// Set the upper bounds for each dimension
    pub fn with_upper_bound(mut self, upper_bound: &[f64]) -> Self {
        self.upper_bound = Some(upper_bound.to_vec());
        self
    }

    /// Set bounds for 1D case
    pub fn with_bounds(mut self, lower: f64, upper: f64) -> Self {
        self.lower_bound = Some(vec![lower]);
        self.upper_bound = Some(vec![upper]);
        self
    }

    /// Set whether to include the endpoint for each dimension
    pub fn with_include_endpoint(mut self, include: &[bool]) -> Self {
        self.include_endpoint = Some(include.to_vec());
        self
    }

    /// Set whether to include the endpoint (single value for all dimensions)
    pub fn include_endpoint(mut self, include: bool) -> Self {
        // Will be expanded to match dimensions during build
        self.include_endpoint = Some(vec![include]);
        self
    }

    /// Set variable names
    pub fn with_variable_names(mut self, names: &[&str]) -> Self {
        self.inner = self.inner.with_variable_names(names);
        self
    }

    /// Set a uniform numeric base (default 2)
    pub fn with_base(mut self, base: usize) -> Self {
        self.inner = self.inner.with_base(base);
        self
    }

    /// Set per-dimension bases
    pub fn with_bases(mut self, bases: &[usize]) -> Self {
        self.inner = self.inner.with_bases(bases);
        self
    }

    /// Set the unfolding scheme
    pub fn with_unfolding_scheme(mut self, scheme: UnfoldingScheme) -> Self {
        self.inner = self.inner.with_unfolding_scheme(scheme);
        self
    }

    /// Build the DiscretizedGrid
    ///
    /// # Errors
    ///
    /// - [`QuanticsGridError::DimensionMismatch`] when `lower_bound`,
    ///   `upper_bound` or `include_endpoint` have a length other than one or the
    ///   number of dimensions.
    /// - [`QuanticsGridError::InvalidBounds`] when a lower bound is not below its
    ///   upper bound.
    /// - [`QuanticsGridError::EndpointWithZeroResolution`] when a dimension with
    ///   resolution 0 includes its endpoint.
    /// - [`QuanticsGridError::GridNotRepresentable`] when the grid cannot be
    ///   represented faithfully in `f64` over the bounds. This check is
    ///   conservative: it requires the step to be at least four ulps of the
    ///   bounds and checks that the last two grid points stay inside the domain
    ///   and map back to their own index. See [`DiscretizedGrid`] for the
    ///   supported resolutions per bound regime.
    /// - [`QuanticsGridError::ResolutionTooLarge`] and
    ///   [`QuanticsGridError::InvalidBase`] from the underlying
    ///   [`InherentDiscreteGrid`].
    pub fn build(self) -> Result<DiscretizedGrid> {
        let discrete_grid = self.inner.build()?;
        let ndims = discrete_grid.ndims();
        let rs = discrete_grid.rs();
        let bases = discrete_grid.bases();

        // Default bounds: [0, 1) for each dimension
        let lower_bound = self.lower_bound.unwrap_or_else(|| vec![0.0; ndims]);
        let upper_bound = self.upper_bound.unwrap_or_else(|| vec![1.0; ndims]);

        // Expand bounds if needed
        let lower_bound = if lower_bound.len() == 1 && ndims > 1 {
            vec![lower_bound[0]; ndims]
        } else {
            lower_bound
        };

        let upper_bound = if upper_bound.len() == 1 && ndims > 1 {
            vec![upper_bound[0]; ndims]
        } else {
            upper_bound
        };

        // Validate dimensions
        if lower_bound.len() != ndims {
            return Err(QuanticsGridError::DimensionMismatch {
                expected: ndims,
                actual: lower_bound.len(),
            });
        }
        if upper_bound.len() != ndims {
            return Err(QuanticsGridError::DimensionMismatch {
                expected: ndims,
                actual: upper_bound.len(),
            });
        }

        // Validate bounds
        for d in 0..ndims {
            if lower_bound[d] >= upper_bound[d] {
                return Err(QuanticsGridError::InvalidBounds {
                    dim: d,
                    lower: lower_bound[d],
                    upper: upper_bound[d],
                });
            }
        }

        // Handle endpoint inclusion
        let include_endpoint = match self.include_endpoint {
            Some(v) if v.len() == 1 && ndims > 1 => vec![v[0]; ndims],
            Some(v) if v.len() == ndims => v,
            Some(v) if v.len() != ndims => {
                return Err(QuanticsGridError::DimensionMismatch {
                    expected: ndims,
                    actual: v.len(),
                });
            }
            _ => vec![false; ndims],
        };

        // The step is the single source of truth for the coordinate conversion,
        // so it is computed and validated once here instead of being re-derived
        // from a rounded upper bound at every call.
        let mut step = Vec::with_capacity(ndims);
        for d in 0..ndims {
            let span = upper_bound[d] - lower_bound[d];
            let n_points = (bases[d] as f64).powi(rs[d] as i32);
            let step_d = if include_endpoint[d] {
                if rs[d] == 0 {
                    return Err(QuanticsGridError::EndpointWithZeroResolution { dim: d });
                }
                span / (n_points - 1.0)
            } else {
                span / n_points
            };

            // Neighbouring grid points must stay distinguishable in f64. Once the
            // step is only a few ulps wide, several indices collapse onto one
            // coordinate and no conversion can recover them.
            step.push(step_d);
        }

        let grid = DiscretizedGrid {
            discrete_grid,
            lower_bound,
            upper_bound,
            step,
            include_endpoint,
        };

        for d in 0..grid.ndims() {
            let n = grid.grid_point_count(d);
            let too_fine =
                grid.step[d] < 4.0 * ulp(grid.lower_bound[d].abs().max(grid.upper_bound[d].abs()));

            // Rounding can also put an extreme grid point outside the domain or
            // make it convert back to a neighbouring index. Checking both ends of
            // the grid covers the indices where the rounding error is largest.
            let broken_ends = [n - 1, n.saturating_sub(2)].iter().any(|&g| {
                let coord = grid.origcoord_at(d, g);
                let inside = coord >= grid.lower_bound[d]
                    && if grid.include_endpoint[d] {
                        coord <= grid.upper_bound[d]
                    } else {
                        coord < grid.upper_bound[d]
                    };
                !inside || grid.grididx_at(d, coord) != g
            });

            if too_fine || broken_ends {
                return Err(QuanticsGridError::GridNotRepresentable {
                    dim: d,
                    step: grid.step[d],
                    lower: grid.lower_bound[d],
                    upper: grid.upper_bound[d],
                });
            }
        }

        Ok(grid)
    }
}

/// Distance from a finite `f64` to the next representable value.
fn ulp(x: f64) -> f64 {
    let x = x.abs();
    if x == 0.0 {
        f64::MIN_POSITIVE
    } else {
        f64::from_bits(x.to_bits() + 1) - x
    }
}

/// Wrap a function to accept quantics indices
pub fn quantics_function<F>(grid: &DiscretizedGrid, f: F) -> impl Fn(&[usize]) -> Result<f64> + '_
where
    F: Fn(&[f64]) -> f64 + 'static,
{
    move |quantics: &[usize]| {
        let coords = grid.quantics_to_origcoord(quantics)?;
        Ok(f(&coords))
    }
}

#[cfg(test)]
mod tests;
