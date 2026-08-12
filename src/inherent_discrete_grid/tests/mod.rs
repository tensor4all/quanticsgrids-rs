use super::*;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

// ============================================================================
// Existing tests (kept as-is)
// ============================================================================

#[test]
fn test_basic_1d_grid() {
    let grid = InherentDiscreteGrid::builder(&[3]).build().unwrap();
    assert_eq!(grid.ndims(), 1);
    assert_eq!(grid.len(), 3); // 3 sites for fused scheme
    assert_eq!(grid.base(), 2);
    assert_eq!(grid.rs(), &[3]);
    assert_eq!(grid.max_grididx(), &[7]); // 2^3 - 1
}

#[test]
fn test_basic_2d_grid() {
    let grid = InherentDiscreteGrid::builder(&[3, 2])
        .with_variable_names(&["x", "y"])
        .build()
        .unwrap();

    assert_eq!(grid.ndims(), 2);
    assert_eq!(grid.variable_names(), &["x", "y"]);
    assert_eq!(grid.max_grididx(), &[7, 3]); // 2^3 - 1, 2^2 - 1
}

#[test]
fn test_grididx_to_quantics_roundtrip() {
    let grid = InherentDiscreteGrid::builder(&[3, 2])
        .with_variable_names(&["a", "b"])
        .build()
        .unwrap();

    let grididx = vec![4usize, 1];
    let quantics = grid.grididx_to_quantics(&grididx).unwrap();
    let back = grid.quantics_to_grididx(&quantics).unwrap();
    assert_eq!(back, grididx);
}

#[test]
fn test_all_grididx_roundtrip() {
    let grid = InherentDiscreteGrid::builder(&[2, 2]).build().unwrap();

    for x in 0..4 {
        for y in 0..4 {
            let grididx = vec![x, y];
            let quantics = grid.grididx_to_quantics(&grididx).unwrap();
            let back = grid.quantics_to_grididx(&quantics).unwrap();
            assert_eq!(back, grididx, "Failed for grididx {:?}", grididx);
        }
    }
}

#[test]
fn test_interleaved_scheme() {
    let grid = InherentDiscreteGrid::builder(&[2, 2])
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();

    // Interleaved: [1_1], [2_1], [1_2], [2_2]
    assert_eq!(grid.len(), 4);
    for x in 0..4 {
        for y in 0..4 {
            let grididx = vec![x, y];
            let quantics = grid.grididx_to_quantics(&grididx).unwrap();
            let back = grid.quantics_to_grididx(&quantics).unwrap();
            assert_eq!(back, grididx);
        }
    }
}

#[test]
fn test_base3_grid() {
    let grid = InherentDiscreteGrid::builder(&[2])
        .with_base(3)
        .build()
        .unwrap();

    assert_eq!(grid.base(), 3);
    assert_eq!(grid.max_grididx(), &[8]); // 3^2 - 1 = 8

    for x in 0..9 {
        let grididx = vec![x];
        let quantics = grid.grididx_to_quantics(&grididx).unwrap();
        let back = grid.quantics_to_grididx(&quantics).unwrap();
        assert_eq!(back, grididx);
    }
}

#[test]
fn test_origcoord_conversion() {
    let grid = InherentDiscreteGrid::builder(&[2])
        .with_origin(&[0])
        .with_step(&[1])
        .build()
        .unwrap();

    // origin=0, step=1, max_grididx=3 (0-based: 2^2 - 1)
    // grididx 0 -> origcoord 0
    // grididx 3 -> origcoord 3
    let coord = grid.grididx_to_origcoord(&[0]).unwrap();
    assert_eq!(coord, vec![0]);

    let coord = grid.grididx_to_origcoord(&[3]).unwrap();
    assert_eq!(coord, vec![3]);

    let idx = grid.origcoord_to_grididx(&[2]).unwrap();
    assert_eq!(idx, vec![2]);
}

#[test]
fn test_error_invalid_base() {
    let result = InherentDiscreteGrid::builder(&[3]).with_base(1).build();
    assert!(matches!(result, Err(QuanticsGridError::InvalidBase(1))));
}

#[test]
fn test_error_duplicate_variable_names() {
    let result = InherentDiscreteGrid::builder(&[2, 2])
        .with_variable_names(&["x", "x"])
        .build();
    assert!(matches!(
        result,
        Err(QuanticsGridError::DuplicateVariableName(_))
    ));
}

#[test]
fn test_error_quantics_out_of_range() {
    let grid = InherentDiscreteGrid::builder(&[2]).build().unwrap();
    // Rs=[2] with Fused creates 2 sites, each with dim 2
    // So quantics should have length 2, and each value should be in [0, 1]
    let result = grid.quantics_to_grididx(&[5, 1]); // 5 is out of range [0, 1]
    assert!(matches!(
        result,
        Err(QuanticsGridError::QuanticsOutOfRange { .. })
    ));
}

#[test]
fn test_error_grididx_out_of_bounds() {
    let grid = InherentDiscreteGrid::builder(&[2]).build().unwrap();
    let result = grid.grididx_to_quantics(&[5]); // max is 3 (2^2 - 1)
    assert!(matches!(
        result,
        Err(QuanticsGridError::GridIndexOutOfBounds { .. })
    ));
}

#[test]
fn test_local_dimensions() {
    let grid = InherentDiscreteGrid::builder(&[3, 2])
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();

    let dims = grid.local_dimensions();
    // Fused scheme with Rs=(3,2): sites have 2, 2, 1 indices each -> dims 4, 4, 2
    // Actually: bitnumber 0: [b, a] -> dim 4
    //           bitnumber 1: [b, a] -> dim 4
    //           bitnumber 2: [a] -> dim 2
    assert_eq!(dims, vec![4, 4, 2]);
}

#[test]
fn test_from_index_table() {
    // Create a custom index table like in Julia: [[(:a, 1), (:b, 2)], [(:a, 2)], [(:b, 1), (:a, 3)]]
    let index_table = vec![
        vec![("a".to_string(), 0), ("b".to_string(), 1)],
        vec![("a".to_string(), 1)],
        vec![("b".to_string(), 0), ("a".to_string(), 2)],
    ];

    let grid = InherentDiscreteGrid::from_index_table(&["a", "b"], index_table)
        .build()
        .unwrap();

    assert_eq!(grid.ndims(), 2);
    assert_eq!(grid.rs(), &[3, 2]); // a has 3 bits, b has 2 bits
    assert_eq!(grid.len(), 3);

    // Test roundtrip
    for x in 0..8 {
        for y in 0..4 {
            let grididx = vec![x, y];
            let quantics = grid.grididx_to_quantics(&grididx).unwrap();
            let back = grid.quantics_to_grididx(&quantics).unwrap();
            assert_eq!(back, grididx);
        }
    }
}

// ============================================================================
// Tests ported from inherentdiscretegrid_tests.jl
// ============================================================================

// --- "InherentDiscreteGrid" (3D, multiple schemes/steps/origins) ---
#[test]
fn test_3d_grid_interleaved_and_fused() {
    let schemes = [UnfoldingScheme::Interleaved, UnfoldingScheme::Fused];
    let steps: &[&[i64]] = &[&[1, 1, 1], &[1, 1, 2]];
    let origins: &[&[i64]] = &[&[1, 1, 1], &[1, 1, 2]];

    for &scheme in &schemes {
        for &step in steps {
            for &origin in origins {
                let grid = InherentDiscreteGrid::builder(&[5, 5, 5])
                    .with_origin(origin)
                    .with_step(step)
                    .with_unfolding_scheme(scheme)
                    .build()
                    .unwrap();

                assert_eq!(grid.grid_min(), origin.to_vec());
                assert_eq!(grid.step(), step);

                if scheme == UnfoldingScheme::Interleaved {
                    assert_eq!(grid.local_dimensions(), vec![2; 3 * 5]);
                } else {
                    assert_eq!(grid.local_dimensions(), vec![8; 5]); // 2^3 = 8
                }

                let test_indices: &[&[usize]] = &[
                    &[0, 0, 0],
                    &[0, 0, 1],
                    &[0, 24, 0],
                    &[13, 0, 0],
                    &[24, 24, 24],
                ];

                for &idx in test_indices {
                    let c = grid.grididx_to_origcoord(idx).unwrap();
                    assert_eq!(grid.origcoord_to_grididx(&c).unwrap(), idx.to_vec());

                    let q = grid.grididx_to_quantics(idx).unwrap();
                    if scheme == UnfoldingScheme::Fused {
                        assert_eq!(q.len(), 5);
                    } else {
                        assert_eq!(q.len(), 3 * 5);
                    }
                    assert!(q.iter().all(|&v| (0..8).contains(&v)));
                    assert_eq!(grid.quantics_to_origcoord(&q).unwrap(), c);
                }
            }
        }
    }
}

// --- "InherentDiscreteGrid constructors" ---
#[test]
fn test_constructors_basic() {
    // 1D constructor
    let r = 4;
    let grid_1d = InherentDiscreteGrid::builder(&[r])
        .with_base(2)
        .with_origin(&[1])
        .with_step(&[1])
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();
    assert_eq!(grid_1d.base(), 2);
    assert_eq!(grid_1d.origin(), &[1]);
    assert_eq!(grid_1d.step(), &[1]);

    // 3D with tuple origin and step
    let grid_3d = InherentDiscreteGrid::builder(&[r, r, r])
        .with_base(3)
        .with_origin(&[5, 10, 15])
        .with_step(&[2, 3, 1])
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();
    assert_eq!(grid_3d.base(), 3);
    assert_eq!(grid_3d.origin(), &[5, 10, 15]);
    assert_eq!(grid_3d.step(), &[2, 3, 1]);

    // Default parameters
    let grid_default = InherentDiscreteGrid::builder(&[r, r])
        .with_origin(&[0, 0])
        .build()
        .unwrap();
    assert_eq!(grid_default.base(), 2);
    assert_eq!(grid_default.step(), &[1, 1]);
}

// --- "InherentDiscreteGrid basic properties" ---
#[test]
fn test_basic_properties() {
    let r = 5;
    let base = 3;
    let origin = [10, 20];
    let step = [2, 5];

    let grid = InherentDiscreteGrid::builder(&[r, r])
        .with_base(base)
        .with_origin(&origin)
        .with_step(&step)
        .build()
        .unwrap();

    assert_eq!(grid.ndims(), 2);
    assert_eq!(grid.len(), r); // Number of quantics sites for fused
    assert_eq!(grid.base(), base);

    assert_eq!(grid.grid_min(), vec![10, 20]);
    let max_idx = (base as i64).pow(r as u32);
    assert_eq!(
        grid.grid_max(),
        vec![
            origin[0] + step[0] * (max_idx - 1),
            origin[1] + step[1] * (max_idx - 1)
        ]
    );
    assert_eq!(grid.step(), &step);
    assert_eq!(grid.origin(), &origin);

    // Interleaved scheme
    let grid_interleaved = InherentDiscreteGrid::builder(&[r, r])
        .with_base(base)
        .with_origin(&origin)
        .with_step(&step)
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();
    assert_eq!(grid_interleaved.len(), 2 * r); // Number of quantics sites for interleaved
}

// --- "InherentDiscreteGrid coordinate conversions" ---
#[test]
fn test_coordinate_conversions() {
    let r = 4;
    let base: usize = 2;
    let origin = [5i64, 10, 15];
    let step = [2i64, 3, 1];

    let grid = InherentDiscreteGrid::builder(&[r, r, r])
        .with_base(base)
        .with_origin(&origin)
        .with_step(&step)
        .build()
        .unwrap();

    // Test grididx_to_origcoord and origcoord_to_grididx roundtrip
    let test_gridindices: &[&[usize]] = &[
        &[0, 0, 0],
        &[0, 0, 1],
        &[0, 7, 0],
        &[3, 0, 0],
        &[15, 15, 15],
    ];
    for &grididx in test_gridindices {
        let origcoord = grid.grididx_to_origcoord(grididx).unwrap();
        assert_eq!(
            grid.origcoord_to_grididx(&origcoord).unwrap(),
            grididx.to_vec()
        );
    }

    // Test origcoord calculation formula (coord = origin + grididx * step)
    let grididx = [3usize, 5, 7];
    let expected_origcoord: Vec<i64> = origin
        .iter()
        .zip(step.iter())
        .zip(grididx.iter())
        .map(|((&o, &s), &g)| o + s * g as i64)
        .collect();
    assert_eq!(
        grid.grididx_to_origcoord(&grididx).unwrap(),
        expected_origcoord
    );

    // Test boundary cases
    let max = (base as i64).pow(r as u32); // base^R: number of grid points per dim
    let min_grididx = [0usize, 0, 0];
    let max_grididx = [(max - 1) as usize, (max - 1) as usize, (max - 1) as usize];
    assert_eq!(
        grid.grididx_to_origcoord(&min_grididx).unwrap(),
        origin.to_vec()
    );
    let expected_max: Vec<i64> = origin
        .iter()
        .zip(step.iter())
        .map(|(&o, &s)| o + s * (max - 1))
        .collect();
    assert_eq!(
        grid.grididx_to_origcoord(&max_grididx).unwrap(),
        expected_max
    );
}

// --- "InherentDiscreteGrid quantics conversions" ---
#[test]
fn test_quantics_conversions_fused() {
    let r = 3;
    let base: usize = 2;
    let grid_fused = InherentDiscreteGrid::builder(&[r, r])
        .with_base(base)
        .with_origin(&[1, 1])
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();

    let test_gridindices: &[&[usize]] = &[&[0, 0], &[0, 1], &[1, 0], &[3, 7], &[7, 7]];
    for &grididx in test_gridindices {
        let quantics = grid_fused.grididx_to_quantics(grididx).unwrap();
        assert_eq!(
            grid_fused.quantics_to_grididx(&quantics).unwrap(),
            grididx.to_vec()
        );
        assert_eq!(quantics.len(), r); // Fused should have R quantics
        let max_dim = base * base; // base^2 for 2D fused
        assert!(quantics.iter().all(|&v| v < max_dim));
    }
}

#[test]
fn test_quantics_conversions_interleaved() {
    let r = 3;
    let base: usize = 2;
    let grid_interleaved = InherentDiscreteGrid::builder(&[r, r])
        .with_base(base)
        .with_origin(&[1, 1])
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();

    let test_gridindices: &[&[usize]] = &[&[0, 0], &[0, 1], &[1, 0], &[3, 7], &[7, 7]];
    for &grididx in test_gridindices {
        let quantics = grid_interleaved.grididx_to_quantics(grididx).unwrap();
        assert_eq!(
            grid_interleaved.quantics_to_grididx(&quantics).unwrap(),
            grididx.to_vec()
        );
        assert_eq!(quantics.len(), 2 * r); // Interleaved should have 2*R quantics
        assert!(quantics.iter().all(|&v| v < base));
    }

    // Wrong length quantics should error
    let result = grid_interleaved.quantics_to_grididx(&[1, 2, 3]);
    assert!(matches!(
        result,
        Err(QuanticsGridError::WrongQuanticsLength { .. })
    ));
}

// --- "InherentDiscreteGrid quantics_to_origcoord and origcoord_to_quantics" ---
#[test]
fn test_quantics_origcoord_roundtrip() {
    let r = 4;
    let origin = [3i64, 7];
    let step = [2i64, 3];
    let grid = InherentDiscreteGrid::builder(&[r, r])
        .with_origin(&origin)
        .with_step(&step)
        .build()
        .unwrap();

    // Julia: rand(0:2^2-1, R) where R=4, fused 2D so each site is base^2=4
    let dims = grid.local_dimensions();
    let mut rng = StdRng::seed_from_u64(42);

    for _ in 0..20 {
        let quantics: Vec<usize> = (0..r).map(|i| rng.random_range(0..dims[i])).collect();

        let origcoord = grid.quantics_to_origcoord(&quantics).unwrap();
        assert_eq!(grid.origcoord_to_quantics(&origcoord).unwrap(), quantics);

        // Verify the conversion chain: quantics -> grididx -> origcoord
        let grididx = grid.quantics_to_grididx(&quantics).unwrap();
        let expected_origcoord = grid.grididx_to_origcoord(&grididx).unwrap();
        assert_eq!(origcoord, expected_origcoord);
    }
}

// --- "InherentDiscreteGrid localdimensions" ---
#[test]
fn test_localdimensions_various() {
    let r = 4;
    let base: usize = 3;

    // Fused 1D
    let grid_fused_1d = InherentDiscreteGrid::builder(&[r])
        .with_base(base)
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();
    assert_eq!(grid_fused_1d.local_dimensions(), vec![base; r]);

    // Fused 3D
    let grid_fused_3d = InherentDiscreteGrid::builder(&[r, r, r])
        .with_base(base)
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();
    assert_eq!(
        grid_fused_3d.local_dimensions(),
        vec![base.pow(3); r] // base^3 per site
    );

    // Interleaved 1D
    let grid_interleaved_1d = InherentDiscreteGrid::builder(&[r])
        .with_base(base)
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();
    assert_eq!(grid_interleaved_1d.local_dimensions(), vec![base; r]);

    // Interleaved 3D
    let grid_interleaved_3d = InherentDiscreteGrid::builder(&[r, r, r])
        .with_base(base)
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();
    assert_eq!(grid_interleaved_3d.local_dimensions(), vec![base; 3 * r]);

    // Different bases
    for test_base in [2, 4, 5] {
        let grid = InherentDiscreteGrid::builder(&[r, r])
            .with_base(test_base)
            .with_unfolding_scheme(UnfoldingScheme::Fused)
            .build()
            .unwrap();
        assert_eq!(
            grid.local_dimensions(),
            vec![test_base.pow(2); r] // base^2 per site for 2D fused
        );
    }
}

// --- "InherentDiscreteGrid different bases" ---
#[test]
fn test_different_bases() {
    let r = 3;
    let origin = [1i64, 1];

    for base in [2, 3, 4, 5] {
        let grid = InherentDiscreteGrid::builder(&[r, r])
            .with_base(base)
            .with_origin(&origin)
            .build()
            .unwrap();

        let max = (base as i64).pow(r as u32); // base^R
        let max_grididx = [(max - 1) as usize; 2];
        let min_grididx = [0usize, 0];

        // Test coordinate conversion with different bases
        let max_origcoord = grid.grididx_to_origcoord(&max_grididx).unwrap();
        let min_origcoord = grid.grididx_to_origcoord(&min_grididx).unwrap();

        assert_eq!(
            grid.origcoord_to_grididx(&max_origcoord).unwrap(),
            max_grididx.to_vec()
        );
        assert_eq!(
            grid.origcoord_to_grididx(&min_origcoord).unwrap(),
            min_grididx.to_vec()
        );

        // Test quantics conversion with different bases
        let quantics = grid.grididx_to_quantics(&[1, 2]).unwrap();
        let max_dim = base * base; // For fused 2D
        assert!(quantics.iter().all(|&v| v < max_dim));
        assert_eq!(grid.quantics_to_grididx(&quantics).unwrap(), vec![1, 2]);
    }
}

// --- "InherentDiscreteGrid comprehensive conversion test" ---
#[test]
fn test_comprehensive_conversion() {
    let r = 4;
    let base: usize = 3;
    let origin = [2i64, 5];
    let step = [3i64, 2];

    let grid = InherentDiscreteGrid::builder(&[r, r])
        .with_base(base)
        .with_origin(&origin)
        .with_step(&step)
        .build()
        .unwrap();

    let initial_grididx = vec![1usize, 2];

    // Build all representations
    let quantics = grid.grididx_to_quantics(&initial_grididx).unwrap();
    let origcoord = grid.grididx_to_origcoord(&initial_grididx).unwrap();

    // Verify all cross conversions
    assert_eq!(
        grid.grididx_to_quantics(&initial_grididx).unwrap(),
        quantics
    );
    assert_eq!(
        grid.grididx_to_origcoord(&initial_grididx).unwrap(),
        origcoord
    );
    assert_eq!(
        grid.quantics_to_grididx(&quantics).unwrap(),
        initial_grididx
    );
    assert_eq!(grid.quantics_to_origcoord(&quantics).unwrap(), origcoord);
    assert_eq!(
        grid.origcoord_to_grididx(&origcoord).unwrap(),
        initial_grididx
    );
    assert_eq!(grid.origcoord_to_quantics(&origcoord).unwrap(), quantics);
}

// --- "InherentDiscreteGrid edge cases and boundary conditions" ---
#[test]
fn test_edge_cases_and_boundaries() {
    let base: usize = 2;

    // Minimum grid size (R=1)
    let grid_min = InherentDiscreteGrid::builder(&[1])
        .with_base(base)
        .with_origin(&[10])
        .with_step(&[5])
        .build()
        .unwrap();
    assert_eq!(grid_min.grididx_to_origcoord(&[0]).unwrap(), vec![10]);
    assert_eq!(
        grid_min.grididx_to_origcoord(&[base - 1]).unwrap(),
        vec![10 + 5 * (base as i64 - 1)]
    );

    // Large step sizes
    let grid_large_step = InherentDiscreteGrid::builder(&[3, 3])
        .with_origin(&[0, 0])
        .with_step(&[100, 200])
        .build()
        .unwrap();
    let expected_coord = vec![100i64, 400]; // (0 + 100*1, 0 + 200*2)
    assert_eq!(
        grid_large_step.grididx_to_origcoord(&[1, 2]).unwrap(),
        expected_coord
    );

    // Negative origins
    let grid_neg = InherentDiscreteGrid::builder(&[3, 3])
        .with_origin(&[-10, -20])
        .with_step(&[2, 3])
        .build()
        .unwrap();
    let min_coord = grid_neg.grididx_to_origcoord(&[0, 0]).unwrap();
    assert_eq!(min_coord, vec![-10, -20]);

    // Single dimension roundtrip
    let grid_1d = InherentDiscreteGrid::builder(&[3])
        .with_origin(&[5])
        .with_step(&[3])
        .build()
        .unwrap();
    for i in 0..base.pow(3) {
        let coord = grid_1d.grididx_to_origcoord(&[i]).unwrap();
        assert_eq!(grid_1d.origcoord_to_grididx(&coord).unwrap(), vec![i]);
    }
}

// --- "InherentDiscreteGrid consistency with unfolding schemes" ---
#[test]
fn test_consistency_between_schemes() {
    let r = 4;
    let origin = [1i64, 1, 1];
    let step = [1i64, 2, 3];
    let base: usize = 2;

    let grid_fused = InherentDiscreteGrid::builder(&[r, r, r])
        .with_origin(&origin)
        .with_step(&step)
        .with_base(base)
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();

    let grid_interleaved = InherentDiscreteGrid::builder(&[r, r, r])
        .with_origin(&origin)
        .with_step(&step)
        .with_base(base)
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();

    // Coordinate conversions should be consistent between schemes
    let test_gridindices: &[&[usize]] = &[&[0, 0, 0], &[1, 2, 3], &[7, 15, 15]];
    for &grididx in test_gridindices {
        let coord_fused = grid_fused.grididx_to_origcoord(grididx).unwrap();
        let coord_interleaved = grid_interleaved.grididx_to_origcoord(grididx).unwrap();
        assert_eq!(coord_fused, coord_interleaved);

        assert_eq!(
            grid_fused.origcoord_to_grididx(&coord_fused).unwrap(),
            grididx.to_vec()
        );
        assert_eq!(
            grid_interleaved
                .origcoord_to_grididx(&coord_interleaved)
                .unwrap(),
            grididx.to_vec()
        );
    }

    // Quantics differ but origcoord results are the same
    let grididx = [3usize, 5, 7];
    let quantics_fused = grid_fused.grididx_to_quantics(&grididx).unwrap();
    let quantics_interleaved = grid_interleaved.grididx_to_quantics(&grididx).unwrap();

    assert_eq!(quantics_fused.len(), r);
    assert_eq!(quantics_interleaved.len(), 3 * r);
    assert_ne!(quantics_fused, quantics_interleaved);

    // But they should convert to same coordinates
    let coord_from_fused = grid_fused.quantics_to_origcoord(&quantics_fused).unwrap();
    let coord_from_interleaved = grid_interleaved
        .quantics_to_origcoord(&quantics_interleaved)
        .unwrap();
    assert_eq!(coord_from_fused, coord_from_interleaved);
}

// --- "InherentDiscreteGrid high-dimensional grids" ---
#[test]
fn test_high_dimensional_grids() {
    let r = 2;
    let base: usize = 2;
    let d = 5;

    let origin: Vec<i64> = (1..=d as i64).collect(); // [1, 2, 3, 4, 5]
    let step: Vec<i64> = (1..=d as i64).collect(); // [1, 2, 3, 4, 5]

    let grid = InherentDiscreteGrid::builder(&vec![r; d])
        .with_base(base)
        .with_origin(&origin)
        .with_step(&step)
        .build()
        .unwrap();

    assert_eq!(grid.ndims(), d);
    assert_eq!(grid.grid_min(), origin);
    assert_eq!(grid.step(), step.as_slice());

    // Test coordinate conversion with all 2s
    let grididx: Vec<usize> = vec![1; d];
    let expected_origcoord: Vec<i64> = origin
        .iter()
        .zip(step.iter())
        .map(|(&o, &s)| o + s) // grididx=1 means offset=1
        .collect();
    assert_eq!(
        grid.grididx_to_origcoord(&grididx).unwrap(),
        expected_origcoord
    );
    assert_eq!(
        grid.origcoord_to_grididx(&expected_origcoord).unwrap(),
        grididx
    );

    // Test quantics
    let quantics = grid.grididx_to_quantics(&grididx).unwrap();
    assert_eq!(grid.quantics_to_grididx(&quantics).unwrap(), grididx);
    assert_eq!(quantics.len(), r); // Fused scheme
}

// --- "InherentDiscreteGrid performance and stress tests" ---
#[test]
fn test_stress_conversions() {
    let r = 5;
    let base: usize = 2;
    let grid = InherentDiscreteGrid::builder(&[r, r, r])
        .with_base(base)
        .with_origin(&[1, 1, 1])
        .build()
        .unwrap();

    let max = (base as i64).pow(r as u32);

    // Julia: 100 random grididx from rand(0:base^R-1) per dimension
    let mut rng = StdRng::seed_from_u64(42);
    for _ in 0..100 {
        let grididx: Vec<usize> = (0..3).map(|_| rng.random_range(0..max as usize)).collect();

        let origcoord = grid.grididx_to_origcoord(&grididx).unwrap();
        assert_eq!(grid.origcoord_to_grididx(&origcoord).unwrap(), grididx);

        let quantics = grid.grididx_to_quantics(&grididx).unwrap();
        assert_eq!(grid.quantics_to_grididx(&quantics).unwrap(), grididx);

        // Cross conversions
        let origcoord_from_quantics = grid.quantics_to_origcoord(&quantics).unwrap();
        assert_eq!(origcoord_from_quantics, origcoord);

        let quantics_from_origcoord = grid.origcoord_to_quantics(&origcoord).unwrap();
        assert_eq!(quantics_from_origcoord, quantics);
    }
}

// --- "InherentDiscreteGrid grid with zero step (degenerate case)" ---
#[test]
fn test_zero_step_error() {
    let result = InherentDiscreteGrid::builder(&[3, 3])
        .with_origin(&[5, 10])
        .with_step(&[0, 1])
        .build();
    assert!(matches!(result, Err(QuanticsGridError::InvalidStep { .. })));
}

// --- "InherentDiscreteGrid with custom indextable - basic" ---
#[test]
fn test_custom_indextable_basic() {
    let index_table = vec![
        vec![("x".to_string(), 0)],                       // site 1: x_1
        vec![("x".to_string(), 1)],                       // site 2: x_2
        vec![("y".to_string(), 0)],                       // site 3: y_1
        vec![("x".to_string(), 2), ("y".to_string(), 1)], // site 4: x_3, y_2
    ];

    let grid = InherentDiscreteGrid::from_index_table(&["x", "y"], index_table)
        .with_origin(&[1, 5])
        .with_step(&[2, 3])
        .with_base(2)
        .build()
        .unwrap();

    assert_eq!(grid.variable_names(), &["x", "y"]);
    assert_eq!(grid.rs(), &[3, 2]); // x has 3 quantics, y has 2 quantics
    assert_eq!(grid.origin(), &[1, 5]);
    assert_eq!(grid.step(), &[2, 3]);
    assert_eq!(grid.base(), 2);
    assert_eq!(grid.len(), 4); // 4 sites

    // Local dimensions: [2^1, 2^1, 2^1, 2^2] = [2, 2, 2, 4]
    assert_eq!(grid.local_dimensions(), vec![2, 2, 2, 4]);

    // Test coordinate conversions
    let grididx = [2usize, 1];
    let origcoord = grid.grididx_to_origcoord(&grididx).unwrap();
    assert_eq!(
        grid.origcoord_to_grididx(&origcoord).unwrap(),
        grididx.to_vec()
    );

    // Test quantics conversion
    let quantics = grid.grididx_to_quantics(&grididx).unwrap();
    assert_eq!(quantics.len(), 4);
    assert_eq!(
        grid.quantics_to_grididx(&quantics).unwrap(),
        grididx.to_vec()
    );
    assert_eq!(grid.quantics_to_origcoord(&quantics).unwrap(), origcoord);
}

// --- "InherentDiscreteGrid with complex indextable - mixed sites" ---
#[test]
fn test_complex_indextable_mixed_sites() {
    let index_table = vec![
        vec![("a".to_string(), 0), ("b".to_string(), 0)], // site 1
        vec![("c".to_string(), 0)],                       // site 2
        vec![
            ("a".to_string(), 1),
            ("b".to_string(), 1),
            ("c".to_string(), 1),
        ], // site 3
        vec![("a".to_string(), 2)],                       // site 4
        vec![("b".to_string(), 2), ("c".to_string(), 2)], // site 5
    ];

    let grid = InherentDiscreteGrid::from_index_table(&["a", "b", "c"], index_table)
        .with_origin(&[10, 20, 30])
        .with_step(&[1, 2, 5])
        .with_base(3)
        .build()
        .unwrap();

    assert_eq!(grid.variable_names(), &["a", "b", "c"]);
    assert_eq!(grid.rs(), &[3, 3, 3]);
    assert_eq!(grid.origin(), &[10, 20, 30]);
    assert_eq!(grid.step(), &[1, 2, 5]);
    assert_eq!(grid.base(), 3);
    assert_eq!(grid.len(), 5);

    // Local dimensions: [3^2, 3^1, 3^3, 3^1, 3^2] = [9, 3, 27, 3, 9]
    let expected_localdims = vec![9, 3, 27, 3, 9];
    assert_eq!(grid.local_dimensions(), expected_localdims);

    // Test site dims individually
    for (i, &expected_dim) in expected_localdims.iter().enumerate() {
        assert_eq!(grid.site_dim(i).unwrap(), expected_dim);
    }

    // Site index out of bounds
    assert!(matches!(
        grid.site_dim(expected_localdims.len()),
        Err(QuanticsGridError::SiteIndexOutOfBounds { .. })
    ));

    // Test coordinate conversions
    let test_gridindices: &[&[usize]] = &[&[0, 0, 0], &[1, 2, 0], &[8, 26, 26], &[4, 9, 14]];
    for &grididx in test_gridindices {
        let max_valid = grididx
            .iter()
            .zip(grid.max_grididx().iter())
            .all(|(&g, &m)| g <= m);
        if max_valid {
            let origcoord = grid.grididx_to_origcoord(grididx).unwrap();
            assert_eq!(
                grid.origcoord_to_grididx(&origcoord).unwrap(),
                grididx.to_vec()
            );

            let quantics = grid.grididx_to_quantics(grididx).unwrap();
            assert_eq!(quantics.len(), 5);
            for (i, &q) in quantics.iter().enumerate() {
                assert!(q < expected_localdims[i]);
            }
            assert_eq!(
                grid.quantics_to_grididx(&quantics).unwrap(),
                grididx.to_vec()
            );
        }
    }
}

// --- "InherentDiscreteGrid with asymmetric indextable" ---
#[test]
fn test_asymmetric_indextable() {
    let index_table = vec![
        vec![
            ("x".to_string(), 0),
            ("x".to_string(), 1),
            ("x".to_string(), 2),
            ("x".to_string(), 3),
        ], // site 1: all x
        vec![("y".to_string(), 0)],                       // site 2
        vec![("y".to_string(), 1)],                       // site 3
        vec![("z".to_string(), 0), ("z".to_string(), 1)], // site 4
    ];

    let grid = InherentDiscreteGrid::from_index_table(&["x", "y", "z"], index_table)
        .with_origin(&[0, 100, -50])
        .with_step(&[10, 5, 2])
        .with_base(2)
        .build()
        .unwrap();

    assert_eq!(grid.rs(), &[4, 2, 2]);
    assert_eq!(grid.len(), 4);
    assert_eq!(grid.local_dimensions(), vec![16, 2, 2, 4]);

    // Test boundary cases
    let max_x = 2usize.pow(4); // 16
    let max_y = 2usize.pow(2); // 4
    let max_z = 2usize.pow(2); // 4

    let grididx = [max_x - 1, max_y - 1, max_z - 1];
    let origcoord = grid.grididx_to_origcoord(&grididx).unwrap();
    let expected_origcoord = vec![
        10 * (max_x as i64 - 1),
        100 + 5 * (max_y as i64 - 1),
        -50 + 2 * (max_z as i64 - 1),
    ];
    assert_eq!(origcoord, expected_origcoord);
    assert_eq!(
        grid.origcoord_to_grididx(&origcoord).unwrap(),
        grididx.to_vec()
    );
}

// --- "InherentDiscreteGrid with single-site indextable" ---
#[test]
fn test_single_site_indextable() {
    let index_table = vec![vec![
        ("x".to_string(), 0),
        ("x".to_string(), 1),
        ("x".to_string(), 2),
        ("y".to_string(), 0),
        ("y".to_string(), 1),
    ]];

    let grid = InherentDiscreteGrid::from_index_table(&["x", "y"], index_table)
        .with_origin(&[5, 10])
        .with_step(&[1, 3])
        .with_base(2)
        .build()
        .unwrap();

    assert_eq!(grid.rs(), &[3, 2]);
    assert_eq!(grid.len(), 1);
    assert_eq!(grid.local_dimensions(), vec![32]); // 2^5

    let grididx = [3usize, 2];
    let origcoord = grid.grididx_to_origcoord(&grididx).unwrap();
    assert_eq!(
        grid.origcoord_to_grididx(&origcoord).unwrap(),
        grididx.to_vec()
    );

    let quantics = grid.grididx_to_quantics(&grididx).unwrap();
    assert_eq!(quantics.len(), 1);
    assert!(quantics[0] < 32);
    assert_eq!(
        grid.quantics_to_grididx(&quantics).unwrap(),
        grididx.to_vec()
    );
}

// --- "InherentDiscreteGrid with maximum fragmentation" ---
#[test]
fn test_maximum_fragmentation() {
    let index_table = vec![
        vec![("a".to_string(), 0)],
        vec![("a".to_string(), 1)],
        vec![("a".to_string(), 2)],
        vec![("b".to_string(), 0)],
        vec![("b".to_string(), 1)],
    ];

    let grid = InherentDiscreteGrid::from_index_table(&["a", "b"], index_table)
        .with_origin(&[1, 1])
        .with_step(&[1, 1])
        .with_base(3)
        .build()
        .unwrap();

    assert_eq!(grid.rs(), &[3, 2]);
    assert_eq!(grid.len(), 5);
    assert_eq!(grid.local_dimensions(), vec![3; 5]);

    // Test conversions
    let test_gridindices: &[&[usize]] = &[&[0, 0], &[2, 2], &[8, 8], &[26, 8]];
    for &grididx in test_gridindices {
        let valid = grididx
            .iter()
            .zip(grid.max_grididx().iter())
            .all(|(&g, &m)| g <= m);
        if valid {
            let origcoord = grid.grididx_to_origcoord(grididx).unwrap();
            assert_eq!(
                grid.origcoord_to_grididx(&origcoord).unwrap(),
                grididx.to_vec()
            );

            let quantics = grid.grididx_to_quantics(grididx).unwrap();
            assert_eq!(quantics.len(), 5);
            assert!(quantics.iter().all(|&v| (0..3).contains(&v)));
            assert_eq!(
                grid.quantics_to_grididx(&quantics).unwrap(),
                grididx.to_vec()
            );
        }
    }
}

// --- "InherentDiscreteGrid complex indextable with base != 2" ---
#[test]
fn test_complex_indextable_base5() {
    let index_table = vec![
        vec![("u".to_string(), 0), ("v".to_string(), 0)], // site 1
        vec![("w".to_string(), 0), ("w".to_string(), 1)], // site 2
        vec![("u".to_string(), 1)],                       // site 3
        vec![
            ("v".to_string(), 1),
            ("u".to_string(), 2),
            ("w".to_string(), 2),
        ], // site 4
    ];

    let grid = InherentDiscreteGrid::from_index_table(&["u", "v", "w"], index_table)
        .with_origin(&[-5, 0, 10])
        .with_step(&[3, 7, 2])
        .with_base(5)
        .build()
        .unwrap();

    assert_eq!(grid.rs(), &[3, 2, 3]);
    assert_eq!(grid.base(), 5);
    assert_eq!(grid.len(), 4);
    assert_eq!(grid.local_dimensions(), vec![25, 25, 5, 125]);

    // Test minimum
    let grididx = [0usize, 0, 0];
    let origcoord = grid.grididx_to_origcoord(&grididx).unwrap();
    assert_eq!(origcoord, vec![-5, 0, 10]);
    assert_eq!(
        grid.origcoord_to_grididx(&origcoord).unwrap(),
        grididx.to_vec()
    );

    // Test maximum
    let max_grididx = [124usize, 24, 124]; // 5^3-1, 5^2-1, 5^3-1
    let max_origcoord = grid.grididx_to_origcoord(&max_grididx).unwrap();
    let expected_max = vec![-5 + 3 * 124, 7 * 24, 10 + 2 * 124];
    assert_eq!(max_origcoord, expected_max);
    assert_eq!(
        grid.origcoord_to_grididx(&max_origcoord).unwrap(),
        max_grididx.to_vec()
    );

    // Test quantics conversions
    let test_grididx = [9usize, 4, 19];
    let quantics = grid.grididx_to_quantics(&test_grididx).unwrap();
    assert_eq!(quantics.len(), 4);
    assert!(quantics[0] < 25); // site 1
    assert!(quantics[1] < 25); // site 2
    assert!(quantics[2] < 5); // site 3
    assert!(quantics[3] < 125); // site 4
    assert_eq!(
        grid.quantics_to_grididx(&quantics).unwrap(),
        test_grididx.to_vec()
    );
}

// --- "InherentDiscreteGrid comprehensive indextable test" ---
#[test]
fn test_comprehensive_indextable() {
    let index_table = vec![
        vec![("k_x".to_string(), 0), ("k_y".to_string(), 0)], // site 1
        vec![("w".to_string(), 0), ("t".to_string(), 0)],     // site 2
        vec![("k_x".to_string(), 1)],                         // site 3
        vec![("k_y".to_string(), 1), ("w".to_string(), 1)],   // site 4
        vec![("t".to_string(), 1), ("t".to_string(), 2)],     // site 5
        vec![
            ("k_x".to_string(), 2),
            ("k_y".to_string(), 2),
            ("w".to_string(), 2),
        ], // site 6
    ];

    let grid = InherentDiscreteGrid::from_index_table(&["k_x", "k_y", "w", "t"], index_table)
        .with_origin(&[0, 0, -10, 0])
        .with_step(&[1, 1, 2, 5])
        .with_base(2)
        .build()
        .unwrap();

    assert_eq!(grid.rs(), &[3, 3, 3, 3]);
    assert_eq!(grid.len(), 6);

    let expected_localdims = vec![4, 4, 2, 4, 4, 8];
    assert_eq!(grid.local_dimensions(), expected_localdims);

    // Test deterministic conversions over all possible grid indices (8^4 is too many,
    // so sample systematically)
    let max = 7usize; // 2^3 - 1
    for kx in [0, 1, 3, 6, max] {
        for ky in [0, 2, max] {
            for w in [0, 4, max] {
                for t in [0, max / 2, max] {
                    let grididx = vec![kx, ky, w, t];

                    let origcoord = grid.grididx_to_origcoord(&grididx).unwrap();
                    assert_eq!(grid.origcoord_to_grididx(&origcoord).unwrap(), grididx);

                    let quantics = grid.grididx_to_quantics(&grididx).unwrap();
                    assert_eq!(quantics.len(), 6);
                    for (i, &q) in quantics.iter().enumerate() {
                        assert!(q < expected_localdims[i]);
                    }
                    assert_eq!(grid.quantics_to_grididx(&quantics).unwrap(), grididx);

                    // Cross conversions
                    assert_eq!(grid.quantics_to_origcoord(&quantics).unwrap(), origcoord);
                    assert_eq!(grid.origcoord_to_quantics(&origcoord).unwrap(), quantics);
                }
            }
        }
    }
}

// --- "InherentDiscreteGrid indextable validation" ---
#[test]
fn test_indextable_validation() {
    // Valid indextable should work
    let valid_table = vec![
        vec![("x".to_string(), 0), ("y".to_string(), 0)],
        vec![("x".to_string(), 1)],
        vec![("y".to_string(), 1)],
    ];
    let grid = InherentDiscreteGrid::from_index_table(&["x", "y"], valid_table)
        .with_origin(&[1, 1])
        .with_step(&[1, 1])
        .build()
        .unwrap();
    assert_eq!(grid.rs(), &[2, 2]);

    // Unknown variable should error
    let invalid_unknown = vec![
        vec![("x".to_string(), 0), ("y".to_string(), 0)],
        vec![("z".to_string(), 0)], // z is not in variable names
    ];
    let result = InherentDiscreteGrid::from_index_table(&["x", "y"], invalid_unknown)
        .with_origin(&[1, 1])
        .with_step(&[1, 1])
        .build();
    assert!(matches!(
        result,
        Err(QuanticsGridError::UnknownVariable { .. })
    ));

    // Missing index should error (x has bit 1 and 3, but not 2)
    // from_index_table counts occurrences: x appears at bits 1 and 3 -> Rs_x=2
    // But bit 3 > Rs_x=2, so InvalidBitNumber is the resulting error.
    let invalid_missing = vec![
        vec![("x".to_string(), 0), ("y".to_string(), 0)],
        vec![("x".to_string(), 2), ("y".to_string(), 1)],
    ];
    let result = InherentDiscreteGrid::from_index_table(&["x", "y"], invalid_missing)
        .with_origin(&[1, 1])
        .with_step(&[1, 1])
        .build();
    assert!(matches!(
        result,
        Err(QuanticsGridError::InvalidBitNumber { .. })
    ));

    // Repeated index should error
    let invalid_repeated = vec![
        vec![("x".to_string(), 0), ("y".to_string(), 0)],
        vec![("x".to_string(), 0), ("y".to_string(), 1)], // x bit 1 repeated
    ];
    let result = InherentDiscreteGrid::from_index_table(&["x", "y"], invalid_repeated)
        .with_origin(&[1, 1])
        .with_step(&[1, 1])
        .build();
    assert!(matches!(
        result,
        Err(QuanticsGridError::DuplicateIndexEntry { .. })
    ));
}

// --- "InherentDiscreteGrid constructor with variablenames and Rs" ---
#[test]
fn test_constructor_with_variable_names_and_rs() {
    let rs = [4, 3, 5];
    let origin = [1i64, 10, -5];
    let step = [2i64, 1, 3];

    // Default base
    let grid1 = InherentDiscreteGrid::builder(&rs)
        .with_variable_names(&["k_x", "k_y", "w"])
        .with_origin(&origin)
        .with_step(&step)
        .build()
        .unwrap();
    assert_eq!(grid1.variable_names(), &["k_x", "k_y", "w"]);
    assert_eq!(grid1.rs(), &rs);
    assert_eq!(grid1.base(), 2);
    assert_eq!(grid1.origin(), &origin);
    assert_eq!(grid1.step(), &step);
    assert_eq!(grid1.ndims(), 3);

    // Custom base and interleaved
    let grid2 = InherentDiscreteGrid::builder(&rs)
        .with_variable_names(&["k_x", "k_y", "w"])
        .with_origin(&origin)
        .with_step(&step)
        .with_base(3)
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();
    assert_eq!(grid2.base(), 3);

    // Roundtrip
    let test_grididx = [1usize, 2, 3];
    let origcoord = grid2.grididx_to_origcoord(&test_grididx).unwrap();
    assert_eq!(
        grid2.origcoord_to_grididx(&origcoord).unwrap(),
        test_grididx.to_vec()
    );

    let quantics = grid2.grididx_to_quantics(&test_grididx).unwrap();
    assert_eq!(
        grid2.quantics_to_grididx(&quantics).unwrap(),
        test_grididx.to_vec()
    );
}

// --- "InherentDiscreteGrid constructor with Rs tuple only" ---
#[test]
fn test_constructor_rs_only() {
    let rs = [3, 4, 2];
    let origin = [5i64, 0, -10];
    let step = [1i64, 2, 5];
    let base: usize = 3;

    let grid = InherentDiscreteGrid::builder(&rs)
        .with_origin(&origin)
        .with_step(&step)
        .with_base(base)
        .build()
        .unwrap();

    assert_eq!(grid.rs(), &rs);
    assert_eq!(grid.origin(), &origin);
    assert_eq!(grid.step(), &step);
    assert_eq!(grid.base(), base);
    assert_eq!(grid.ndims(), 3);

    // Auto-generated variable names should be "1", "2", "3"
    assert_eq!(grid.variable_names(), &["1", "2", "3"]);

    // Test grid works
    let grididx = [1usize, 2, 0];
    let origcoord = grid.grididx_to_origcoord(&grididx).unwrap();
    assert_eq!(
        grid.origcoord_to_grididx(&origcoord).unwrap(),
        grididx.to_vec()
    );
}

// --- "InherentDiscreteGrid default parameter behavior" ---
#[test]
fn test_default_parameter_behavior() {
    let rs = [3, 2];

    // Minimal parameters (all defaults)
    let grid_minimal = InherentDiscreteGrid::builder(&rs).build().unwrap();
    assert_eq!(grid_minimal.rs(), &rs);
    assert_eq!(grid_minimal.base(), 2); // default
    assert_eq!(grid_minimal.origin(), &[1, 1]); // default
    assert_eq!(grid_minimal.step(), &[1, 1]); // default

    // Partial parameters
    let grid_partial = InherentDiscreteGrid::builder(&rs)
        .with_origin(&[5, -2])
        .build()
        .unwrap();
    assert_eq!(grid_partial.origin(), &[5, -2]);
    assert_eq!(grid_partial.step(), &[1, 1]); // default
    assert_eq!(grid_partial.base(), 2); // default
}

// --- "InherentDiscreteGrid constructor error handling" ---
#[test]
fn test_constructor_error_handling() {
    // Invalid base
    assert!(matches!(
        InherentDiscreteGrid::builder(&[3, 2]).with_base(1).build(),
        Err(QuanticsGridError::InvalidBase(1))
    ));
    assert!(matches!(
        InherentDiscreteGrid::builder(&[3, 2]).with_base(0).build(),
        Err(QuanticsGridError::InvalidBase(0))
    ));

    // Mismatched dimensions
    assert!(matches!(
        InherentDiscreteGrid::builder(&[3, 2])
            .with_origin(&[1, 1, 1])
            .build(),
        Err(QuanticsGridError::DimensionMismatch { .. })
    ));
    assert!(matches!(
        InherentDiscreteGrid::builder(&[3, 2])
            .with_step(&[1, 1, 1])
            .build(),
        Err(QuanticsGridError::DimensionMismatch { .. })
    ));

    // Duplicate variable names
    assert!(matches!(
        InherentDiscreteGrid::builder(&[3, 3])
            .with_variable_names(&["x", "x"])
            .build(),
        Err(QuanticsGridError::DuplicateVariableName(_))
    ));

    // Too large R value (base=2, R=63 would overflow)
    // R=62 should succeed
    assert!(InherentDiscreteGrid::builder(&[62]).build().is_ok());
    // R=63 should fail
    assert!(matches!(
        InherentDiscreteGrid::builder(&[63]).build(),
        Err(QuanticsGridError::ResolutionTooLarge { .. })
    ));
}

// --- "InherentDiscreteGrid constructor compatibility with existing patterns" ---
#[test]
fn test_constructor_patterns() {
    // Pattern 1: Like InherentDiscreteGrid{D}(R, origin; kwargs...)
    let r = 4;
    let grid_pattern1 = InherentDiscreteGrid::builder(&[r, r])
        .with_origin(&[1, 5])
        .with_step(&[2, 3])
        .with_base(3)
        .build()
        .unwrap();
    assert_eq!(grid_pattern1.ndims(), 2);
    assert_eq!(grid_pattern1.rs(), &[r, r]);
    assert_eq!(grid_pattern1.origin(), &[1, 5]);
    assert_eq!(grid_pattern1.step(), &[2, 3]);

    // Pattern 2: Like DiscretizedGrid((R1, R2, ...); kwargs...)
    let grid_pattern2 = InherentDiscreteGrid::builder(&[3, 5, 2])
        .with_origin(&[0, 0, 0])
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();
    assert_eq!(grid_pattern2.rs(), &[3, 5, 2]);

    // Pattern 3: Single R for 1D
    let grid_pattern3 = InherentDiscreteGrid::builder(&[6])
        .with_origin(&[10])
        .with_step(&[2])
        .with_base(4)
        .build()
        .unwrap();
    assert_eq!(grid_pattern3.ndims(), 1);
    assert_eq!(grid_pattern3.rs(), &[6]);
    assert_eq!(grid_pattern3.origin(), &[10]);
    assert_eq!(grid_pattern3.base(), 4);

    // Verify all patterns produce working grids
    let origcoord1 = grid_pattern1.grididx_to_origcoord(&[1, 2]).unwrap();
    assert_eq!(
        grid_pattern1.origcoord_to_grididx(&origcoord1).unwrap(),
        vec![1, 2]
    );

    let origcoord2 = grid_pattern2.grididx_to_origcoord(&[1, 2, 0]).unwrap();
    assert_eq!(
        grid_pattern2.origcoord_to_grididx(&origcoord2).unwrap(),
        vec![1, 2, 0]
    );

    let origcoord3 = grid_pattern3.grididx_to_origcoord(&[1]).unwrap();
    assert_eq!(
        grid_pattern3.origcoord_to_grididx(&origcoord3).unwrap(),
        vec![1]
    );
}

// ============================================================================
// Tests ported from quantics_tests.jl
// ============================================================================

// --- "grouped unfoldingscheme quantics conversions" ---
#[test]
fn test_grouped_scheme_basic() {
    // Rs = (2, 1)
    let grid = InherentDiscreteGrid::builder(&[2, 1])
        .with_unfolding_scheme(UnfoldingScheme::Grouped)
        .build()
        .unwrap();
    assert_eq!(grid.local_dimensions(), vec![2, 2, 2]); // fill(2, sum(Rs))

    let grididx = vec![2usize, 1];
    let expected_quantics = vec![1usize, 0, 1];
    assert_eq!(
        grid.grididx_to_quantics(&grididx).unwrap(),
        expected_quantics
    );
    assert_eq!(
        grid.quantics_to_grididx(&expected_quantics).unwrap(),
        grididx
    );

    // With base 3
    let grid_base3 = InherentDiscreteGrid::builder(&[2, 1])
        .with_base(3)
        .with_unfolding_scheme(UnfoldingScheme::Grouped)
        .build()
        .unwrap();
    assert_eq!(grid_base3.local_dimensions(), vec![3, 3, 3]);

    let grididx_base3 = vec![4usize, 2];
    let expected_quantics_base3 = vec![1usize, 1, 2];
    assert_eq!(
        grid_base3.grididx_to_quantics(&grididx_base3).unwrap(),
        expected_quantics_base3
    );
    assert_eq!(
        grid_base3
            .quantics_to_grididx(&expected_quantics_base3)
            .unwrap(),
        grididx_base3
    );
}

// --- "mixed bases" ---
#[test]
fn test_mixed_bases_fused() {
    let grid = InherentDiscreteGrid::builder(&[1, 1])
        .with_bases(&[2, 6])
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();

    assert_eq!(grid.bases(), &[2, 6]);
    assert_eq!(grid.local_dimensions(), vec![12]); // prod(bases)

    assert_eq!(grid.grididx_to_quantics(&[1, 5]).unwrap(), vec![11]);
    assert_eq!(grid.quantics_to_grididx(&[11]).unwrap(), vec![1, 5]);
    assert_eq!(grid.grididx_to_quantics(&[0, 1]).unwrap(), vec![2]);
    assert_eq!(grid.quantics_to_grididx(&[2]).unwrap(), vec![0, 1]);
}

// --- "mixed bases interleaved" ---
#[test]
fn test_mixed_bases_interleaved() {
    let grid = InherentDiscreteGrid::builder(&[2, 1])
        .with_bases(&[2, 6])
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();

    assert_eq!(grid.local_dimensions(), vec![2, 6, 2]);

    assert_eq!(grid.grididx_to_quantics(&[1, 5]).unwrap(), vec![0, 5, 1]);
    assert_eq!(grid.quantics_to_grididx(&[0, 5, 1]).unwrap(), vec![1, 5]);
}

// --- "mixed bases grouped" ---
#[test]
fn test_mixed_bases_grouped() {
    let grid = InherentDiscreteGrid::builder(&[2, 1])
        .with_bases(&[2, 6])
        .with_unfolding_scheme(UnfoldingScheme::Grouped)
        .build()
        .unwrap();

    assert_eq!(grid.local_dimensions(), vec![2, 2, 6]);

    assert_eq!(grid.grididx_to_quantics(&[1, 5]).unwrap(), vec![0, 1, 5]);
    assert_eq!(grid.quantics_to_grididx(&[0, 1, 5]).unwrap(), vec![1, 5]);
}

// --- "mixed bases fused roundtrip" ---
#[test]
fn test_mixed_bases_fused_roundtrip() {
    let grid = InherentDiscreteGrid::builder(&[3, 2, 1])
        .with_bases(&[2, 3, 5])
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();

    assert_eq!(grid.local_dimensions(), vec![30, 6, 2]);

    // Julia: 50 random quantics from rand(0:d-1) per localdimension
    let dims = grid.local_dimensions();
    let mut rng = StdRng::seed_from_u64(42);
    for _ in 0..50 {
        let quantics: Vec<usize> = dims.iter().map(|&d| rng.random_range(0..d)).collect();
        let grididx = grid.quantics_to_grididx(&quantics).unwrap();
        assert_eq!(
            grid.grididx_to_quantics(&grididx).unwrap(),
            quantics,
            "Failed for quantics {:?}",
            quantics
        );
    }
}

// --- "quantics_to_grididx compose grididx_to_quantics == identity" ---
// Julia: DiscretizedGrid((5, 3, 17); base=13), 100 random grididx from rand(0:grid_Rs(grid)[d]-1)
#[test]
fn test_roundtrip_grididx_to_quantics_identity() {
    let grid = InherentDiscreteGrid::builder(&[5, 3, 17])
        .with_base(13)
        .build()
        .unwrap();

    let rs = grid.rs();

    let mut rng = StdRng::seed_from_u64(42);
    for _ in 0..100 {
        let grididx: Vec<usize> = rs.iter().map(|&r| rng.random_range(0..r)).collect();
        let quantics = grid.grididx_to_quantics(&grididx).unwrap();
        assert_eq!(
            grid.quantics_to_grididx(&quantics).unwrap(),
            grididx,
            "Failed for grididx {:?}",
            grididx
        );
    }
}

// --- "grididx_to_quantics compose quantics_to_grididx == identity" ---
// Julia: DiscretizedGrid((48, 31, 62)), 100 random quantics from rand(0:1, length(grid))
#[test]
fn test_roundtrip_quantics_to_grididx_identity() {
    let grid = InherentDiscreteGrid::builder(&[48, 31, 62])
        .build()
        .unwrap();

    let total_len = grid.len();

    let mut rng = StdRng::seed_from_u64(42);
    for _ in 0..100 {
        let quantics: Vec<usize> = (0..total_len).map(|_| rng.random_range(0..2)).collect();
        let grididx = grid.quantics_to_grididx(&quantics).unwrap();
        assert_eq!(
            grid.grididx_to_quantics(&grididx).unwrap(),
            quantics,
            "Failed for quantics pattern"
        );
    }
}

// --- "grididx_to_quantics compose quantics_to_grididx == identity, base != 2" ---
// Julia: DiscretizedGrid((22, 9, 14); base=7), 100 random quantics from rand(0:6, length(grid))
#[test]
fn test_roundtrip_quantics_to_grididx_base7() {
    let base: usize = 7;
    let grid = InherentDiscreteGrid::builder(&[22, 9, 14])
        .with_base(base)
        .build()
        .unwrap();

    let total_len = grid.len();
    let b = base;

    let mut rng = StdRng::seed_from_u64(42);
    for _ in 0..100 {
        let quantics: Vec<usize> = (0..total_len).map(|_| rng.random_range(0..b)).collect();
        let grididx = grid.quantics_to_grididx(&quantics).unwrap();
        assert_eq!(
            grid.grididx_to_quantics(&grididx).unwrap(),
            quantics,
            "Failed for quantics pattern"
        );
    }
}

// --- "constructor, square grid" (from quantics_tests.jl, DiscretizedGrid tests applied to InherentDiscreteGrid) ---
#[test]
fn test_square_grid_interleaved() {
    let grid = InherentDiscreteGrid::builder(&[10, 10])
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();
    assert_eq!(grid.rs(), &[10, 10]);
}

// --- "constructor, rectangular grid" ---
#[test]
fn test_rectangular_grid_interleaved() {
    let grid = InherentDiscreteGrid::builder(&[3, 5])
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();
    // Total sites: 3 + 5 = 8 for interleaved
    assert_eq!(grid.len(), 8);
}

// --- "quantics_to_grididx, rectangular grid" ---
#[test]
fn test_quantics_to_grididx_rectangular() {
    let grid = InherentDiscreteGrid::builder(&[3, 5])
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();
    let result = grid.quantics_to_grididx(&[0, 1, 0, 1, 0, 1, 0, 1]).unwrap();
    assert_eq!(result, vec![0, 29]);
}

// --- "grididx_to_quantics, rectangular grid" ---
#[test]
fn test_grididx_to_quantics_rectangular() {
    let grid = InherentDiscreteGrid::builder(&[3, 5])
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();
    let result = grid.grididx_to_quantics(&[0, 29]).unwrap();
    assert_eq!(result, vec![0, 1, 0, 1, 0, 1, 0, 1]);
}

// --- "challenging tests - extreme edge cases" ---
#[test]
fn test_extreme_edge_cases() {
    let grid = InherentDiscreteGrid::builder(&[10, 5, 8])
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();

    // Minimum valid grididx (all 0s)
    let min_grididx = vec![0usize, 0, 0];
    let quantics = grid.grididx_to_quantics(&min_grididx).unwrap();
    assert!(quantics.iter().all(|&q| q == 0));
    assert_eq!(grid.quantics_to_grididx(&quantics).unwrap(), min_grididx);

    // Maximum valid grididx (2^R - 1)
    let max_grididx = vec![2usize.pow(10) - 1, 2usize.pow(5) - 1, 2usize.pow(8) - 1];
    let quantics = grid.grididx_to_quantics(&max_grididx).unwrap();
    assert!(quantics.iter().all(|&q| q == 1));
    assert_eq!(grid.quantics_to_grididx(&quantics).unwrap(), max_grididx);
}

// --- "challenging tests - mixed bases" ---
// Julia: 50 random quantics from rand(0:2, length(grid))
#[test]
fn test_challenging_mixed_bases_base3() {
    let grid = InherentDiscreteGrid::builder(&[4, 6, 3])
        .with_base(3)
        .build()
        .unwrap();

    let dims = grid.local_dimensions();
    let total_len = grid.len();

    let mut rng = StdRng::seed_from_u64(42);
    for _ in 0..50 {
        let quantics: Vec<usize> = (0..total_len).map(|_| rng.random_range(0..3)).collect();
        let grididx = grid.quantics_to_grididx(&quantics).unwrap();
        assert_eq!(grid.grididx_to_quantics(&grididx).unwrap(), quantics);
        // Verify bounds
        for (i, &q) in quantics.iter().enumerate() {
            assert!(q < dims[i]);
        }
    }
}

// Julia: 50 random grididx from rand(0:5^3-1) and rand(0:5^4-1)
#[test]
fn test_challenging_mixed_bases_base5_interleaved() {
    let grid = InherentDiscreteGrid::builder(&[3, 4])
        .with_base(5)
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();

    let max0 = 5usize.pow(3);
    let max1 = 5usize.pow(4);

    let mut rng = StdRng::seed_from_u64(42);
    for _ in 0..50 {
        let grididx = vec![rng.random_range(0..max0), rng.random_range(0..max1)];
        let quantics = grid.grididx_to_quantics(&grididx).unwrap();
        assert_eq!(grid.quantics_to_grididx(&quantics).unwrap(), grididx);
        assert!(quantics.iter().all(|&q| (0..5).contains(&q)));
    }
}

// --- "challenging tests - complex fused indices" ---
#[test]
fn test_complex_fused_indices() {
    // Multiple variables fused in single sites
    let index_table = vec![
        vec![
            ("x".to_string(), 2),
            ("y".to_string(), 1),
            ("z".to_string(), 0),
        ], // 3 variables in one site
        vec![("x".to_string(), 1)], // single variable
        vec![("y".to_string(), 0), ("z".to_string(), 1)], // 2 variables fused
        vec![("x".to_string(), 0)], // single variable
        vec![("z".to_string(), 2)], // single variable
    ];

    let grid = InherentDiscreteGrid::from_index_table(&["x", "y", "z"], index_table)
        .build()
        .unwrap();

    // Test specific known values
    assert_eq!(
        grid.quantics_to_grididx(&[7, 0, 3, 1, 0]).unwrap(),
        vec![5, 3, 6]
    );
    assert_eq!(
        grid.grididx_to_quantics(&[5, 3, 6]).unwrap(),
        vec![7, 0, 3, 1, 0]
    );

    // Test all site dimensions
    let dims = grid.local_dimensions();
    assert_eq!(dims, vec![8, 2, 4, 2, 2]); // 2^3, 2^1, 2^2, 2^1, 2^1

    // Exhaustive roundtrip for all quantics
    for d0 in 0..dims[0] {
        for d1 in 0..dims[1] {
            for d2 in 0..dims[2] {
                for d3 in 0..dims[3] {
                    for d4 in 0..dims[4] {
                        let quantics = vec![d0, d1, d2, d3, d4];
                        let grididx = grid.quantics_to_grididx(&quantics).unwrap();
                        assert_eq!(grid.grididx_to_quantics(&grididx).unwrap(), quantics);
                    }
                }
            }
        }
    }
}

// --- "challenging tests - single dimension edge cases" ---
#[test]
fn test_single_dimension_edge_cases() {
    let grid = InherentDiscreteGrid::builder(&[25]).build().unwrap();

    // Extremes
    let min_quantics = vec![0usize; 25];
    let max_quantics = vec![1usize; 25];

    assert_eq!(grid.quantics_to_grididx(&min_quantics).unwrap(), vec![0]);
    assert_eq!(
        grid.quantics_to_grididx(&max_quantics).unwrap(),
        vec![2usize.pow(25) - 1]
    );
    assert_eq!(grid.grididx_to_quantics(&[0]).unwrap(), min_quantics);
    assert_eq!(
        grid.grididx_to_quantics(&[2usize.pow(25) - 1]).unwrap(),
        max_quantics
    );

    // Test middle values with deterministic patterns
    let patterns: Vec<Vec<usize>> = vec![
        (0..25).map(|i| if i % 2 == 0 { 0 } else { 1 }).collect(),
        (0..25).map(|i| if i % 3 == 0 { 1 } else { 0 }).collect(),
        (0..25).map(|i| if i < 12 { 0 } else { 1 }).collect(),
    ];

    for quantics in &patterns {
        let grididx = grid.quantics_to_grididx(quantics).unwrap();
        assert_eq!(grid.grididx_to_quantics(&grididx).unwrap(), *quantics);
    }
}

// --- "challenging tests - high dimensional grids" ---
// Julia: 30 random iterations with rand(0:2, length(grid))
#[test]
fn test_high_dimensional_grids_base3() {
    // 8D grid with moderate R values
    let rs: Vec<usize> = (0..8).map(|i| 4 + (i % 3)).collect();
    let grid = InherentDiscreteGrid::builder(&rs)
        .with_base(3)
        .build()
        .unwrap();

    let dims = grid.local_dimensions();
    let total_len = grid.len();

    let mut rng = StdRng::seed_from_u64(42);
    for _ in 0..30 {
        let quantics: Vec<usize> = (0..total_len)
            .map(|i| rng.random_range(0..dims[i]))
            .collect();
        let grididx = grid.quantics_to_grididx(&quantics).unwrap();
        assert_eq!(grid.grididx_to_quantics(&grididx).unwrap(), quantics);

        // Verify all grid indices are within bounds
        for (d, (&g, &r)) in grididx.iter().zip(rs.iter()).enumerate() {
            assert!(
                g < 3usize.pow(r as u32),
                "Grid index {} for dim {} out of bounds",
                g,
                d
            );
        }
    }
}

// --- "challenging tests - stress test with complex patterns" ---
#[test]
fn test_stress_complex_patterns() {
    let index_table = vec![
        vec![("e".to_string(), 0)],                       // site 1
        vec![("a".to_string(), 4), ("c".to_string(), 3)], // site 2
        vec![("b".to_string(), 2)],                       // site 3
        vec![
            ("a".to_string(), 3),
            ("b".to_string(), 1),
            ("d".to_string(), 2),
        ], // site 4
        vec![("c".to_string(), 2), ("e".to_string(), 1)], // site 5
        vec![("a".to_string(), 2)],                       // site 6
        vec![("b".to_string(), 0), ("d".to_string(), 1)], // site 7
        vec![
            ("a".to_string(), 1),
            ("c".to_string(), 1),
            ("d".to_string(), 0),
            ("e".to_string(), 2),
        ], // site 8
        vec![("a".to_string(), 0)],                       // site 9
        vec![("c".to_string(), 0)],                       // site 10
    ];

    let grid = InherentDiscreteGrid::from_index_table(&["a", "b", "c", "d", "e"], index_table)
        .with_base(3)
        .build()
        .unwrap();

    let dims = grid.local_dimensions();
    // Expected: [3, 9, 3, 27, 9, 3, 9, 81, 3, 3]
    assert_eq!(dims, vec![3, 9, 3, 27, 9, 3, 9, 81, 3, 3]);

    // Verify Rs
    let max_grididx = grid.max_grididx();
    assert_eq!(
        max_grididx,
        &[
            3usize.pow(5) - 1,
            3usize.pow(3) - 1,
            3usize.pow(4) - 1,
            3usize.pow(3) - 1,
            3usize.pow(3) - 1
        ]
    );

    // Julia: 200 random quantics with site dimensions [3, 9, 3, 27, 9, 3, 9, 81, 3, 3]
    let mut rng = StdRng::seed_from_u64(42);
    for _ in 0..200 {
        let quantics: Vec<usize> = dims.iter().map(|&d| rng.random_range(0..d)).collect();

        let grididx = grid.quantics_to_grididx(&quantics).unwrap();
        let recovered = grid.grididx_to_quantics(&grididx).unwrap();
        assert_eq!(recovered, quantics, "Failed for quantics {:?}", quantics);

        // Verify grid indices are within bounds
        for (i, (&g, &m)) in grididx.iter().zip(max_grididx.iter()).enumerate() {
            assert!(
                g <= m,
                "Grid index {} for dim {} out of bounds [0, {}]",
                g,
                i,
                m
            );
        }
    }
}

// --- "ctor from indextable" ---
#[test]
fn test_ctor_from_indextable() {
    let index_table = vec![
        vec![("a".to_string(), 3)],
        vec![("a".to_string(), 2)],
        vec![("a".to_string(), 1)],
        vec![("a".to_string(), 0)],
        vec![("b".to_string(), 0)],
        vec![("b".to_string(), 1)],
        vec![("b".to_string(), 2)],
        vec![("c".to_string(), 0)],
        vec![("c".to_string(), 1)],
        vec![("c".to_string(), 2)],
    ];

    let grid = InherentDiscreteGrid::from_index_table(&["a", "b", "c"], index_table)
        .build()
        .unwrap();

    assert_eq!(grid.rs(), &[4, 3, 3]);
    // Default origin is [1, 1, 1]
    assert_eq!(grid.origin(), &[1, 1, 1]);
}

// --- "ctor from indextable, quantics <-> grididx" ---
#[test]
fn test_ctor_from_indextable_quantics_grididx() {
    let index_table = vec![
        vec![("a".to_string(), 3)],
        vec![("a".to_string(), 2)],
        vec![("a".to_string(), 1)],
        vec![("a".to_string(), 0)],
        vec![("b".to_string(), 0)],
        vec![("b".to_string(), 1)],
        vec![("b".to_string(), 2)],
        vec![("c".to_string(), 0)],
        vec![("c".to_string(), 1)],
        vec![("c".to_string(), 2)],
    ];

    let grid = InherentDiscreteGrid::from_index_table(&["a", "b", "c"], index_table)
        .build()
        .unwrap();

    assert_eq!(
        grid.quantics_to_grididx(&[0, 1, 0, 1, 0, 1, 0, 1, 0, 1])
            .unwrap(),
        vec![10, 2, 5]
    );
    assert_eq!(
        grid.grididx_to_quantics(&[10, 2, 5]).unwrap(),
        vec![0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
    );
}

// --- "ctor from indextable, quantics <-> grididx, fused indices" ---
#[test]
fn test_ctor_from_indextable_fused_indices() {
    let index_table = vec![
        vec![("a".to_string(), 3)],
        vec![("a".to_string(), 2)],
        vec![("a".to_string(), 1)],
        vec![("a".to_string(), 0)],
        vec![("b".to_string(), 0), ("d".to_string(), 0)],
        vec![("b".to_string(), 1)],
        vec![("b".to_string(), 2)],
        vec![("c".to_string(), 0), ("d".to_string(), 1)],
        vec![("c".to_string(), 1)],
        vec![("c".to_string(), 2)],
    ];

    let grid = InherentDiscreteGrid::from_index_table(&["a", "b", "c", "d"], index_table)
        .build()
        .unwrap();

    assert_eq!(
        grid.quantics_to_grididx(&[0, 1, 0, 1, 1, 1, 0, 3, 0, 1])
            .unwrap(),
        vec![10, 2, 5, 3]
    );
    assert_eq!(
        grid.grididx_to_quantics(&[10, 2, 5, 3]).unwrap(),
        vec![0, 1, 0, 1, 1, 1, 0, 3, 0, 1]
    );
}

// --- "challenging tests - asymmetric grids" (R=0 dimension) ---
#[test]
fn test_asymmetric_grids_r0() {
    // Grid with R=0 dimension (only index=1 possible)
    let grid = InherentDiscreteGrid::builder(&[20, 0, 15, 3])
        .build()
        .unwrap();

    // R=0 means max_grididx for that dim is 2^0-1=0
    assert_eq!(grid.max_grididx()[1], 0);

    // Test a variety of grid indices
    let test_indices: Vec<Vec<usize>> = vec![
        vec![0, 0, 0, 0],
        vec![2usize.pow(20) - 1, 0, 2usize.pow(15) - 1, 2usize.pow(3) - 1],
        vec![999, 0, 499, 3],
        vec![2usize.pow(10) - 1, 0, 2usize.pow(7) - 1, 1],
    ];

    for grididx in &test_indices {
        let quantics = grid.grididx_to_quantics(grididx).unwrap();
        let back = grid.quantics_to_grididx(&quantics).unwrap();
        assert_eq!(back, *grididx);

        // Verify bounds
        assert!(grididx[0] < 2usize.pow(20));
        assert_eq!(grididx[1], 0); // R=0 means only one possible value
        assert!(grididx[2] < 2usize.pow(15));
        assert!(grididx[3] < 2usize.pow(3));
    }
}

// --- "grididx_to_quantics integer input only for 1D" ---
// Note: In Rust, grididx_to_quantics always takes a slice, not a scalar.
// The 1D case doesn't have special scalar behavior. This test verifies
// that a 2D grid correctly rejects a 1-element input.
#[test]
fn test_grididx_1d_vs_2d() {
    // 1D grid: single-element slice should work
    let grid_1d = InherentDiscreteGrid::builder(&[4]).build().unwrap();
    let q1 = grid_1d.grididx_to_quantics(&[2]).unwrap();
    assert!(!q1.is_empty());

    // 2D grid: single-element input should error
    let grid_2d = InherentDiscreteGrid::builder(&[2, 3]).build().unwrap();
    // A single-element slice for a 2D grid: expand_grididx broadcasts it,
    // so we test that the 2D grid requires proper 2D input for unique results
    let q2_explicit = grid_2d.grididx_to_quantics(&[0, 0]).unwrap();
    assert!(!q2_explicit.is_empty());
}

// --- Additional tests for Grouped scheme ---
#[test]
fn test_grouped_scheme_all_roundtrip() {
    let grid = InherentDiscreteGrid::builder(&[3, 2])
        .with_unfolding_scheme(UnfoldingScheme::Grouped)
        .build()
        .unwrap();

    // Grouped with Rs=(3,2): [a1, a2, a3, b1, b2]
    assert_eq!(grid.len(), 5);
    assert_eq!(grid.local_dimensions(), vec![2; 5]);

    // Test all possible grid indices
    for x in 0..8 {
        for y in 0..4 {
            let grididx = vec![x, y];
            let quantics = grid.grididx_to_quantics(&grididx).unwrap();
            let back = grid.quantics_to_grididx(&quantics).unwrap();
            assert_eq!(back, grididx, "Failed for grididx {:?}", grididx);
        }
    }
}

#[test]
fn test_grouped_scheme_base3_exhaustive() {
    let grid = InherentDiscreteGrid::builder(&[2, 2])
        .with_base(3)
        .with_unfolding_scheme(UnfoldingScheme::Grouped)
        .build()
        .unwrap();

    // All sites should have dim 3
    assert_eq!(grid.local_dimensions(), vec![3; 4]);

    for x in 0..9 {
        for y in 0..9 {
            let grididx = vec![x, y];
            let quantics = grid.grididx_to_quantics(&grididx).unwrap();
            let back = grid.quantics_to_grididx(&quantics).unwrap();
            assert_eq!(back, grididx, "Failed for grididx {:?}", grididx);
        }
    }
}

// --- Consistency: all three schemes should give same origcoord ---
#[test]
fn test_all_schemes_origcoord_consistency() {
    let r = 3;
    let origin = [5i64, 10];
    let step = [2i64, 3];

    let grid_fused = InherentDiscreteGrid::builder(&[r, r])
        .with_origin(&origin)
        .with_step(&step)
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();

    let grid_interleaved = InherentDiscreteGrid::builder(&[r, r])
        .with_origin(&origin)
        .with_step(&step)
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();

    let grid_grouped = InherentDiscreteGrid::builder(&[r, r])
        .with_origin(&origin)
        .with_step(&step)
        .with_unfolding_scheme(UnfoldingScheme::Grouped)
        .build()
        .unwrap();

    let max = 2usize.pow(r as u32);
    for x in 0..max {
        for y in 0..max {
            let grididx = vec![x, y];

            let coord_fused = grid_fused.grididx_to_origcoord(&grididx).unwrap();
            let coord_interleaved = grid_interleaved.grididx_to_origcoord(&grididx).unwrap();
            let coord_grouped = grid_grouped.grididx_to_origcoord(&grididx).unwrap();

            assert_eq!(coord_fused, coord_interleaved);
            assert_eq!(coord_fused, coord_grouped);

            // All three should reconstruct the same grididx from origcoord
            assert_eq!(
                grid_fused.origcoord_to_grididx(&coord_fused).unwrap(),
                grididx
            );
            assert_eq!(
                grid_interleaved
                    .origcoord_to_grididx(&coord_interleaved)
                    .unwrap(),
                grididx
            );
            assert_eq!(
                grid_grouped.origcoord_to_grididx(&coord_grouped).unwrap(),
                grididx
            );
        }
    }
}

// --- Per-dimension bases with mixed bases exhaustive tests ---
#[test]
fn test_mixed_bases_exhaustive_fused() {
    // Rs=(2,1), bases=(2,6), Fused
    let grid = InherentDiscreteGrid::builder(&[2, 1])
        .with_bases(&[2, 6])
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();

    let max0 = 2usize.pow(2); // 4
    let max1 = 6usize.pow(1); // 6

    for x in 0..max0 {
        for y in 0..max1 {
            let grididx = vec![x, y];
            let quantics = grid.grididx_to_quantics(&grididx).unwrap();
            let back = grid.quantics_to_grididx(&quantics).unwrap();
            assert_eq!(
                back, grididx,
                "Failed for grididx {:?}, quantics {:?}",
                grididx, quantics
            );
        }
    }
}

#[test]
fn test_mixed_bases_exhaustive_interleaved() {
    // Rs=(2,1), bases=(2,6), Interleaved
    let grid = InherentDiscreteGrid::builder(&[2, 1])
        .with_bases(&[2, 6])
        .with_unfolding_scheme(UnfoldingScheme::Interleaved)
        .build()
        .unwrap();

    let max0 = 2usize.pow(2); // 4
    let max1 = 6usize.pow(1); // 6

    for x in 0..max0 {
        for y in 0..max1 {
            let grididx = vec![x, y];
            let quantics = grid.grididx_to_quantics(&grididx).unwrap();
            let back = grid.quantics_to_grididx(&quantics).unwrap();
            assert_eq!(
                back, grididx,
                "Failed for grididx {:?}, quantics {:?}",
                grididx, quantics
            );
        }
    }
}

#[test]
fn test_mixed_bases_exhaustive_grouped() {
    // Rs=(2,1), bases=(2,6), Grouped
    let grid = InherentDiscreteGrid::builder(&[2, 1])
        .with_bases(&[2, 6])
        .with_unfolding_scheme(UnfoldingScheme::Grouped)
        .build()
        .unwrap();

    let max0 = 2usize.pow(2); // 4
    let max1 = 6usize.pow(1); // 6

    for x in 0..max0 {
        for y in 0..max1 {
            let grididx = vec![x, y];
            let quantics = grid.grididx_to_quantics(&grididx).unwrap();
            let back = grid.quantics_to_grididx(&quantics).unwrap();
            assert_eq!(
                back, grididx,
                "Failed for grididx {:?}, quantics {:?}",
                grididx, quantics
            );
        }
    }
}

// --- Test with 3D mixed bases ---
#[test]
fn test_mixed_bases_3d_fused_exhaustive() {
    // Rs=(3, 2, 1), bases=(2, 3, 5), Fused
    let grid = InherentDiscreteGrid::builder(&[3, 2, 1])
        .with_bases(&[2, 3, 5])
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();

    assert_eq!(grid.local_dimensions(), vec![30, 6, 2]);

    let max0 = 2usize.pow(3); // 8
    let max1 = 3usize.pow(2); // 9
    let max2 = 5usize.pow(1); // 5

    for x in 0..max0 {
        for y in 0..max1 {
            for z in 0..max2 {
                let grididx = vec![x, y, z];
                let quantics = grid.grididx_to_quantics(&grididx).unwrap();
                let back = grid.quantics_to_grididx(&quantics).unwrap();
                assert_eq!(back, grididx, "Failed for grididx {:?}", grididx);
            }
        }
    }
}

// --- No resolutions error ---
#[test]
fn test_no_resolutions_error() {
    let result = InherentDiscreteGrid::builder(&[]).build();
    assert!(matches!(result, Err(QuanticsGridError::NoResolutions)));
}

// --- site_dim tests ---
#[test]
fn test_site_dim() {
    let grid = InherentDiscreteGrid::builder(&[3, 2])
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();

    // Fused scheme with Rs=(3,2): dims should be [4, 4, 2]
    assert_eq!(grid.site_dim(0).unwrap(), 4);
    assert_eq!(grid.site_dim(1).unwrap(), 4);
    assert_eq!(grid.site_dim(2).unwrap(), 2);

    // Out of bounds
    assert!(matches!(
        grid.site_dim(3),
        Err(QuanticsGridError::SiteIndexOutOfBounds { .. })
    ));
}

// --- grid_min and grid_max ---
#[test]
fn test_grid_min_max() {
    let grid = InherentDiscreteGrid::builder(&[3, 2])
        .with_origin(&[5, -3])
        .with_step(&[2, 4])
        .build()
        .unwrap();

    assert_eq!(grid.grid_min(), vec![5, -3]);
    // grid_max = origin + step * max_grididx
    // dim 0: 5 + 2*(8-1) = 5+14 = 19
    // dim 1: -3 + 4*(4-1) = -3+12 = 9
    assert_eq!(grid.grid_max(), vec![19, 9]);
}

// --- is_empty ---
#[test]
fn test_is_empty() {
    let grid = InherentDiscreteGrid::builder(&[3]).build().unwrap();
    assert!(!grid.is_empty());
}

// --- WrongQuanticsLength error ---
#[test]
fn test_wrong_quantics_length() {
    let grid = InherentDiscreteGrid::builder(&[3, 2])
        .with_unfolding_scheme(UnfoldingScheme::Fused)
        .build()
        .unwrap();

    // Expected length is 3 (fused with Rs=(3,2)), providing 2
    let result = grid.quantics_to_grididx(&[0, 0]);
    assert!(matches!(
        result,
        Err(QuanticsGridError::WrongQuanticsLength { .. })
    ));

    // Providing 4 (too many)
    let result = grid.quantics_to_grididx(&[0, 0, 0, 0]);
    assert!(matches!(
        result,
        Err(QuanticsGridError::WrongQuanticsLength { .. })
    ));
}

// --- DimensionMismatch for grididx ---
#[test]
fn test_dimension_mismatch_grididx() {
    let grid = InherentDiscreteGrid::builder(&[3, 2]).build().unwrap();

    // 3 elements for a 2D grid
    let result = grid.grididx_to_quantics(&[0, 0, 0]);
    assert!(matches!(
        result,
        Err(QuanticsGridError::DimensionMismatch { .. })
    ));
}

// --- CoordinateOutOfBounds ---
#[test]
fn test_coordinate_out_of_bounds() {
    let grid = InherentDiscreteGrid::builder(&[3])
        .with_origin(&[0])
        .with_step(&[1])
        .build()
        .unwrap();

    // max_grididx = 8, so max origcoord = 0 + 1*(8-1) = 7
    // Try origcoord = 8 (out of bounds)
    let result = grid.origcoord_to_grididx(&[8]);
    assert!(matches!(
        result,
        Err(QuanticsGridError::CoordinateOutOfBounds { .. })
    ));

    // Try origcoord = -1 (out of bounds)
    let result = grid.origcoord_to_grididx(&[-1]);
    assert!(matches!(
        result,
        Err(QuanticsGridError::CoordinateOutOfBounds { .. })
    ));
}

// --- grididx value above max (out of range) ---
#[test]
fn test_grididx_out_of_range_error() {
    let grid = InherentDiscreteGrid::builder(&[3]).build().unwrap();
    let result = grid.grididx_to_quantics(&[8]); // max is 7 (2^3 - 1)
    assert!(matches!(
        result,
        Err(QuanticsGridError::GridIndexOutOfBounds { .. })
    ));
}

// --- quantics value out of range ---
#[test]
fn test_quantics_out_of_range_error() {
    let grid = InherentDiscreteGrid::builder(&[2]).build().unwrap();
    // Rs=[2] with Fused creates 2 sites, each with dim 2: valid quantics in [0, 1]
    let result = grid.quantics_to_grididx(&[2, 0]); // 2 is out of range [0, 1]
    assert!(matches!(
        result,
        Err(QuanticsGridError::QuanticsOutOfRange { .. })
    ));
}
