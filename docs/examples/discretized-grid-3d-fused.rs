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
