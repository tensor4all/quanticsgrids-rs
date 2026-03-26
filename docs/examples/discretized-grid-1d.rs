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
