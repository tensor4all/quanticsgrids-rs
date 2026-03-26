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
