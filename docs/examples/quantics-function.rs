use quanticsgrids::{DiscretizedGrid, quantics_function};

fn main() -> quanticsgrids::Result<()> {
    let grid = DiscretizedGrid::builder(&[4, 4]).build()?;

    let f = |coords: &[f64]| coords[0].sin() + coords[1];
    let fq = quantics_function(&grid, f);

    let quantics = grid.origcoord_to_quantics(&[0.0, 0.5])?;
    let coords = grid.quantics_to_origcoord(&quantics)?;
    let value = fq(&quantics)?;

    assert!((value - (coords[0].sin() + coords[1])).abs() < 1e-12);

    Ok(())
}
