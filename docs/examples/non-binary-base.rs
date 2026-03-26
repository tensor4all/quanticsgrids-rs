use quanticsgrids::InherentDiscreteGrid;

fn main() -> quanticsgrids::Result<()> {
    let r = 4;
    let base = 10;
    let grid = InherentDiscreteGrid::builder(&[r])
        .with_origin(&[0])
        .with_step(&[1])
        .with_base(base)
        .build()?;

    let quantics = vec![base as i64; r];
    let grididx = vec![(base as i64).pow(r as u32)];
    let origcoord = vec![(base as i64).pow(r as u32) - 1];

    assert_eq!(grid.quantics_to_grididx(&quantics)?, grididx);
    assert_eq!(grid.grididx_to_quantics(&grididx)?, quantics);
    assert_eq!(grid.grididx_to_origcoord(&grididx)?, origcoord);
    assert_eq!(grid.origcoord_to_grididx(&origcoord)?, grididx);
    assert_eq!(grid.origcoord_to_quantics(&origcoord)?, quantics);

    Ok(())
}
