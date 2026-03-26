use quanticsgrids::InherentDiscreteGrid;

fn main() -> quanticsgrids::Result<()> {
    let r = 4;
    let grid = InherentDiscreteGrid::builder(&[r])
        .with_origin(&[0])
        .with_step(&[1])
        .build()?;

    let quantics = vec![1; r];
    let grididx = vec![1];
    let origcoord = vec![0];

    assert_eq!(grid.quantics_to_grididx(&quantics)?, grididx);
    assert_eq!(grid.grididx_to_quantics(&grididx)?, quantics);
    assert_eq!(grid.grididx_to_origcoord(&grididx)?, origcoord);
    assert_eq!(grid.origcoord_to_grididx(&origcoord)?, grididx);
    assert_eq!(grid.origcoord_to_quantics(&origcoord)?, quantics);

    let quantics = vec![2; r];
    let grididx = vec![2_i64.pow(r as u32)];
    let origcoord = vec![2_i64.pow(r as u32) - 1];

    assert_eq!(grid.quantics_to_grididx(&quantics)?, grididx);
    assert_eq!(grid.grididx_to_quantics(&grididx)?, quantics);
    assert_eq!(grid.grididx_to_origcoord(&grididx)?, origcoord);
    assert_eq!(grid.origcoord_to_grididx(&origcoord)?, grididx);
    assert_eq!(grid.origcoord_to_quantics(&origcoord)?, quantics);

    Ok(())
}
