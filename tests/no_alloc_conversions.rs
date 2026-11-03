//! `origcoord_to_grididx_into` must not allocate: it is the hot-path form of
//! the coordinate conversion (one call per function evaluation in the TCI
//! drivers), so this test counts the allocations made by the call itself.
//!
//! The counting allocator belongs to this test binary only.

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use quanticsgrids::{DiscretizedGrid, QuanticsGridError};

static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);
static COUNTING: AtomicBool = AtomicBool::new(false);

struct CountingAllocator;

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if COUNTING.load(Ordering::Relaxed) {
            ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        }
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        if COUNTING.load(Ordering::Relaxed) {
            ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        }
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

#[test]
fn origcoord_to_grididx_into_does_not_allocate() {
    let grid = DiscretizedGrid::builder(&[4, 3])
        .with_lower_bound(&[0.0, -1.0])
        .with_upper_bound(&[1.0, 1.0])
        .build()
        .unwrap();
    let coord = [0.6, 0.1];
    let mut out = vec![0usize; grid.ndims()];

    // warm the code paths up before counting
    grid.origcoord_to_grididx_into(&coord, &mut out).unwrap();

    ALLOCATIONS.store(0, Ordering::Relaxed);
    COUNTING.store(true, Ordering::Relaxed);
    for _ in 0..100 {
        grid.origcoord_to_grididx_into(&coord, &mut out).unwrap();
    }
    COUNTING.store(false, Ordering::Relaxed);

    assert_eq!(
        ALLOCATIONS.load(Ordering::Relaxed),
        0,
        "origcoord_to_grididx_into allocated"
    );
    let expected = grid.origcoord_to_grididx(&coord).unwrap();
    assert_eq!(out, expected);
}

#[test]
fn error_paths_do_not_allocate_either() {
    let grid = DiscretizedGrid::builder(&[4, 4]).build().unwrap();
    let mut out = vec![0usize; grid.ndims()];

    ALLOCATIONS.store(0, Ordering::Relaxed);
    COUNTING.store(true, Ordering::Relaxed);
    let result = grid.origcoord_to_grididx_into(&[0.5, 2.0], &mut out);
    COUNTING.store(false, Ordering::Relaxed);

    assert!(matches!(
        result,
        Err(QuanticsGridError::CoordinateOutOfBounds { dim: 1, .. })
    ));
    assert_eq!(ALLOCATIONS.load(Ordering::Relaxed), 0);
}
