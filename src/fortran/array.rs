//! 2-D arrays with Fortran 1-based indexing.

// ---------------------------------------------------------------------------
// 2-D array helper with Fortran 1-based indexing (row 0 / column 0 unused).
// ---------------------------------------------------------------------------
pub struct Arr2<T: Copy> {
    pub data: Vec<T>,
    pub stride: usize,
}

impl<T: Copy> Arr2<T> {
    /// Plumbing, not part of the method itself. `Arr2` is a plain 2-D table (rows × columns), standing
    /// in for the original Fortran's statically-sized 2-D arrays — e.g. one row per site, one column
    /// per species, holding whether that species was recorded there (`idata`), whether it counts as a
    /// benchmark species there (`ibench`), or its pooled local frequency (`ffij`). `new` allocates a
    /// table of a given size, pre-filled with a starting value; `at` reads one cell.
    /// Allocate rows 0..=nrows and columns 0..=ncols.
    pub fn new(nrows: usize, ncols: usize, init: T) -> Self {
        Arr2 {
            data: vec![init; (nrows + 1) * (ncols + 1)],
            stride: ncols + 1,
        }
    }
    #[inline(always)]
    pub fn at(&self, i: usize, j: usize) -> T {
        self.data[i * self.stride + j]
    }
    /// Plumbing, not part of the method itself. `set` overwrites one cell of an `Arr2` table; `add`
    /// accumulates a value into it. `add` is how `frescalo` builds up each site's pooled local
    /// frequency for each species — every neighbour's weighted contribution is added into the
    /// site×species table one neighbour at a time as the weights file is read, which is the direct
    /// computation of the neighbourhood-frequency formula (a weighted average of presence/absence
    /// across the neighbourhood) that the method is built on (Hill 2012, p.197–198).
    #[inline(always)]
    pub fn set(&mut self, i: usize, j: usize, v: T) {
        self.data[i * self.stride + j] = v;
    }
    #[inline(always)]
    pub fn add(&mut self, i: usize, j: usize, v: T)
    where
        T: std::ops::AddAssign,
    {
        self.data[i * self.stride + j] += v;
    }
}
