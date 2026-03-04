use std::collections::HashMap;

/// Unbounded grid-based spatial hash. Cell size is set to the largest bonding
/// cutoff so a single-ring of cells always covers the full interaction range.
pub struct SpatialHash {
    cell_size: f64,
    cells: HashMap<(i64, i64), Vec<usize>>,
}

impl SpatialHash {
    pub fn new(cell_size: f64) -> Self {
        Self { cell_size, cells: HashMap::new() }
    }

    fn key(&self, x: f64, y: f64) -> (i64, i64) {
        (
            x.div_euclid(self.cell_size).floor() as i64,
            y.div_euclid(self.cell_size).floor() as i64,
        )
    }

    pub fn insert(&mut self, idx: usize, x: f64, y: f64) {
        self.cells.entry(self.key(x, y)).or_default().push(idx);
    }

    pub fn remove(&mut self, idx: usize, x: f64, y: f64) {
        let k = self.key(x, y);
        if let Some(v) = self.cells.get_mut(&k) {
            if let Some(p) = v.iter().position(|&i| i == idx) {
                v.swap_remove(p);
            }
        }
    }

    /// All particle indices whose cell overlaps a circle of `radius` around (cx, cy).
    /// Caller must do final distance filtering.
    pub fn query(&self, cx: f64, cy: f64, radius: f64) -> Vec<usize> {
        let min_gx = ((cx - radius) / self.cell_size).floor() as i64;
        let max_gx = ((cx + radius) / self.cell_size).floor() as i64;
        let min_gy = ((cy - radius) / self.cell_size).floor() as i64;
        let max_gy = ((cy + radius) / self.cell_size).floor() as i64;

        let mut result = Vec::new();
        for gx in min_gx..=max_gx {
            for gy in min_gy..=max_gy {
                if let Some(v) = self.cells.get(&(gx, gy)) {
                    result.extend_from_slice(v);
                }
            }
        }
        result
    }
}
