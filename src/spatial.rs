use std::collections::HashMap;
use glam::DVec2;

use crate::particle::Particle;

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

    /// Remove multiple (idx, x, y) entries in one pass, batching by cell key
    /// to avoid repeated HashMap lookups for entries in the same cell.
    pub fn batch_remove(&mut self, entries: &[(usize, f64, f64)]) {
        for &(idx, x, y) in entries {
            let k = self.key(x, y);
            if let Some(v) = self.cells.get_mut(&k) {
                if let Some(p) = v.iter().position(|&i| i == idx) {
                    v.swap_remove(p);
                }
            }
        }
    }

    /// All particle indices whose cell overlaps a circle of `radius` around (cx, cy).
    /// Caller must do final distance filtering.
    pub fn query(&self, cx: f64, cy: f64, radius: f64) -> Vec<usize> {
        let mut result = Vec::new();
        self.query_into(cx, cy, radius, &mut result);
        result
    }

    /// Like `query`, but appends into a caller-provided buffer (cleared first).
    pub fn query_into(&self, cx: f64, cy: f64, radius: f64, out: &mut Vec<usize>) {
        out.clear();
        let min_gx = ((cx - radius) / self.cell_size).floor() as i64;
        let max_gx = ((cx + radius) / self.cell_size).floor() as i64;
        let min_gy = ((cy - radius) / self.cell_size).floor() as i64;
        let max_gy = ((cy + radius) / self.cell_size).floor() as i64;

        for gx in min_gx..=max_gx {
            for gy in min_gy..=max_gy {
                if let Some(v) = self.cells.get(&(gx, gy)) {
                    out.extend_from_slice(v);
                }
            }
        }
    }
}

// ── ParticleGrid ─────────────────────────────────────────────────────────────

/// Unified particle store: particles live directly inside grid cells.
///
/// Each particle is assigned a "slot" index in `[0, len)` that matches
/// `RateCatalog` indexing.  `swap_remove` mirrors `Vec::swap_remove` and
/// `RateCatalog::swap_remove_rate` so the three stay in sync automatically.
///
/// Internally every cell holds a `Vec<(slot, Particle)>` and `slot_loc[slot]`
/// records `(cell_key, within_cell_index)` for O(1) random access and removal.
pub struct ParticleGrid {
    cell_size: f64,
    /// Particles stored directly in grid cells as (slot, Particle) pairs.
    cells: HashMap<(i64, i64), Vec<(usize, Particle)>>,
    /// slot_loc[slot] = (cell_key, within-cell index).
    /// Always `slot_loc.len() == number of live particles`.
    slot_loc: Vec<((i64, i64), usize)>,
}

impl ParticleGrid {
    pub fn new(cell_size: f64) -> Self {
        Self { cell_size, cells: HashMap::new(), slot_loc: Vec::new() }
    }

    pub fn len(&self) -> usize {
        self.slot_loc.len()
    }

    pub fn is_empty(&self) -> bool {
        self.slot_loc.is_empty()
    }

    fn cell_key(&self, x: f64, y: f64) -> (i64, i64) {
        (
            x.div_euclid(self.cell_size).floor() as i64,
            y.div_euclid(self.cell_size).floor() as i64,
        )
    }

    /// Insert a particle, assigning it slot = current `len()`.  Returns the slot.
    pub fn insert(&mut self, p: Particle) -> usize {
        let slot = self.slot_loc.len();
        let key = self.cell_key(p.pos.x, p.pos.y);
        let cell = self.cells.entry(key).or_default();
        let within = cell.len();
        cell.push((slot, p));
        self.slot_loc.push((key, within));
        slot
    }

    /// Swap-remove: remove the particle at `slot`, moving the particle at
    /// `last` into its place.  Mirrors `Vec::swap_remove` and
    /// `RateCatalog::swap_remove_rate` so all three stay index-aligned.
    ///
    /// Returns the removed particle.
    pub fn swap_remove(&mut self, slot: usize) -> Particle {
        let last = self.slot_loc.len() - 1;
        let (s_key, s_within) = self.slot_loc[slot];

        // ── 1. Remove from the cell via within-cell swap_remove ───────────
        let cell = self.cells.get_mut(&s_key).unwrap();
        let (_, removed_p) = cell.swap_remove(s_within);
        // Fix slot_loc for any entry that was moved within the cell
        if s_within < cell.len() {
            let moved_slot = cell[s_within].0;
            self.slot_loc[moved_slot].1 = s_within;
        }

        if slot == last {
            self.slot_loc.pop();
            return removed_p;
        }

        // ── 2. Move the particle at `last` into slot `slot` ──────────────
        let (l_key, l_within) = self.slot_loc[last];
        self.cells.get_mut(&l_key).unwrap()[l_within].0 = slot;
        self.slot_loc[slot] = (l_key, l_within);
        self.slot_loc.pop();

        removed_p
    }

    /// Move particle at `slot` to `new_pos`, updating cell membership if needed.
    /// Replaces the separate `spatial.remove` + update pos + `spatial.insert` pattern.
    pub fn move_to(&mut self, slot: usize, new_pos: DVec2) {
        let (old_key, old_within) = self.slot_loc[slot];
        let new_key = self.cell_key(new_pos.x, new_pos.y);

        if old_key == new_key {
            // Same cell: just update the position in-place
            self.cells.get_mut(&old_key).unwrap()[old_within].1.pos = new_pos;
        } else {
            // Different cell: remove from old, insert into new
            let old_cell = self.cells.get_mut(&old_key).unwrap();
            let (_, mut p) = old_cell.swap_remove(old_within);
            if old_within < old_cell.len() {
                let moved_slot = old_cell[old_within].0;
                self.slot_loc[moved_slot].1 = old_within;
            }
            p.pos = new_pos;
            let new_cell = self.cells.entry(new_key).or_default();
            let new_within = new_cell.len();
            new_cell.push((slot, p));
            self.slot_loc[slot] = (new_key, new_within);
        }
    }

    /// Return all slots whose cell overlaps the query circle.
    /// Caller must distance-filter the results.
    pub fn query(&self, cx: f64, cy: f64, radius: f64) -> Vec<usize> {
        let mut result = Vec::new();
        self.query_into(cx, cy, radius, &mut result);
        result
    }

    /// Like `query` but appends into a caller-supplied buffer (cleared first).
    pub fn query_into(&self, cx: f64, cy: f64, radius: f64, out: &mut Vec<usize>) {
        out.clear();
        let min_gx = ((cx - radius) / self.cell_size).floor() as i64;
        let max_gx = ((cx + radius) / self.cell_size).floor() as i64;
        let min_gy = ((cy - radius) / self.cell_size).floor() as i64;
        let max_gy = ((cy + radius) / self.cell_size).floor() as i64;

        for gx in min_gx..=max_gx {
            for gy in min_gy..=max_gy {
                if let Some(v) = self.cells.get(&(gx, gy)) {
                    out.extend(v.iter().map(|(slot, _)| slot));
                }
            }
        }
    }

    /// Iterate over `(slot, &Particle)` pairs in unspecified order.
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Particle)> {
        self.cells.values().flat_map(|v| v.iter().map(|(s, p)| (*s, p)))
    }
}

impl std::ops::Index<usize> for ParticleGrid {
    type Output = Particle;
    fn index(&self, slot: usize) -> &Particle {
        let (key, within) = self.slot_loc[slot];
        &self.cells[&key][within].1
    }
}

impl std::ops::IndexMut<usize> for ParticleGrid {
    fn index_mut(&mut self, slot: usize) -> &mut Particle {
        let (key, within) = self.slot_loc[slot];
        &mut self.cells.get_mut(&key).unwrap()[within].1
    }
}

// Support `for p in &particle_grid { ... }` yielding `&Particle`
pub struct ParticleIter<'a> {
    cells_iter: std::collections::hash_map::Values<'a, (i64, i64), Vec<(usize, Particle)>>,
    current: Option<std::slice::Iter<'a, (usize, Particle)>>,
}

impl<'a> Iterator for ParticleIter<'a> {
    type Item = &'a Particle;
    fn next(&mut self) -> Option<&'a Particle> {
        loop {
            if let Some(slice) = &mut self.current {
                if let Some((_, p)) = slice.next() {
                    return Some(p);
                }
            }
            self.current = Some(self.cells_iter.next()?.iter());
        }
    }
}

impl<'a> IntoIterator for &'a ParticleGrid {
    type Item = &'a Particle;
    type IntoIter = ParticleIter<'a>;
    fn into_iter(self) -> ParticleIter<'a> {
        ParticleIter { cells_iter: self.cells.values(), current: None }
    }
}
