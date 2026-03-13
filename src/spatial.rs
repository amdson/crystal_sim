use std::collections::{HashMap, HashSet};
use glam::DVec2;

use crate::forces::lj_force_vec;
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

    pub fn query(&self, cx: f64, cy: f64, radius: f64) -> Vec<usize> {
        let mut result = Vec::new();
        self.query_into(cx, cy, radius, &mut result);
        result
    }

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

/// Per-cell struct-of-arrays.  All Vecs inside a `Cell` are kept the same
/// length; `swap_remove` on one must be mirrored across all others.
pub struct Cell {
    pub particles:  Vec<Particle>,
    pub rates:      Vec<f64>,
    pub aggr_rate:  f64,
    pub velocities: Vec<DVec2>,
    pub new_pos:    Vec<DVec2>,
    pub is_active:  Vec<bool>,
    pub is_frozen:  Vec<bool>,
}

impl Cell {
    fn new() -> Self {
        Self {
            particles:  Vec::new(),
            rates:      Vec::new(),
            aggr_rate:  0.0,
            velocities: Vec::new(),
            new_pos:    Vec::new(),
            is_active:  Vec::new(),
            is_frozen:  Vec::new(),
        }
    }

    fn push(&mut self, p: Particle, rate: f64) {
        let pos = p.pos;
        self.aggr_rate += rate;
        self.particles.push(p);
        self.rates.push(rate);
        self.velocities.push(DVec2::ZERO);
        self.new_pos.push(pos);
        self.is_active.push(false);
        self.is_frozen.push(false);
    }

    fn swap_remove(&mut self, ind: usize) -> (Particle, f64) {
        let p    = self.particles.swap_remove(ind);
        let rate = self.rates.swap_remove(ind);
        self.aggr_rate -= rate;
        self.velocities.swap_remove(ind);
        self.new_pos.swap_remove(ind);
        self.is_active.swap_remove(ind);
        self.is_frozen.swap_remove(ind);
        (p, rate)
    }
}

/// Particle grid with stable cell indices.
///
/// Cells are stored in an append-only `Vec<Cell>`; `cell_map` maps grid keys to
/// Vec indices that are permanent for the lifetime of the simulation.  This
/// lets the relaxation hot loop use `cells[idx]` (one array offset) instead of
/// a HashMap lookup.
///
/// `cell_forces` is a separate parallel `Vec` so that `accumulate_lj_forces`
/// can hold `&cells` and `&mut cell_forces` simultaneously via field splitting.
pub struct ParticleGrid {
    pub cell_size:   f64,
    pub cells:       Vec<Cell>,           // append-only; indices are stable forever
    pub cell_forces: Vec<Vec<DVec2>>,     // parallel to `cells`
    pub cell_map:    HashMap<(i64, i64), usize>,
}

impl ParticleGrid {
    pub fn new(cell_size: f64) -> Self {
        Self {
            cell_size,
            cells:       Vec::new(),
            cell_forces: Vec::new(),
            cell_map:    HashMap::new(),
        }
    }

    pub fn cell_key(&self, x: f64, y: f64) -> (i64, i64) {
        (
            x.div_euclid(self.cell_size).floor() as i64,
            y.div_euclid(self.cell_size).floor() as i64,
        )
    }

    /// Get or create a cell for `key`, returning its stable Vec index.
    fn get_or_create(&mut self, key: (i64, i64)) -> usize {
        if let Some(&idx) = self.cell_map.get(&key) {
            return idx;
        }
        let idx = self.cells.len();
        self.cells.push(Cell::new());
        self.cell_forces.push(Vec::new());
        self.cell_map.insert(key, idx);
        idx
    }

    pub fn insert(&mut self, p: Particle, rate: f64) {
        let key = self.cell_key(p.pos.x, p.pos.y);
        let idx = self.get_or_create(key);
        self.cells[idx].push(p, rate);
        self.cell_forces[idx].push(DVec2::ZERO);
    }

    /// Remove the particle at `(cell_idx, ind)`.  The `Cell` struct stays; only
    /// its internal Vecs shrink (swap_remove).
    pub fn swap_remove(&mut self, cell_idx: usize, ind: usize) -> (Particle, f64) {
        let (p, rate) = self.cells[cell_idx].swap_remove(ind);
        self.cell_forces[cell_idx].swap_remove(ind);
        (p, rate)
    }

    /// Sample a particle proportional to its detachment rate.
    /// Returns `((cell_idx, particle_ind), rate)`.
    pub fn sample_rate(&self, u: f64) -> Option<((usize, usize), f64)> {
        let total: f64 = self.cells.iter().map(|c| c.aggr_rate).sum();
        if total == 0.0 { return None; }
        let mut threshold = u * total;
        for (idx, cell) in self.cells.iter().enumerate() {
            if threshold < cell.aggr_rate {
                for (i, &rate) in cell.rates.iter().enumerate() {
                    if threshold < rate {
                        return Some(((idx, i), rate));
                    }
                    threshold -= rate;
                }
            } else {
                threshold -= cell.aggr_rate;
            }
        }
        None
    }

    pub fn query(&self, cx: f64, cy: f64, radius: f64) -> Vec<Particle> {
        let mut result = Vec::new();
        self.query_into(cx, cy, radius, &mut result);
        result
    }

    pub fn query_into(&self, cx: f64, cy: f64, radius: f64, out: &mut Vec<Particle>) {
        out.clear();
        let min_gx = ((cx - radius) / self.cell_size).floor() as i64;
        let max_gx = ((cx + radius) / self.cell_size).floor() as i64;
        let min_gy = ((cy - radius) / self.cell_size).floor() as i64;
        let max_gy = ((cy + radius) / self.cell_size).floor() as i64;

        for gx in min_gx..=max_gx {
            for gy in min_gy..=max_gy {
                if let Some(&idx) = self.cell_map.get(&(gx, gy)) {
                    out.extend_from_slice(&self.cells[idx].particles);
                }
            }
        }
    }

    pub fn total_rate(&self) -> f64 {
        self.cells.iter().map(|c| c.aggr_rate).sum()
    }

    pub fn set_rate(&mut self, pos: DVec2, new_rate: f64) {
        let key = self.cell_key(pos.x, pos.y);
        if let Some(&idx) = self.cell_map.get(&key) {
            let cell = &mut self.cells[idx];
            if let Some(ind) = cell.particles.iter().position(|p| p.pos == pos) {
                cell.aggr_rate += new_rate - cell.rates[ind];
                cell.rates[ind] = new_rate;
            }
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = &Particle> {
        self.cells.iter().flat_map(|c| c.particles.iter())
    }

    pub fn len(&self) -> usize {
        self.cells.iter().map(|c| c.particles.len()).sum()
    }

    pub fn query_iter<F>(&self, cx: f64, cy: f64, radius: f64, mut f: F)
    where
        F: FnMut(&Particle),
    {
        let min_gx = ((cx - radius) / self.cell_size).floor() as i64;
        let max_gx = ((cx + radius) / self.cell_size).floor() as i64;
        let min_gy = ((cy - radius) / self.cell_size).floor() as i64;
        let max_gy = ((cy + radius) / self.cell_size).floor() as i64;

        for gx in min_gx..=max_gx {
            for gy in min_gy..=max_gy {
                if let Some(&idx) = self.cell_map.get(&(gx, gy)) {
                    for p in &self.cells[idx].particles { f(p); }
                }
            }
        }
    }

    /// Convenience: borrow a cell by grid key.
    pub fn get_cell(&self, key: (i64, i64)) -> Option<&Cell> {
        self.cell_map.get(&key).map(|&i| &self.cells[i])
    }

    // ── Physics state helpers ─────────────────────────────────────────────────

    pub fn mark_active(&mut self, key: (i64, i64), ind: usize, initial_vel: DVec2) {
        if let Some(&idx) = self.cell_map.get(&key) {
            let cell = &mut self.cells[idx];
            cell.is_active[ind] = true;
            cell.velocities[ind] = initial_vel;
        }
    }

    pub fn init_physics_for_cells(&mut self, cells_iter: impl IntoIterator<Item = (i64, i64)>) {
        for key in cells_iter {
            let idx = self.get_or_create(key);
            let n = {
                let cell = &mut self.cells[idx];
                for i in 0..cell.particles.len() {
                    cell.velocities[i] = DVec2::ZERO;
                    cell.new_pos[i]    = cell.particles[i].pos;
                    cell.is_active[i]  = false;
                    cell.is_frozen[i]  = false;
                }
                cell.particles.len()
            };
            let fs = &mut self.cell_forces[idx];
            fs.resize(n, DVec2::ZERO);
            for f in fs.iter_mut() { *f = DVec2::ZERO; }
        }
    }

    pub fn commit_positions(&mut self, active_cells: &HashSet<(i64, i64)>) -> HashSet<(i64, i64)> {
        let mut moved_to: HashSet<(i64, i64)> = HashSet::new();

        for &key in active_cells {
            let Some(&cell_idx) = self.cell_map.get(&key) else { continue };
            let n = self.cells[cell_idx].particles.len();

            let mut cross_inds: Vec<usize> = Vec::new();
            for ind in 0..n {
                if !self.cells[cell_idx].is_active[ind] { continue; }
                let new_pos = self.cells[cell_idx].new_pos[ind];
                if self.cell_key(new_pos.x, new_pos.y) == key {
                    self.cells[cell_idx].particles[ind].pos = new_pos;
                } else {
                    cross_inds.push(ind);
                }
            }

            cross_inds.sort_unstable_by(|a, b| b.cmp(a));
            for ind in cross_inds {
                let (new_pos, vel, is_act, is_frz, rate) = {
                    let c = &self.cells[cell_idx];
                    (c.new_pos[ind], c.velocities[ind], c.is_active[ind], c.is_frozen[ind], c.rates[ind])
                };
                let (mut p, _) = self.swap_remove(cell_idx, ind);
                p.pos = new_pos;
                let target_key = self.cell_key(new_pos.x, new_pos.y);
                self.insert(p, rate); // get_or_create + push
                let target_idx = self.cell_map[&target_key];
                let last = self.cells[target_idx].particles.len() - 1;
                let tc = &mut self.cells[target_idx];
                tc.velocities[last] = vel;
                tc.is_active[last]  = is_act;
                tc.is_frozen[last]  = is_frz;
                moved_to.insert(target_key);
            }
        }
        moved_to
    }

    pub fn clear_physics_for_cells(&mut self, cells_set: &HashSet<(i64, i64)>) {
        for &key in cells_set {
            let Some(&idx) = self.cell_map.get(&key) else { continue };
            let cell = &mut self.cells[idx];
            for i in 0..cell.particles.len() {
                cell.velocities[i] = DVec2::ZERO;
                cell.new_pos[i]    = cell.particles[i].pos;
                cell.is_active[i]  = false;
                cell.is_frozen[i]  = false;
            }
            for f in self.cell_forces[idx].iter_mut() { *f = DVec2::ZERO; }
        }
    }
}

// ── Relaxation free functions ─────────────────────────────────────────────────

/// Accumulate LJ forces on every active particle from all neighbours.
/// Also accumulates reaction forces on non-active, non-frozen neighbours
/// for activation checking.
///
/// `cells` and `cell_forces` are parallel Vecs from distinct fields of
/// `ParticleGrid`, so they can be split-borrowed simultaneously.
/// Sequential `&mut cell_forces[idx]` accesses are safe Rust — each index
/// operation creates a temporary borrow that is released before the next.
pub fn accumulate_lj_forces(
    cells:        &[Cell],
    cell_forces:  &mut [Vec<DVec2>],
    cell_map:     &HashMap<(i64, i64), usize>,
    active_cells: &HashSet<(i64, i64)>,
    cell_size:    f64,
    max_radius:   f64,
    lj_cutoff_factor: f64,
) {
    let query_r = max_radius * lj_cutoff_factor * 2.0 + 1e-6;
    let cell_r  = (query_r / cell_size).ceil() as i64;

    for &act_key in active_cells {
        let Some(&ai_idx) = cell_map.get(&act_key) else { continue };
        let act_cell = &cells[ai_idx];

        for ai in 0..act_cell.particles.len() {
            if !act_cell.is_active[ai] { continue; }
            let a_pos = act_cell.particles[ai].pos;
            let a_r   = act_cell.particles[ai].radius;

            for dx in -cell_r..=cell_r {
                for dy in -cell_r..=cell_r {
                    let nb_key = (act_key.0 + dx, act_key.1 + dy);
                    let Some(&nb_idx) = cell_map.get(&nb_key) else { continue };
                    let nb_cell = &cells[nb_idx];

                    for bi in 0..nb_cell.particles.len() {
                        if nb_idx == ai_idx && bi == ai { continue; }
                        let b_pos = nb_cell.particles[bi].pos;
                        let b_r   = nb_cell.particles[bi].radius;

                        let r_contact = a_r + b_r;
                        let f = -lj_force_vec(b_pos - a_pos, r_contact, 1.0,
                                              r_contact * lj_cutoff_factor);
                        if f == DVec2::ZERO { continue; }

                        // Sequential borrows — each is a temporary &mut released immediately.
                        if ai < cell_forces[ai_idx].len() { cell_forces[ai_idx][ai] += f; }

                        if !nb_cell.is_active[bi] && !nb_cell.is_frozen[bi]
                            && bi < cell_forces[nb_idx].len()
                        {
                            cell_forces[nb_idx][bi] -= f;
                        }
                    }
                }
            }
        }
    }
}

/// Velocity update (FIRE-style damped) for all active particles; writes target
/// positions into each cell's `new_pos`.  Checks neighbour cells for activation.
/// Returns `(converged, newly_activated_cells)`.
pub fn step_velocities_and_activate(
    cells:          &mut [Cell],
    cell_forces:    &[Vec<DVec2>],
    cell_map:       &HashMap<(i64, i64), usize>,
    active_cells:   &HashSet<(i64, i64)>,
    neighbor_cells: &HashSet<(i64, i64)>,
    alpha:           f64,
    damping:         f64,
    static_friction: f64,
    max_active_cells: usize,
) -> (bool, HashSet<(i64, i64)>) {
    let mut converged = true;

    for &key in active_cells {
        let Some(&idx) = cell_map.get(&key) else { continue };
        let forces = &cell_forces[idx];
        let cell   = &mut cells[idx];

        for ind in 0..cell.particles.len() {
            if !cell.is_active[ind] { continue; }

            let f_raw     = forces.get(ind).copied().unwrap_or(DVec2::ZERO);
            let f_clamped = f_raw.clamp_length_max(1.0);
            let f_eff = f_clamped; 
            // let f_mag     = f_clamped.length();
            // let f_eff     = if f_mag > static_friction {
            //     f_clamped * ((f_mag - static_friction) / f_mag)
            // } else { DVec2::ZERO };

            let v_raw = damping * cell.velocities[ind] + alpha * f_eff;
            let v_new = if v_raw.dot(f_eff) >= 0.0 { v_raw } else { DVec2::ZERO };

            cell.new_pos[ind]    = cell.particles[ind].pos + v_new;
            cell.velocities[ind] = v_new;

            if f_eff.length_squared() >= static_friction * static_friction
                || v_new.length_squared() >= 1e-12
            {
                converged = false;
            }
        }
    }

    let mut newly_activated: HashSet<(i64, i64)> = HashSet::new();
    for &key in neighbor_cells {
        if active_cells.contains(&key) { continue; }
        if active_cells.len() + newly_activated.len() >= max_active_cells { break; }

        let Some(&idx) = cell_map.get(&key) else { continue };
        let forces = &cell_forces[idx];
        let cell   = &mut cells[idx];

        for ind in 0..cell.particles.len() {
            if cell.is_frozen[ind] || cell.is_active[ind] { continue; }
            let f     = forces.get(ind).copied().unwrap_or(DVec2::ZERO);
            let f_mag = f.length();
            if f_mag > static_friction {
                let vel = alpha * f * ((f_mag - static_friction) / f_mag);
                cell.is_active[ind]  = true;
                cell.velocities[ind] = vel;
                cell.new_pos[ind]    = cell.particles[ind].pos + vel;
                newly_activated.insert(key);
                converged = false;
            }
        }
    }

    (converged, newly_activated)
}
