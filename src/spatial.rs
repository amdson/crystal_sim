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

pub struct ParticleGrid {
    pub cell_size: f64,
    pub cells: HashMap<(i64, i64), Vec<Particle>>,
    pub cell_rates: HashMap<(i64, i64), Vec<f64>>,
    pub cell_aggr_rates: HashMap<(i64, i64), f64>,
    // Per-particle physics state for relaxation (parallel to `cells`)
    pub cell_forces:     HashMap<(i64, i64), Vec<DVec2>>,
    pub cell_velocities: HashMap<(i64, i64), Vec<DVec2>>,
    pub cell_new_pos:    HashMap<(i64, i64), Vec<DVec2>>,
    pub cell_is_active:  HashMap<(i64, i64), Vec<bool>>,
    pub cell_is_frozen:  HashMap<(i64, i64), Vec<bool>>,
}

impl ParticleGrid {
    pub fn new(cell_size: f64) -> Self {
        Self {
            cell_size,
            cells: HashMap::new(),
            cell_rates: HashMap::new(),
            cell_aggr_rates: HashMap::new(),
            cell_forces: HashMap::new(),
            cell_velocities: HashMap::new(),
            cell_new_pos: HashMap::new(),
            cell_is_active: HashMap::new(),
            cell_is_frozen: HashMap::new(),
        }
    }

    pub fn cell_key(&self, x: f64, y: f64) -> (i64, i64) {
        (
            x.div_euclid(self.cell_size).floor() as i64,
            y.div_euclid(self.cell_size).floor() as i64,
        )
    }

    pub fn insert(&mut self, p: Particle, rate: f64) {
        let key = self.cell_key(p.pos.x, p.pos.y);
        let pos = p.pos;
        self.cells.entry(key).or_default().push(p);
        self.cell_rates.entry(key).or_default().push(rate);
        self.cell_aggr_rates.entry(key).and_modify(|v| *v += rate).or_insert(rate);
        self.cell_forces.entry(key).or_default().push(DVec2::ZERO);
        self.cell_velocities.entry(key).or_default().push(DVec2::ZERO);
        self.cell_new_pos.entry(key).or_default().push(pos);
        self.cell_is_active.entry(key).or_default().push(false);
        self.cell_is_frozen.entry(key).or_default().push(false);
    }

    pub fn swap_remove(&mut self, cell: (i64, i64), ind: usize) -> (Particle, f64) {
        let removed_p = self.cells.get_mut(&cell).unwrap().swap_remove(ind);
        let particle_rate = self.cell_rates.get_mut(&cell).unwrap().swap_remove(ind);
        self.cell_aggr_rates.entry(cell).and_modify(|v| *v -= particle_rate).or_default();
        self.cell_forces.get_mut(&cell).unwrap().swap_remove(ind);
        self.cell_velocities.get_mut(&cell).unwrap().swap_remove(ind);
        self.cell_new_pos.get_mut(&cell).unwrap().swap_remove(ind);
        self.cell_is_active.get_mut(&cell).unwrap().swap_remove(ind);
        self.cell_is_frozen.get_mut(&cell).unwrap().swap_remove(ind);
        (removed_p, particle_rate)
    }

    pub fn sample_rate(&self, u: f64) -> Option<(((i64, i64), usize), f64)> {
        let total_rate: f64 = self.cell_aggr_rates.values().sum();
        if total_rate == 0.0 {
            return None;
        }
        let mut threshold = u * total_rate;
        for (cell, &cell_rate) in &self.cell_aggr_rates {
            if threshold < cell_rate {
                let rates = self.cell_rates.get(cell).unwrap();
                for (i, &rate) in rates.iter().enumerate() {
                    if threshold < rate {
                        return Some(((*cell, i), rate));
                    }
                    threshold -= rate;
                }
            } else {
                threshold -= cell_rate;
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
                if let Some(v) = self.cells.get(&(gx, gy)) {
                    out.extend(v.iter().cloned());
                }
            }
        }
    }

    /// Total detachment rate across all particles.
    pub fn total_rate(&self) -> f64 {
        self.cell_aggr_rates.values().sum()
    }

    /// Update the stored rate for the particle at `pos`, adjusting the cell
    /// aggregate accordingly.  No-op if no particle with that position exists.
    pub fn set_rate(&mut self, pos: DVec2, new_rate: f64) {
        let cell = self.cell_key(pos.x, pos.y);
        if let Some(particles) = self.cells.get(&cell) {
            if let Some(ind) = particles.iter().position(|p| p.pos == pos) {
                let rates = self.cell_rates.get_mut(&cell).unwrap();
                let old_rate = rates[ind];
                rates[ind] = new_rate;
                if let Some(agg) = self.cell_aggr_rates.get_mut(&cell) {
                    *agg = *agg - old_rate + new_rate;
                }
            }
        }
    }

    /// Iterate over every particle in the grid regardless of position.
    pub fn iter(&self) -> impl Iterator<Item = &Particle> {
        self.cells.values().flat_map(|v| v.iter())
    }

    /// Total number of particles across all cells.
    pub fn len(&self) -> usize {
        self.cells.values().map(|v| v.len()).sum()
    }

    /// Execute a closure on all particles whose cell overlaps the query circle,
    /// avoiding intermediate allocations.
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
                if let Some(v) = self.cells.get(&(gx, gy)) {
                    for p in v {
                        f(p);
                    }
                }
            }
        }
    }

    // ── Physics state helpers ─────────────────────────────────────────────────

    /// Mark a particle active and set its initial velocity.
    pub fn mark_active(&mut self, cell: (i64, i64), ind: usize, initial_vel: DVec2) {
        if let Some(a) = self.cell_is_active.get_mut(&cell) { a[ind] = true; }
        if let Some(v) = self.cell_velocities.get_mut(&cell) { v[ind] = initial_vel; }
    }

    /// Reset physics state to clean defaults for all particles in the given cells.
    /// Forces and velocities → ZERO; new_pos → current particle pos; flags → false.
    pub fn init_physics_for_cells(&mut self, cells_iter: impl IntoIterator<Item = (i64, i64)>) {
        for cell in cells_iter {
            let positions: Vec<DVec2> = match self.cells.get(&cell) {
                Some(v) => v.iter().map(|p| p.pos).collect(),
                None => continue,
            };
            let n = positions.len();

            let fs = self.cell_forces.entry(cell).or_default();
            fs.resize(n, DVec2::ZERO);
            for f in fs.iter_mut() { *f = DVec2::ZERO; }

            let vs = self.cell_velocities.entry(cell).or_default();
            vs.resize(n, DVec2::ZERO);
            for v in vs.iter_mut() { *v = DVec2::ZERO; }

            let nps = self.cell_new_pos.entry(cell).or_default();
            nps.resize(n, DVec2::ZERO);
            for (np, &pos) in nps.iter_mut().zip(positions.iter()) { *np = pos; }

            let acts = self.cell_is_active.entry(cell).or_default();
            acts.resize(n, false);
            for a in acts.iter_mut() { *a = false; }

            let frzs = self.cell_is_frozen.entry(cell).or_default();
            frzs.resize(n, false);
            for f in frzs.iter_mut() { *f = false; }
        }
    }

    /// Apply `cell_new_pos` to particle positions. Particles that cross a cell
    /// boundary are relocated with all parallel state carried over.
    /// Returns the set of cells that received cross-boundary arrivals.
    pub fn commit_positions(&mut self, active_cells: &HashSet<(i64, i64)>) -> HashSet<(i64, i64)> {
        let mut moved_to: HashSet<(i64, i64)> = HashSet::new();

        for &cell in active_cells {
            let n = match self.cells.get(&cell) { Some(v) => v.len(), None => continue };

            // Identify cross-cell movers; update same-cell movers in place.
            let mut cross_inds: Vec<usize> = Vec::new();
            for ind in 0..n {
                let is_active = self.cell_is_active.get(&cell).map_or(false, |a| a[ind]);
                if !is_active { continue; }
                let new_pos = self.cell_new_pos.get(&cell).map_or(DVec2::ZERO, |np| np[ind]);
                if self.cell_key(new_pos.x, new_pos.y) == cell {
                    self.cells.get_mut(&cell).unwrap()[ind].pos = new_pos;
                } else {
                    cross_inds.push(ind);
                }
            }

            // Process cross-cell moves highest-index-first to keep swap_remove valid.
            cross_inds.sort_unstable_by(|a, b| b.cmp(a));
            for ind in cross_inds {
                let new_pos  = self.cell_new_pos[&cell][ind];
                let vel      = self.cell_velocities[&cell][ind];
                let is_act   = self.cell_is_active[&cell][ind];
                let is_frz   = self.cell_is_frozen[&cell][ind];
                let rate     = self.cell_rates[&cell][ind];

                let (mut p, _) = self.swap_remove(cell, ind);
                p.pos = new_pos;
                let target = self.cell_key(new_pos.x, new_pos.y);
                self.insert(p, rate); // pushes with zeroed physics
                let last = self.cells[&target].len() - 1;
                self.cell_velocities.get_mut(&target).unwrap()[last] = vel;
                self.cell_is_active.get_mut(&target).unwrap()[last]  = is_act;
                self.cell_is_frozen.get_mut(&target).unwrap()[last]  = is_frz;
                moved_to.insert(target);
            }
        }
        moved_to
    }

    /// Reset physics state after relaxation completes.
    pub fn clear_physics_for_cells(&mut self, cells_set: &HashSet<(i64, i64)>) {
        for &cell in cells_set {
            let positions: Vec<DVec2> = match self.cells.get(&cell) {
                Some(v) => v.iter().map(|p| p.pos).collect(),
                None => continue,
            };
            if let Some(fs)   = self.cell_forces.get_mut(&cell)     { for f in fs.iter_mut()  { *f = DVec2::ZERO; } }
            if let Some(vs)   = self.cell_velocities.get_mut(&cell) { for v in vs.iter_mut()  { *v = DVec2::ZERO; } }
            if let Some(nps)  = self.cell_new_pos.get_mut(&cell)    { for (np, &pos) in nps.iter_mut().zip(positions.iter()) { *np = pos; } }
            if let Some(acts) = self.cell_is_active.get_mut(&cell)  { for a in acts.iter_mut() { *a = false; } }
            if let Some(frzs) = self.cell_is_frozen.get_mut(&cell)  { for f in frzs.iter_mut() { *f = false; } }
        }
    }
}

// ── Relaxation free functions ─────────────────────────────────────────────────

/// Accumulate LJ forces on every active particle from all neighbours.
/// Also accumulates reaction forces on non-active, non-frozen neighbours
/// for activation checking.  Uses `cells[cell][ind].pos` (committed positions).
/// Split-borrow safe: `cells` and `cell_forces` are distinct fields.
pub fn accumulate_lj_forces(
    cells:          &HashMap<(i64, i64), Vec<Particle>>,
    cell_forces:    &mut HashMap<(i64, i64), Vec<DVec2>>,
    cell_is_active: &HashMap<(i64, i64), Vec<bool>>,
    cell_is_frozen: &HashMap<(i64, i64), Vec<bool>>,
    active_cells:   &HashSet<(i64, i64)>,
    cell_size:      f64,
    max_radius:     f64,
    lj_cutoff_factor: f64,
) {
    let query_r = max_radius * lj_cutoff_factor * 2.0 + 1e-6;
    let cell_r  = (query_r / cell_size).ceil() as i64;

    for &act_cell in active_cells {
        let Some(act_ps)      = cells.get(&act_cell)          else { continue };
        let Some(act_active)  = cell_is_active.get(&act_cell) else { continue };

        for ai in 0..act_ps.len() {
            if !act_active[ai] { continue; }
            let a_pos = act_ps[ai].pos;
            let a_r   = act_ps[ai].radius;

            for dx in -cell_r..=cell_r {
                for dy in -cell_r..=cell_r {
                    let nb_cell = (act_cell.0 + dx, act_cell.1 + dy);
                    let Some(nb_ps) = cells.get(&nb_cell) else { continue };

                    for bi in 0..nb_ps.len() {
                        if nb_cell == act_cell && bi == ai { continue; }
                        let b_pos = nb_ps[bi].pos;
                        let b_r   = nb_ps[bi].radius;

                        let r_vec     = b_pos - a_pos;
                        let r_contact = a_r + b_r;
                        let r_cut     = r_contact * lj_cutoff_factor;
                        let f = -lj_force_vec(r_vec, r_contact, 1.0, r_cut);
                        if f == DVec2::ZERO { continue; }

                        // Force on active particle a
                        if let Some(fs) = cell_forces.get_mut(&act_cell) {
                            if ai < fs.len() { fs[ai] += f; }
                        }

                        // Reaction on b if non-active and non-frozen
                        let b_active = cell_is_active.get(&nb_cell).map_or(false, |a| a[bi]);
                        let b_frozen = cell_is_frozen.get(&nb_cell).map_or(false, |a| a[bi]);
                        if !b_active && !b_frozen {
                            if nb_cell == act_cell {
                                if let Some(fs) = cell_forces.get_mut(&act_cell) {
                                    if bi < fs.len() { fs[bi] -= f; }
                                }
                            } else if let Some(fs) = cell_forces.get_mut(&nb_cell) {
                                if bi < fs.len() { fs[bi] -= f; }
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Velocity update (FIRE-style damped) for all active particles; writes target
/// positions to `cell_new_pos`.  Checks non-active particles in `neighbor_cells`
/// for activation via accumulated reaction forces.
/// Returns `(converged, newly_activated_cells)`.
pub fn step_velocities_and_activate(
    cells:           &HashMap<(i64, i64), Vec<Particle>>,
    cell_forces:     &HashMap<(i64, i64), Vec<DVec2>>,
    cell_velocities: &mut HashMap<(i64, i64), Vec<DVec2>>,
    cell_new_pos:    &mut HashMap<(i64, i64), Vec<DVec2>>,
    cell_is_active:  &mut HashMap<(i64, i64), Vec<bool>>,
    cell_is_frozen:  &HashMap<(i64, i64), Vec<bool>>,
    active_cells:    &HashSet<(i64, i64)>,
    neighbor_cells:  &HashSet<(i64, i64)>,
    alpha:           f64,
    damping:         f64,
    static_friction: f64,
    max_active_cells: usize,
) -> (bool, HashSet<(i64, i64)>) {
    let mut converged = true;

    for &cell in active_cells {
        let Some(particles) = cells.get(&cell)           else { continue };
        let Some(forces)    = cell_forces.get(&cell)     else { continue };
        let n = particles.len();

        for ind in 0..n {
            let is_active = cell_is_active.get(&cell).map_or(false, |a| a[ind]);
            if !is_active { continue; }

            let f_raw     = forces.get(ind).copied().unwrap_or(DVec2::ZERO);
            let f_clamped = f_raw.clamp_length_max(3.0);
            let f_mag     = f_clamped.length();
            let f_eff     = if f_mag > static_friction {
                f_clamped * ((f_mag - static_friction) / f_mag)
            } else {
                DVec2::ZERO
            };

            let v_old = cell_velocities.get(&cell).map_or(DVec2::ZERO, |vs| vs[ind]);
            let v_raw = damping * v_old + alpha * f_eff;
            let v_new = if v_raw.dot(f_eff) >= 0.0 { v_raw } else { DVec2::ZERO };

            cell_velocities.get_mut(&cell).unwrap()[ind] = v_new;
            cell_new_pos.get_mut(&cell).unwrap()[ind] = particles[ind].pos + v_new;

            if f_eff.length_squared() >= static_friction * static_friction
                || v_new.length_squared() >= 1e-12
            {
                converged = false;
            }
        }
    }

    // Check neighbor cells for particles that exceeded friction threshold.
    let mut newly_activated: HashSet<(i64, i64)> = HashSet::new();
    for &cell in neighbor_cells {
        if active_cells.contains(&cell) { continue; }
        if active_cells.len() + newly_activated.len() >= max_active_cells { break; }

        let Some(particles) = cells.get(&cell)       else { continue };
        let Some(forces)    = cell_forces.get(&cell) else { continue };

        for ind in 0..particles.len() {
            let is_frozen = cell_is_frozen.get(&cell).map_or(false, |a| a[ind]);
            let is_active = cell_is_active.get(&cell).map_or(false, |a| a[ind]);
            if is_frozen || is_active { continue; }

            let f     = forces.get(ind).copied().unwrap_or(DVec2::ZERO);
            let f_mag = f.length();
            if f_mag > static_friction {
                let vel = alpha * f * ((f_mag - static_friction) / f_mag);
                cell_is_active.entry(cell).or_default()[ind]  = true;
                cell_velocities.entry(cell).or_default()[ind] = vel;
                cell_new_pos.entry(cell).or_default()[ind]    = particles[ind].pos + vel;
                newly_activated.insert(cell);
                converged = false;
            }
        }
    }

    (converged, newly_activated)
}
