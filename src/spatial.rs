use std::collections::{HashMap, HashSet};
use glam::Vec2;

use crate::forces::lj_force_vec;
use crate::particle::Particle;

// ── Chunk constants ───────────────────────────────────────────────────────────

pub const CHUNK_W:         usize = 16;
pub const CHUNK_H:         usize = 16;
pub const CELLS_PER_CHUNK: usize = CHUNK_W * CHUNK_H; // 256
pub const CELL_CAP:        usize = 8;  // max particles per cell

// ── Chunk indexing helpers ────────────────────────────────────────────────────

#[inline]
fn chunk_key(gx: i64, gy: i64) -> (i64, i64) {
    (gx.div_euclid(CHUNK_W as i64), gy.div_euclid(CHUNK_H as i64))
}

#[inline]
fn local_cell_offset(gx: i64, gy: i64) -> usize {
    let lx = gx.rem_euclid(CHUNK_W as i64) as usize;
    let ly = gy.rem_euclid(CHUNK_H as i64) as usize;
    ly * CHUNK_W + lx
}

// ── ParticleGrid ──────────────────────────────────────────────────────────────

/// Flat struct-of-arrays spatial grid.
///
/// Cells are identified by a stable `cell_idx`.  All per-particle data lives in
/// parallel flat `Vec`s; a particle at slot `s` in cell `cidx` is at flat index
/// `cidx * CELL_CAP + s`.
///
/// Chunks (16×16 cells) are the allocation unit.  When a new chunk is created,
/// all flat Vecs are extended by `CELLS_PER_CHUNK * CELL_CAP` in one shot —
/// no per-cell heap allocation ever occurs after that.
///
/// `cell_forces` and `cell_torques` are separate flat Vecs so that the force
/// accumulation loop can hold a shared borrow of particle data and a mutable
/// borrow of forces simultaneously via field splitting.
pub struct ParticleGrid {
    pub cell_size:    f32,
    pub total_particles: usize,

    /// chunk key → base *cell index* (not byte index) for that chunk's first cell.
    pub chunk_map: HashMap<(i64, i64), usize>,

    // ── per-cell flat arrays (indexed by cell_idx) ────────────────────────
    pub cell_len:         Vec<u8>,   // number of live particles in this cell
    pub cell_aggr_rate:   Vec<f64>,  // sum of detach rates in this cell
    pub cell_aggr_attach: Vec<f64>,  // sum of attach rates in this cell

    // ── per-particle flat arrays (indexed by cell_idx * CELL_CAP + slot) ─
    pub particles:          Vec<Particle>,
    pub rates:              Vec<f64>,
    pub attach_rates:       Vec<f64>,
    pub velocities:         Vec<Vec2>,
    pub new_pos:            Vec<Vec2>,
    pub ang_velocities:     Vec<f32>,
    pub new_orientation:    Vec<Vec2>,
    pub is_active:          Vec<bool>,
    pub is_frozen:          Vec<bool>,
    /// Persistent frozen flag: survives physics resets; set for preset frozen particles.
    pub permanently_frozen: Vec<bool>,
    /// If true, this particle's detach rate is always kept at 0 by set_rate.
    pub no_detach:          Vec<bool>,
    /// Force accumulator for the relaxation pass (per-particle, same flat index).
    pub cell_forces:        Vec<Vec2>,
    /// Torque accumulator for the relaxation pass.
    pub cell_torques:       Vec<f32>,
}

impl ParticleGrid {
    pub fn new(cell_size: f32) -> Self {
        Self {
            cell_size,
            total_particles: 0,
            chunk_map: HashMap::new(),
            cell_len:         Vec::new(),
            cell_aggr_rate:   Vec::new(),
            cell_aggr_attach: Vec::new(),
            particles:          Vec::new(),
            rates:              Vec::new(),
            attach_rates:       Vec::new(),
            velocities:         Vec::new(),
            new_pos:            Vec::new(),
            ang_velocities:     Vec::new(),
            new_orientation:    Vec::new(),
            is_active:          Vec::new(),
            is_frozen:          Vec::new(),
            permanently_frozen: Vec::new(),
            no_detach:          Vec::new(),
            cell_forces:        Vec::new(),
            cell_torques:       Vec::new(),
        }
    }

    /// Convert world position to grid cell coordinates.
    pub fn cell_key(&self, x: f32, y: f32) -> (i64, i64) {
        (
            x.div_euclid(self.cell_size).floor() as i64,
            y.div_euclid(self.cell_size).floor() as i64,
        )
    }

    /// Look up the flat cell index for grid coordinates `(gx, gy)`.
    /// Returns `None` if the chunk hasn't been allocated yet.
    pub fn cell_idx_lookup(&self, gx: i64, gy: i64) -> Option<usize> {
        self.chunk_map
            .get(&chunk_key(gx, gy))
            .map(|&base| base + local_cell_offset(gx, gy))
    }

    /// Like `cell_idx_lookup`, but allocates the chunk if it doesn't exist.
    pub fn cell_idx_or_create(&mut self, gx: i64, gy: i64) -> usize {
        let ck = chunk_key(gx, gy);
        let base = if let Some(&b) = self.chunk_map.get(&ck) {
            b
        } else {
            let b = self.cell_len.len(); // current cell count = new chunk's base cell index
            let np = CELLS_PER_CHUNK * CELL_CAP;
            self.particles.resize_with(self.particles.len() + np, || Particle::new(0.0, 0.0, 0, 0.0, Vec2::X));
            self.rates.resize(self.rates.len() + np, 0.0f64);
            self.attach_rates.resize(self.attach_rates.len() + np, 0.0f64);
            self.velocities.resize(self.velocities.len() + np, Vec2::ZERO);
            self.new_pos.resize(self.new_pos.len() + np, Vec2::ZERO);
            self.ang_velocities.resize(self.ang_velocities.len() + np, 0.0f32);
            self.new_orientation.resize(self.new_orientation.len() + np, Vec2::X);
            self.is_active.resize(self.is_active.len() + np, false);
            self.is_frozen.resize(self.is_frozen.len() + np, false);
            self.permanently_frozen.resize(self.permanently_frozen.len() + np, false);
            self.no_detach.resize(self.no_detach.len() + np, false);
            self.cell_forces.resize(self.cell_forces.len() + np, Vec2::ZERO);
            self.cell_torques.resize(self.cell_torques.len() + np, 0.0f32);
            self.cell_len.resize(self.cell_len.len() + CELLS_PER_CHUNK, 0u8);
            self.cell_aggr_rate.resize(self.cell_aggr_rate.len() + CELLS_PER_CHUNK, 0.0f64);
            self.cell_aggr_attach.resize(self.cell_aggr_attach.len() + CELLS_PER_CHUNK, 0.0f64);
            self.chunk_map.insert(ck, b);
            b
        };
        base + local_cell_offset(gx, gy)
    }

    pub fn insert(&mut self, p: Particle, rate: f64) {
        let (gx, gy) = self.cell_key(p.pos.x, p.pos.y);
        let cidx = self.cell_idx_or_create(gx, gy);
        let len = self.cell_len[cidx] as usize;
        debug_assert!(len < CELL_CAP, "cell overflow at ({},{})", gx, gy);
        let flat = cidx * CELL_CAP + len;
        self.particles[flat]          = p.clone();
        self.rates[flat]              = rate;
        self.attach_rates[flat]       = 0.0;
        self.velocities[flat]         = Vec2::ZERO;
        self.new_pos[flat]            = p.pos;
        self.ang_velocities[flat]     = 0.0;
        self.new_orientation[flat]    = p.orientation;
        self.is_active[flat]          = false;
        self.is_frozen[flat]          = false;
        self.permanently_frozen[flat] = false;
        self.no_detach[flat]          = false;
        self.cell_forces[flat]        = Vec2::ZERO;
        self.cell_torques[flat]       = 0.0;
        self.cell_aggr_rate[cidx]   += rate;
        self.cell_len[cidx]          += 1;
        self.total_particles         += 1;
    }

    /// Remove the particle at `(cidx, slot)` using a swap with the last slot.
    /// Returns the removed `(Particle, rate)`.
    pub fn swap_remove(&mut self, cidx: usize, slot: usize) -> (Particle, f64) {
        let len = self.cell_len[cidx] as usize;
        debug_assert!(slot < len, "swap_remove: slot {} out of range {}", slot, len);
        let last = len - 1;
        let base = cidx * CELL_CAP;
        let flat      = base + slot;
        let flat_last = base + last;

        // Decrement aggregates using the to-be-removed values.
        self.cell_aggr_rate[cidx]   -= self.rates[flat];
        self.cell_aggr_attach[cidx] -= self.attach_rates[flat];

        let p    = self.particles[flat].clone();
        let rate = self.rates[flat];

        // Swap every parallel array so the removed slot holds the former last element.
        macro_rules! swap_field {
            ($field:ident) => { self.$field.swap(flat, flat_last); };
        }
        swap_field!(particles);
        swap_field!(rates);
        swap_field!(attach_rates);
        swap_field!(velocities);
        swap_field!(new_pos);
        swap_field!(ang_velocities);
        swap_field!(new_orientation);
        swap_field!(is_active);
        swap_field!(is_frozen);
        swap_field!(permanently_frozen);
        swap_field!(no_detach);
        swap_field!(cell_forces);
        swap_field!(cell_torques);

        self.cell_len[cidx]  -= 1;
        self.total_particles -= 1;
        (p, rate)
    }

    /// Sample a particle proportional to its detachment rate.
    /// Returns `((cell_idx, slot), rate)`.
    pub fn sample_rate(&self, u: f64) -> Option<((usize, usize), f64)> {
        let total: f64 = self.cell_aggr_rate.iter().sum();
        if total == 0.0 { return None; }
        let mut threshold = u * total;
        for (cidx, &aggr) in self.cell_aggr_rate.iter().enumerate() {
            if aggr == 0.0 { continue; }
            if threshold < aggr {
                let base = cidx * CELL_CAP;
                let len  = self.cell_len[cidx] as usize;
                for s in 0..len {
                    let rate = self.rates[base + s];
                    if threshold < rate {
                        return Some(((cidx, s), rate));
                    }
                    threshold -= rate;
                }
            } else {
                threshold -= aggr;
            }
        }
        None
    }

    /// Sample a particle proportional to its attachment rate sum.
    /// Returns `(cell_idx, slot)`.
    pub fn sample_attach_rate(&self, u: f64) -> Option<(usize, usize)> {
        let total: f64 = self.cell_aggr_attach.iter().sum();
        if total == 0.0 { return None; }
        let mut threshold = u * total;
        for (cidx, &aggr) in self.cell_aggr_attach.iter().enumerate() {
            if aggr == 0.0 { continue; }
            if threshold < aggr {
                let base = cidx * CELL_CAP;
                let len  = self.cell_len[cidx] as usize;
                for s in 0..len {
                    let rate = self.attach_rates[base + s];
                    if threshold < rate {
                        return Some((cidx, s));
                    }
                    threshold -= rate;
                }
            } else {
                threshold -= aggr;
            }
        }
        None
    }

    pub fn total_rate(&self) -> f64 {
        self.cell_aggr_rate.iter().sum()
    }

    pub fn total_attach_rate(&self) -> f64 {
        self.cell_aggr_attach.iter().sum()
    }

    pub fn query(&self, cx: f32, cy: f32, radius: f32) -> Vec<Particle> {
        let mut result = Vec::new();
        self.query_into(cx, cy, radius, &mut result);
        result
    }

    pub fn query_into(&self, cx: f32, cy: f32, radius: f32, out: &mut Vec<Particle>) {
        out.clear();
        let min_gx = ((cx - radius) / self.cell_size).floor() as i64;
        let max_gx = ((cx + radius) / self.cell_size).floor() as i64;
        let min_gy = ((cy - radius) / self.cell_size).floor() as i64;
        let max_gy = ((cy + radius) / self.cell_size).floor() as i64;
        for gx in min_gx..=max_gx {
            for gy in min_gy..=max_gy {
                if let Some(cidx) = self.cell_idx_lookup(gx, gy) {
                    let base = cidx * CELL_CAP;
                    let len  = self.cell_len[cidx] as usize;
                    for s in 0..len { out.push(self.particles[base + s].clone()); }
                }
            }
        }
    }

    pub fn query_iter<F>(&self, cx: f32, cy: f32, radius: f32, mut f: F)
    where F: FnMut(&Particle) {
        let min_gx = ((cx - radius) / self.cell_size).floor() as i64;
        let max_gx = ((cx + radius) / self.cell_size).floor() as i64;
        let min_gy = ((cy - radius) / self.cell_size).floor() as i64;
        let max_gy = ((cy + radius) / self.cell_size).floor() as i64;
        for gx in min_gx..=max_gx {
            for gy in min_gy..=max_gy {
                if let Some(cidx) = self.cell_idx_lookup(gx, gy) {
                    let base = cidx * CELL_CAP;
                    let len  = self.cell_len[cidx] as usize;
                    for s in 0..len { f(&self.particles[base + s]); }
                }
            }
        }
    }

    pub fn set_rate(&mut self, pos: Vec2, new_rate: f64) {
        let (gx, gy) = self.cell_key(pos.x, pos.y);
        if let Some(cidx) = self.cell_idx_lookup(gx, gy) {
            let base = cidx * CELL_CAP;
            let len  = self.cell_len[cidx] as usize;
            if let Some(s) = (0..len).find(|&s| self.particles[base + s].pos == pos) {
                if self.no_detach[base + s] { return; }
                self.cell_aggr_rate[cidx] += new_rate - self.rates[base + s];
                self.rates[base + s] = new_rate;
            }
        }
    }

    /// Mark the particle at `pos` as permanently frozen (never moves during relaxation).
    pub fn mark_permanently_frozen_at(&mut self, pos: Vec2) {
        let (gx, gy) = self.cell_key(pos.x, pos.y);
        if let Some(cidx) = self.cell_idx_lookup(gx, gy) {
            let base = cidx * CELL_CAP;
            let len  = self.cell_len[cidx] as usize;
            if let Some(s) = (0..len).find(|&s| self.particles[base + s].pos == pos) {
                self.permanently_frozen[base + s] = true;
                self.is_frozen[base + s]          = true;
            }
        }
    }

    /// Mark the particle at `pos` as no-detach (detach rate kept at 0).
    pub fn mark_no_detach_at(&mut self, pos: Vec2) {
        let (gx, gy) = self.cell_key(pos.x, pos.y);
        if let Some(cidx) = self.cell_idx_lookup(gx, gy) {
            let base = cidx * CELL_CAP;
            let len  = self.cell_len[cidx] as usize;
            if let Some(s) = (0..len).find(|&s| self.particles[base + s].pos == pos) {
                self.no_detach[base + s] = true;
                let old = self.rates[base + s];
                self.cell_aggr_rate[cidx] -= old;
                self.rates[base + s] = 0.0;
            }
        }
    }

    pub fn mark_active(&mut self, key: (i64, i64), ind: usize, initial_vel: Vec2) {
        if let Some(cidx) = self.cell_idx_lookup(key.0, key.1) {
            let flat = cidx * CELL_CAP + ind;
            self.is_active[flat]  = true;
            self.velocities[flat] = initial_vel;
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = &Particle> {
        let cell_len = &self.cell_len;
        let particles = &self.particles;
        (0..cell_len.len()).flat_map(move |cidx| {
            let base = cidx * CELL_CAP;
            let len  = cell_len[cidx] as usize;
            (0..len).map(move |s| &particles[base + s])
        })
    }

    pub fn len(&self) -> usize {
        self.total_particles
    }

    /// Iterate over all non-empty (cell_idx, cell_len) pairs.
    pub fn nonempty_cell_iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        (0..self.cell_len.len()).filter_map(|cidx| {
            let len = self.cell_len[cidx] as usize;
            if len > 0 { Some((cidx, len)) } else { None }
        })
    }

    // ── Physics state helpers ─────────────────────────────────────────────────

    pub fn init_physics_for_cells(&mut self, cells_iter: impl IntoIterator<Item = (i64, i64)>) {
        for (gx, gy) in cells_iter {
            // Do not create new chunks here; only init cells that already exist.
            let cidx = match self.cell_idx_lookup(gx, gy) {
                Some(c) => c,
                None => continue,
            };
            let base = cidx * CELL_CAP;
            let len  = self.cell_len[cidx] as usize;
            for s in 0..len {
                let flat = base + s;
                self.velocities[flat]      = Vec2::ZERO;
                self.new_pos[flat]         = self.particles[flat].pos;
                self.ang_velocities[flat]  = 0.0;
                self.new_orientation[flat] = self.particles[flat].orientation;
                self.is_active[flat]       = false;
                self.is_frozen[flat]       = self.permanently_frozen[flat];
                self.cell_forces[flat]     = Vec2::ZERO;
                self.cell_torques[flat]    = 0.0;
            }
        }
    }

    pub fn commit_positions(&mut self, active_cells: &HashSet<(i64, i64)>) -> HashSet<(i64, i64)> {
        let mut moved_to: HashSet<(i64, i64)> = HashSet::new();

        for &(gx, gy) in active_cells {
            let Some(cidx) = self.cell_idx_lookup(gx, gy) else { continue };
            let len = self.cell_len[cidx] as usize;

            let mut cross_inds: Vec<usize> = Vec::new();
            for s in 0..len {
                if !self.is_active[cidx * CELL_CAP + s] { continue; }
                let new_pos = self.new_pos[cidx * CELL_CAP + s];
                if self.cell_key(new_pos.x, new_pos.y) == (gx, gy) {
                    let flat = cidx * CELL_CAP + s;
                    self.particles[flat].pos         = new_pos;
                    self.particles[flat].orientation = self.new_orientation[flat];
                } else {
                    cross_inds.push(s);
                }
            }

            // Process in reverse order so swap_remove doesn't invalidate earlier indices.
            cross_inds.sort_unstable_by(|a, b| b.cmp(a));
            for s in cross_inds {
                let flat = cidx * CELL_CAP + s;
                let new_pos    = self.new_pos[flat];
                let new_ori    = self.new_orientation[flat];
                let vel        = self.velocities[flat];
                let omega      = self.ang_velocities[flat];
                let is_act     = self.is_active[flat];
                let is_frz     = self.is_frozen[flat];
                let perm_frz   = self.permanently_frozen[flat];
                let no_det     = self.no_detach[flat];
                let rate       = self.rates[flat];
                let attach_r   = self.attach_rates[flat];

                let (mut p, _) = self.swap_remove(cidx, s);
                p.pos         = new_pos;
                p.orientation = new_ori;
                let target_key = self.cell_key(new_pos.x, new_pos.y);
                self.insert(p, rate);

                let target_cidx = self.cell_idx_lookup(target_key.0, target_key.1)
                    .expect("target cell must exist after insert");
                let last_slot = self.cell_len[target_cidx] as usize - 1;
                let tf = target_cidx * CELL_CAP + last_slot;
                self.velocities[tf]      = vel;
                self.ang_velocities[tf]  = omega;
                self.is_active[tf]       = is_act;
                self.is_frozen[tf]       = is_frz;
                self.permanently_frozen[tf] = perm_frz;
                self.no_detach[tf]       = no_det;
                // Restore attach_rate that was zeroed by insert.
                self.cell_aggr_attach[target_cidx] += attach_r;
                self.attach_rates[tf] = attach_r;

                moved_to.insert(target_key);
            }
        }
        moved_to
    }

    pub fn clear_physics_for_cells(&mut self, cells_set: &HashSet<(i64, i64)>) {
        for &(gx, gy) in cells_set {
            let Some(cidx) = self.cell_idx_lookup(gx, gy) else { continue };
            let base = cidx * CELL_CAP;
            let len  = self.cell_len[cidx] as usize;
            for s in 0..len {
                let flat = base + s;
                self.velocities[flat]      = Vec2::ZERO;
                self.new_pos[flat]         = self.particles[flat].pos;
                self.ang_velocities[flat]  = 0.0;
                self.new_orientation[flat] = self.particles[flat].orientation;
                self.is_active[flat]       = false;
                self.is_frozen[flat]       = self.permanently_frozen[flat];
                self.cell_forces[flat]     = Vec2::ZERO;
                self.cell_torques[flat]    = 0.0;
            }
        }
    }
}

// ── Relaxation free functions ─────────────────────────────────────────────────

/// Accumulate LJ forces on every active particle from all neighbours.
/// Also accumulates reaction forces on non-active, non-frozen neighbours
/// for activation checking.
///
/// Takes flat array slices from `ParticleGrid` fields; the caller is
/// responsible for field-splitting to satisfy the borrow checker.
pub fn accumulate_lj_forces(
    particles:    &[Particle],
    is_active:    &[bool],
    is_frozen:    &[bool],
    cell_len:     &[u8],
    cell_forces:  &mut [Vec2],
    chunk_map:    &HashMap<(i64, i64), usize>,
    cell_size:    f32,
    active_cells: &HashSet<(i64, i64)>,
    max_radius:   f32,
    lj_cutoff_factor: f32,
) {
    let query_r = max_radius * lj_cutoff_factor * 2.0 + 1e-6;
    let cell_r  = (query_r / cell_size).ceil() as i64;

    let lookup = |gx: i64, gy: i64| -> Option<usize> {
        chunk_map
            .get(&chunk_key(gx, gy))
            .map(|&base| base + local_cell_offset(gx, gy))
    };

    for &(ax, ay) in active_cells {
        let Some(ai_cidx) = lookup(ax, ay) else { continue };
        let ai_base = ai_cidx * CELL_CAP;
        let ai_len  = cell_len[ai_cidx] as usize;

        for ai in 0..ai_len {
            if !is_active[ai_base + ai] { continue; }
            let a_pos = particles[ai_base + ai].pos;
            let a_r   = particles[ai_base + ai].radius;

            for dx in -cell_r..=cell_r {
                for dy in -cell_r..=cell_r {
                    let Some(nb_cidx) = lookup(ax + dx, ay + dy) else { continue };
                    let nb_base = nb_cidx * CELL_CAP;
                    let nb_len  = cell_len[nb_cidx] as usize;

                    for bi in 0..nb_len {
                        if nb_cidx == ai_cidx && bi == ai { continue; }
                        let b_pos = particles[nb_base + bi].pos;
                        let b_r   = particles[nb_base + bi].radius;

                        let r_contact = a_r + b_r;
                        let f = -lj_force_vec(b_pos - a_pos, r_contact, 1.0,
                                              r_contact * lj_cutoff_factor);
                        if f == Vec2::ZERO { continue; }

                        cell_forces[ai_base + ai] += f;

                        if !is_active[nb_base + bi] && !is_frozen[nb_base + bi] {
                            cell_forces[nb_base + bi] -= f;
                        }
                    }
                }
            }
        }
    }
}

/// Velocity update (FIRE-style damped) for all active particles.
/// Returns `(converged, newly_activated_cells)`.
pub fn step_velocities_and_activate(
    particles:       &mut [Particle],
    velocities:      &mut [Vec2],
    new_pos:         &mut [Vec2],
    ang_velocities:  &mut [f32],
    new_orientation: &mut [Vec2],
    is_active:       &mut [bool],
    is_frozen:       &[bool],
    cell_len:        &[u8],
    cell_forces:     &[Vec2],
    cell_torques:    &[f32],
    chunk_map:       &HashMap<(i64, i64), usize>,
    active_cells:    &HashSet<(i64, i64)>,
    neighbor_cells:  &HashSet<(i64, i64)>,
    alpha:           f32,
    damping:         f32,
    static_friction: f32,
    max_active_cells: usize,
) -> (bool, HashSet<(i64, i64)>) {
    let mut converged = true;

    let lookup = |gx: i64, gy: i64| -> Option<usize> {
        chunk_map
            .get(&chunk_key(gx, gy))
            .map(|&base| base + local_cell_offset(gx, gy))
    };

    for &(gx, gy) in active_cells {
        let Some(cidx) = lookup(gx, gy) else { continue };
        let base = cidx * CELL_CAP;
        let len  = cell_len[cidx] as usize;

        for ind in 0..len {
            if !is_active[base + ind] { continue; }
            let flat = base + ind;

            let f_raw     = cell_forces[flat];
            let f_clamped = f_raw.clamp_length_max(1.0);
            let f_eff     = f_clamped;

            let v_raw = damping * velocities[flat] + alpha * f_eff;
            let v_new = if v_raw.dot(f_eff) >= 0.0 { v_raw } else { Vec2::ZERO };

            new_pos[flat]    = particles[flat].pos + v_new;
            velocities[flat] = v_new;

            if f_eff.length_squared() >= static_friction * static_friction
                || v_new.length_squared() >= 1e-12
            {
                converged = false;
            }
        }
    }

    let mut newly_activated: HashSet<(i64, i64)> = HashSet::new();
    for &(gx, gy) in neighbor_cells {
        if active_cells.contains(&(gx, gy)) { continue; }
        if active_cells.len() + newly_activated.len() >= max_active_cells { break; }

        let Some(cidx) = lookup(gx, gy) else { continue };
        let base = cidx * CELL_CAP;
        let len  = cell_len[cidx] as usize;

        for ind in 0..len {
            let flat = base + ind;
            if is_frozen[flat] || is_active[flat] { continue; }
            let f     = cell_forces[flat];
            let f_mag = f.length();
            if f_mag > static_friction {
                let vel = alpha * f * ((f_mag - static_friction) / f_mag);
                is_active[flat]  = true;
                velocities[flat] = vel;
                new_pos[flat]    = particles[flat].pos + vel;
                newly_activated.insert((gx, gy));
                converged = false;
            }
        }
    }

    (converged, newly_activated)
}
