use std::collections::HashSet;
use std::f64::consts::TAU;
use glam::DVec2;

use crate::candidates::{CandidateSite, circle_intersections, site_has_overlap};
use crate::config::SimConfig;
// use crate::forces::lj_force_vec;
use crate::particle::Particle;
use crate::rates::RateCatalog;
use crate::rng::Rng;
use crate::forces::lj_force_vec;
use crate::spatial::{accumulate_lj_forces, step_velocities_and_activate, Cell, ParticleGrid, SpatialHash};

/// Read the CPU timestamp counter — ~4 cycles, no syscall.
#[inline(always)]
fn rdtsc() -> u64 {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: RDTSC is universally supported on x86_64.
    unsafe { core::arch::x86_64::_rdtsc() }
    #[cfg(not(target_arch = "x86_64"))]
    { 0 }
}

/// Measure TSC frequency by comparing against a wall-clock sleep.
fn calibrate_tsc_ghz() -> f64 {
    #[cfg(target_arch = "x86_64")] {
        use std::time::Instant;
        let t0 = Instant::now();
        let c0 = rdtsc();
        std::thread::sleep(std::time::Duration::from_millis(20));
        let elapsed_ns = t0.elapsed().as_nanos() as f64;
        let dc = rdtsc().wrapping_sub(c0) as f64;
        dc / elapsed_ns
    }
    #[cfg(not(target_arch = "x86_64"))]
    { 1.0 }
}

#[repr(usize)]
enum TimerEntry {
    Select,
    Attach,
    Detach,
    Rebuild,
    Relax,
    RegenCandidates,
    Count,
}

static ENUM_STR: [&str; TimerEntry::Count as usize] = [
    "Select", "Attach", "Detach", "Rebuild", "Relax", "RegenCandidates",
];

pub struct Simulation {
    /// Particles stored directly in the spatial grid; indexed by slot.
    pub particle_grid: ParticleGrid,

    /// candidates[type_id] = current valid attachment sites for that type
    pub candidates: Vec<CandidateSite>,
    pub candidates_spatial: SpatialHash,
    pub attach_rates: RateCatalog, // indexed by candidate index

    pub static_friction: f64,

    pub config: SimConfig,
    pub time: f64,
    pub rng: Rng,

    /// Accumulated TSC cycles per timer bucket.
    pub timers: Vec<u64>,
    pub usages: Vec<u64>,
    cpu_ghz: f64,
    
    /// Reusable buffer for spatial hash queries (avoids per-query allocation).
    query_buf: Vec<usize>,
    particle_buf: Vec<f32>, // Flat buffer of (x, y, type_id, radius) for WASM export.
}

/*
Simulation steps
- Select candidate (add particle t at x, y or / remove particle i)
- If attach:
    - Add particle
    - Run relaxation (Nx steps of gradient descent in neighborhood)
    - Replace influenced candidates
- If remove:
    - Remove particle
    - Replace influenced candidates
*/

impl Simulation {
    pub fn new(mut config: SimConfig) -> Self {
        config.init_cache();
        let cell_size = config.max_cutoff().max(0.1)*1.1; //Large cell size

        let mut sim = Self {
            particle_grid: ParticleGrid::new(cell_size),
            candidates: Vec::new(),
            candidates_spatial: SpatialHash::new(cell_size),
            attach_rates: RateCatalog::new(),
            static_friction: config.static_friction,
            time: 0.0,
            rng: Rng::new(config.seed),
            config,
            timers: vec![0u64; TimerEntry::Count as usize],
            usages: vec![0u64; TimerEntry::Count as usize],
            cpu_ghz: calibrate_tsc_ghz(),
            query_buf: Vec::new(),
            particle_buf: Vec::new(),
        };

        if sim.config.n_types() >= 1 {
            let r0 = sim.config.radius(0);
            sim.add_particle_raw(Particle::new(0.0, 0.0, 0, r0));
        }

        let regen_r = sim.config.max_cutoff() * 4.0;
        sim.regen_candidates_near(DVec2::ZERO, regen_r);
        sim
    }

    // ── KMC step loop ────────────────────────────────────────────────────────

    pub fn step(&mut self, n: u32) {
        for _ in 0..n {
            let t0_sel = rdtsc();
            let attach_total = self.attach_rates.total();
            let detach_total = self.particle_grid.total_rate();
            let r_total = attach_total + detach_total;
            if r_total <= 0.0 {
                break;
            }

            let u1 = self.rng.next_f64();
            let u2 = self.rng.next_f64();

            self.time += -(u2.ln()) / r_total;

            if u1 < attach_total / r_total {
                let attach_u = u1 * r_total / attach_total;
                let candidate_ind = self.attach_rates.select(attach_u);
                self.timers[TimerEntry::Select as usize] += rdtsc().wrapping_sub(t0_sel);
                self.usages[TimerEntry::Select as usize] += 1;

                let t0 = rdtsc();
                let site = self.candidates[candidate_ind].clone();
                self.attach(site);
                self.timers[TimerEntry::Attach as usize] += rdtsc().wrapping_sub(t0);
                self.usages[TimerEntry::Attach as usize] += 1;
            } else {
                let detach_u = (u1 - attach_total / r_total) * r_total / detach_total;
                self.timers[TimerEntry::Select as usize] += rdtsc().wrapping_sub(t0_sel);
                self.usages[TimerEntry::Select as usize] += 1;

                let t0 = rdtsc();
                self.detach(detach_u);
                self.timers[TimerEntry::Detach as usize] += rdtsc().wrapping_sub(t0);
                self.usages[TimerEntry::Detach as usize] += 1;
            }
        }

        let t0 = rdtsc();
        self.rebuild_particle_buf();
        self.timers[TimerEntry::Rebuild as usize] += rdtsc().wrapping_sub(t0);
        self.usages[TimerEntry::Rebuild as usize] += 1;

        let ghz = self.cpu_ghz;
        for i in 0..(TimerEntry::Count as usize) {
            let cycles = self.timers[i];
            let calls = self.usages[i];
            let label = ENUM_STR[i];
            let total_ms = cycles as f64 / ghz / 1e6;
            if calls > 0 {
                let avg_ms = cycles as f64 / ghz / 1e6 / calls as f64;
                println!("{label:18} calls={calls:6}  avg={avg_ms:8.3}ms  total={total_ms:.3}ms");
            } else {
                println!("{label:18} calls=     0  avg=       -    total=  0.000ms");
            }
        }
    }

    fn del_candidate(&mut self, site_idx: usize) {
        if site_idx >= self.candidates.len() {
            return;
        }
        let site = &self.candidates[site_idx];
        self.candidates_spatial.remove(site_idx, site.pos.x, site.pos.y);
        self.candidates.swap_remove(site_idx);
        self.attach_rates.swap_remove_rate(site_idx);

        if site_idx < self.candidates.len() {
            let moved_site = &self.candidates[site_idx];
            self.candidates_spatial.remove(self.candidates.len(), moved_site.pos.x, moved_site.pos.y);
            self.candidates_spatial.insert(site_idx, moved_site.pos.x, moved_site.pos.y);
        }
    }

    fn del_candidates_batch(&mut self, mut indices: Vec<usize>) {
        indices.sort_unstable_by(|a, b| b.cmp(a));
        indices.dedup();
        for site_idx in indices {
            self.del_candidate(site_idx);
        }
    }

    // ── Attach ───────────────────────────────────────────────────────────────

    fn attach(&mut self, site: CandidateSite) {
        let type_c = site.type_id;
        let rc = self.config.radius(type_c);
        let pos = site.pos;

        let query_r = rc + self.config.max_radius() + 1e-6;
        let delta = self.config.delta;
        let mut has_overlap = false;
        self.particle_grid.query_iter(pos.x, pos.y, query_r, |p| {
            if has_overlap {
                return;
            }
            let overlap_dist = rc + p.radius - delta;
            if overlap_dist > 0.0 && pos.distance_squared(p.pos) < overlap_dist * overlap_dist {
                has_overlap = true;
            }
        });
        if has_overlap {
            return;
        }

        let new_p = Particle::new(pos.x, pos.y, type_c, rc);
        self.add_particle_raw(new_p.clone());

        if self.config.relax_steps > 0 {
            let t0 = rdtsc();
            self.relax_new_particle_fixed_grid(new_p.clone());
            self.timers[TimerEntry::Relax as usize] += rdtsc().wrapping_sub(t0);
            self.usages[TimerEntry::Relax as usize] += 1;
        } else {
            self.update_detach_rate(&new_p);
            let nbs = self.get_neighbors(&new_p);
            for nb in &nbs {
                self.update_detach_rate(nb);
            }
            let regen_r = (self.config.max_cutoff() + rc) * 2.0;
            let t0 = rdtsc();
            self.regen_candidates_near(pos, regen_r);
            self.timers[TimerEntry::RegenCandidates as usize] += rdtsc().wrapping_sub(t0);
            self.usages[TimerEntry::RegenCandidates as usize] += 1;
        }
    }

    // ── Detach ───────────────────────────────────────────────────────────────

    fn detach(&mut self, u: f64) {
        let ((cell, ind), _) = match self.particle_grid.sample_rate(u) {
            Some(x) => x,
            None => return,
        };

        // Snapshot the particle and its neighbours before removal.
        let detach_p = self.particle_grid.cells[cell].particles[ind].clone();
        let nbs = self.get_neighbors(&detach_p);
        let pos = detach_p.pos;

        self.particle_grid.swap_remove(cell, ind);

        // Old neighbours lost a bond → recompute their detach rates.
        for nb in &nbs {
            self.update_detach_rate(nb);
        }

        let regen_r = self.config.max_cutoff() * 2.0;
        let t0 = rdtsc();
        self.regen_candidates_near(pos, regen_r);
        self.timers[TimerEntry::RegenCandidates as usize] += rdtsc().wrapping_sub(t0);
        self.usages[TimerEntry::RegenCandidates as usize] += 1;
    }

    // ── Internal helpers ─────────────────────────────────────────────────────

    /// Push a particle into the data structures without updating rates/candidates.
    fn add_particle_raw(&mut self, p: Particle) {
        self.particle_grid.insert(p, 0.0); // placeholder rate; caller must call update_detach_rate
    }

    /// Compute E_bind for `p` and refresh its detachment rate in the grid.
    pub fn update_detach_rate(&mut self, p: &Particle) {
        let e_bind = self.binding_energy(p);
        let rate = self.config.detach_rate(e_bind);
        self.particle_grid.set_rate(p.pos, rate);
    }

    pub fn get_neighbors(&self, p: &Particle) -> Vec<Particle> {
        let cutoff = p.radius + self.config.max_radius() + self.config.delta;
        let mut neighbors = Vec::new();
        self.particle_grid.query_iter(p.pos.x, p.pos.y, cutoff, |neighbor| {
            if neighbor.pos != p.pos && p.bonds_to(neighbor, self.config.delta) {
                neighbors.push(neighbor.clone());
            }
        });
        neighbors
    }

    fn binding_energy(&self, p: &Particle) -> f64 {
        let cutoff = p.radius + self.config.max_radius() + self.config.delta;
        let config = &self.config;
        let mut energy = 0.0;
        self.particle_grid.query_iter(p.pos.x, p.pos.y, cutoff, |neighbor| {
            if neighbor.pos != p.pos && p.bonds_to(neighbor, config.delta) {
                energy += config.epsilon(p.type_id, neighbor.type_id);
            }
        });
        energy
    }

    fn site_binding_energy(&self, site: &CandidateSite) -> f64 {
        let type_c = site.type_id;
        let rc = self.config.radius(type_c);
        let pos = site.pos;
        let query_r = rc + self.config.max_radius() + self.config.delta;
        let config = &self.config;
        let mut energy = 0.0;
        self.particle_grid.query_iter(pos.x, pos.y, query_r, |neighbor| {
            let contact = neighbor.radius + rc;
            let d = neighbor.pos.distance(pos);
            if d >= contact - config.delta && d <= contact + config.delta {
                energy += config.epsilon(type_c, neighbor.type_id);
            }
        });
        energy
    }

    fn add_candidate(&mut self, site: CandidateSite) {
        let type_c = site.type_id;
        let pos = site.pos;
        let binding_energy = self.site_binding_energy(&site);

        self.candidates.push(site);
        let idx = self.candidates.len() - 1;
        self.candidates_spatial.insert(idx, pos.x, pos.y);
        let rate = self.config.attach_rate(type_c, binding_energy);
        self.attach_rates.add_rate(rate);
    }

    fn regen_candidates_near(&mut self, center: DVec2, radius: f64) {
        self.regen_candidates_batch(&[(center, radius)]);
    }

    fn regen_candidates_batch(&mut self, regions: &[(DVec2, f64)]) {
        if regions.is_empty() {
            return;
        }

        // ── 1. Batch-delete stale candidates across all regions ───────────
        let mut to_delete: Vec<usize> = Vec::new();
        for &(center, radius) in regions {
            let r_sq = radius * radius;
            self.candidates_spatial.query_into(
                center.x, center.y,
                radius + self.config.max_cutoff(),
                &mut self.query_buf,
            );
            for &site_idx in &self.query_buf {
                if site_idx < self.candidates.len()
                    && self.candidates[site_idx].pos.distance_squared(center) < r_sq
                {
                    to_delete.push(site_idx);
                }
            }
        }
        to_delete.sort_unstable();
        to_delete.dedup();
        self.del_candidates_batch(to_delete);

        // ── 2. Build merged focal set (particles near any region) ─────────
        let mut particle_buf: Vec<Particle> = Vec::new();
        let mut focal_particles: Vec<Particle> = Vec::new();
        for &(center, radius) in regions {
            let search_r = radius + self.config.max_cutoff() * 2.0;
            self.particle_grid.query_into(center.x, center.y, search_r, &mut particle_buf);
            focal_particles.extend_from_slice(&particle_buf);
        }
        focal_particles.sort_unstable_by(|a, b| {
            a.pos.x.total_cmp(&b.pos.x).then(a.pos.y.total_cmp(&b.pos.y))
        });
        focal_particles.dedup_by(|a, b| a.pos == b.pos);

        // ── 3. Build merged overlap-check arrays ──────────────────────────
        let mut overlap_particles: Vec<Particle> = Vec::new();
        for &(center, radius) in regions {
            let search_r = radius + self.config.max_cutoff() * 2.0;
            let overlap_r = search_r + 4.0 * self.config.max_radius();
            self.particle_grid.query_into(center.x, center.y, overlap_r, &mut particle_buf);
            overlap_particles.extend_from_slice(&particle_buf);
        }
        overlap_particles.sort_unstable_by(|a, b| {
            a.pos.x.total_cmp(&b.pos.x).then(a.pos.y.total_cmp(&b.pos.y))
        });
        overlap_particles.dedup_by(|a, b| a.pos == b.pos);
        let overlap_pos: Vec<DVec2> = overlap_particles.iter().map(|p| p.pos).collect();
        let overlap_rad: Vec<f64> = overlap_particles.iter().map(|p| p.radius).collect();
        let overlap_idx: Vec<usize> = (0..overlap_particles.len()).collect();

        // ── 4. For each focal particle A and each candidate type C, generate sites.
        let delta = self.config.delta;
        let max_radius = self.config.max_radius();
        for a in focal_particles {
            let a_pos = a.pos;
            let a_r = a.radius;

            let bonded_cutoff = a_r + max_radius + delta;
            let mut bonded_angles: Vec<f64> = Vec::new();
            self.particle_grid.query_iter(a_pos.x, a_pos.y, bonded_cutoff, |other| {
                if other.pos == a_pos { return; }
                let d = a_pos.distance(other.pos);
                if d < a_r + other.radius + delta {
                    bonded_angles.push(f64::atan2(other.pos.y - a_pos.y, other.pos.x - a_pos.x));
                }
            });

            for type_c in 0..self.config.n_types() {
                let rc = self.config.radius(type_c);

                // ── 4a. Intersection sites ────────────────────────────────
                let pair_search = a_r + rc + max_radius + rc + delta;
                self.particle_grid.query_into(a_pos.x, a_pos.y, pair_search, &mut particle_buf);
                for b in &particle_buf {
                    if (b.pos.x, b.pos.y) <= (a_pos.x, a_pos.y) { continue; }
                    let b_pos = b.pos;
                    let b_r = b.radius;

                    let r_ac = a_r + rc;
                    let r_bc = b_r + rc;

                    for maybe_pt in circle_intersections(a_pos, r_ac, b_pos, r_bc) {
                        let Some(site_pos) = maybe_pt else { continue };

                        let da = a_pos.distance(site_pos);
                        let db = b_pos.distance(site_pos);
                        if da < r_ac - delta || da > r_ac + delta { continue; }
                        if db < r_bc - delta || db > r_bc + delta { continue; }

                        if site_has_overlap(site_pos, rc, delta, &overlap_idx, &overlap_pos, &overlap_rad) {
                            continue;
                        }
                        self.add_candidate(CandidateSite { pos: site_pos, type_id: type_c });
                    }
                }

                // ── 4b. Arc sites ─────────────────────────────────────────
                let n_ang = self.config.num_isolated_angles;
                let exclusion = TAU / n_ang as f64;
                let r_contact = a_r + rc;
                let max_arc = self.config.max_arc_sites_per_type;

                let mut arc_candidates: Vec<DVec2> = Vec::new();
                for k in 0..n_ang {
                    let theta = (k as f64) * TAU / (n_ang as f64);
                    let blocked = bonded_angles.iter().any(|&phi| {
                        let diff = (theta - phi + TAU).rem_euclid(TAU);
                        diff < exclusion || diff > TAU - exclusion
                    });
                    if blocked { continue; }
                    let site_pos = a_pos + DVec2::new(theta.cos(), theta.sin()) * r_contact;
                    if site_has_overlap(site_pos, rc, delta, &overlap_idx, &overlap_pos, &overlap_rad) {
                        continue;
                    }
                    arc_candidates.push(site_pos);
                }

                if arc_candidates.len() <= max_arc {
                    for pos in arc_candidates {
                        self.add_candidate(CandidateSite { pos, type_id: type_c });
                    }
                } else {
                    let step = arc_candidates.len() as f64 / max_arc as f64;
                    for i in 0..max_arc {
                        let idx = (i as f64 * step) as usize;
                        self.add_candidate(CandidateSite { pos: arc_candidates[idx], type_id: type_c });
                    }
                }
            }
        }
    }

    /// Rebuild the flat f32 particle buffer for the WASM export.
    fn rebuild_particle_buf(&mut self) {
        self.particle_buf.clear();
        for p in self.particle_grid.iter() {
            self.particle_buf.push(p.pos.x as f32);
            self.particle_buf.push(p.pos.y as f32);
            self.particle_buf.push(p.type_id as f32);
            self.particle_buf.push(p.radius as f32);
        }
    }

    // ── Public getters ────────────────────────────────────────────────────────

    pub fn particle_count(&self) -> u32 {
        self.particle_grid.len() as u32
    }

    pub fn simulation_time(&self) -> f64 {
        self.time
    }

    pub fn particle_buf_ptr(&self) -> *const f32 {
        self.particle_buf.as_ptr()
    }

    // ── Post-attachment relaxation ────────────────────────────────────────────

    // /// Steepest-descent relaxation of the neighbourhood of a newly attached particle.
    // fn relax_new_particle(&mut self, new_p: Particle) {
    //     if self.config.relax_steps == 0 { return; }

    //     let alpha            = self.config.relax_alpha;
    //     let damping          = self.config.relax_damping;
    //     let static_friction  = self.static_friction;
    //     let lj_cutoff_factor = self.config.lj_cutoff_factor;
    //     let max_radius       = self.config.max_radius();
    //     let max_active_cells: usize = 32;

    //     // The new particle is the last entry in its cell (inserted by add_particle_raw).
    //     let new_cell = self.particle_grid.cell_key(new_p.pos.x, new_p.pos.y);
    //     let new_ind  = self.particle_grid.cells[new_cell].particles.len() - 1;

    //     // Compute the cell-neighborhood radius for LJ interactions.
    //     let query_r = max_radius * lj_cutoff_factor * 2.0 + 1e-6;
    //     let cell_r  = (query_r / self.particle_grid.cell_size).ceil() as i64;

    //     // Helper: expand a cell set by `cell_r` in each direction.
    //     let expand_cells = |src: &HashSet<(i64, i64)>| -> HashSet<(i64, i64)> {
    //         let mut out = HashSet::new();
    //         for &(cx, cy) in src {
    //             for dx in -cell_r..=cell_r {
    //                 for dy in -cell_r..=cell_r {
    //                     out.insert((cx + dx, cy + dy));
    //                 }
    //             }
    //         }
    //         out
    //     };

    //     let mut active_cells: HashSet<(i64, i64)> = std::iter::once(new_cell).collect();
    //     let initial_cells = active_cells.clone();
    //     let mut initialized_cells = expand_cells(&active_cells);
    //     initialized_cells.extend(active_cells.iter());

    //     self.particle_grid.init_physics_for_cells(initialized_cells.iter().copied());
    //     self.particle_grid.mark_active(new_cell, new_ind, DVec2::ZERO);

    //     for _ in 0..self.config.relax_steps {
    //         // Zero forces for all initialized cells.
    //         for &cell in &initialized_cells {
    //             if let Some(fs) = self.particle_grid.cell_forces.get_mut(&cell) {
    //                 for f in fs.iter_mut() { *f = DVec2::ZERO; }
    //             }
    //         }

    //         // Accumulate LJ forces (split borrow: &cells / &mut cell_forces).
    //         {
    //             let cells       = &self.particle_grid.cells;
    //             let cell_forces = &mut self.particle_grid.cell_forces;
    //             accumulate_lj_forces(
    //                 cells, cell_forces,
    //                 &active_cells, self.particle_grid.cell_size, max_radius, lj_cutoff_factor,
    //             );
    //         }

    //         // Velocity step + activation check (split borrow: &mut cells / &cell_forces).
    //         let neighbor_only: HashSet<(i64, i64)> =
    //             initialized_cells.difference(&active_cells).copied().collect();
    //         let (converged, new_cells) = {
    //             let cells       = &mut self.particle_grid.cells;
    //             let cell_forces = &self.particle_grid.cell_forces;
    //             step_velocities_and_activate(
    //                 cells, cell_forces,
    //                 &active_cells, &neighbor_only,
    //                 alpha, damping, static_friction, max_active_cells,
    //             )
    //         };

    //         // Expand active + initialized sets if new cells were activated.
    //         if !new_cells.is_empty() {
    //             let new_neighbors: HashSet<(i64, i64)> = expand_cells(&new_cells)
    //                 .into_iter()
    //                 .filter(|c| !initialized_cells.contains(c))
    //                 .collect();
    //             if !new_neighbors.is_empty() {
    //                 self.particle_grid.init_physics_for_cells(new_neighbors.iter().copied());
    //                 initialized_cells.extend(new_neighbors);
    //             }
    //             active_cells.extend(new_cells);
    //         }

    //         // Commit positions; incorporate any cells particles moved into.
    //         let moved_to = self.particle_grid.commit_positions(&active_cells);
    //         for cell in &moved_to {
    //             if !initialized_cells.contains(cell) {
    //                 self.particle_grid.init_physics_for_cells(std::iter::once(*cell));
    //                 initialized_cells.insert(*cell);
    //             }
    //         }
    //         active_cells.extend(moved_to);

    //         if converged { break; }
    //     }

    //     // ── Post-relaxation: update rates for moved particles and their bonds ──
    //     let delta = self.config.delta;
    //     let mut rate_particles: Vec<Particle> = Vec::new();
    //     for &cell in &active_cells {
    //         if let Some(c) = self.particle_grid.cells.get(&cell) {
    //             rate_particles.extend_from_slice(&c.particles);
    //         }
    //     }
    //     // Include bonded neighbours of active particles.
    //     let mut bonded: Vec<Particle> = Vec::new();
    //     for p in &rate_particles {
    //         self.particle_grid.query_iter(p.pos.x, p.pos.y, p.radius + max_radius + delta, |nb| {
    //             if nb.pos != p.pos && p.bonds_to(nb, delta) { bonded.push(nb.clone()); }
    //         });
    //     }
    //     rate_particles.extend(bonded);
    //     rate_particles.sort_unstable_by(|a, b| a.pos.x.total_cmp(&b.pos.x).then(a.pos.y.total_cmp(&b.pos.y)));
    //     rate_particles.dedup_by(|a, b| a.pos == b.pos);
    //     for p in &rate_particles { self.update_detach_rate(p); }

    //     // Regen candidates for current and initial cell positions.
    //     let regen_r = self.config.max_cutoff() * 2.0;
    //     let mut regen_regions: Vec<(DVec2, f64)> = Vec::new();
    //     for &cell in active_cells.iter().chain(initial_cells.iter()) {
    //         if let Some(c) = self.particle_grid.cells.get(&cell) {
    //             for p in &c.particles { regen_regions.push((p.pos, regen_r)); }
    //         }
    //     }
    //     self.regen_candidates_batch(&regen_regions);

    //     self.particle_grid.clear_physics_for_cells(&initialized_cells);
    // }


    /// Fixed-grid relaxation variant.
    ///
    /// Instead of dynamically growing an active set, this activates **all** particles
    /// in a 5×5 cell neighbourhood around the new particle immediately, backed by a
    /// `(5 + 2*cell_r)²` initialised border.  Raw `Cell` pointers are cached before
    /// the loop so the hot inner loops use array indexing instead of HashMap lookups.
    ///
    /// Swap for `relax_new_particle` at the call site in `attach` to compare.
    fn relax_new_particle_fixed_grid(&mut self, new_p: Particle) {
        if self.config.relax_steps == 0 { return; }

        let alpha            = self.config.relax_alpha;
        let damping          = self.config.relax_damping;
        let static_friction  = self.static_friction;
        let lj_cutoff_factor = self.config.lj_cutoff_factor;
        let max_radius       = self.config.max_radius();
        let cell_size        = self.particle_grid.cell_size;

        let center    = self.particle_grid.cell_key(new_p.pos.x, new_p.pos.y);
        let query_r   = max_radius * lj_cutoff_factor * 2.0 + 1e-6;
        let cell_r    = (query_r / cell_size).ceil() as i64;
        let act_half  = 5i64;              // 5×5 active zone
        let init_half = act_half + cell_r; // margin so every active particle sees full LJ range
        let w         = (2 * init_half + 1) as usize;
        let n         = w * w;

        // ── Build flat key array (row-major, origin = top-left corner) ────────
        let keys: Vec<(i64, i64)> = {
            let mut v = vec![(0i64, 0i64); n];
            for gy in -init_half..=init_half {
                for gx in -init_half..=init_half {
                    v[((gy + init_half) as usize) * w + (gx + init_half) as usize] =
                        (center.0 + gx, center.1 + gy);
                }
            }
            v
        };

        let cell_ind: Vec<Option<usize>> = keys.iter().map(|k| {
            self.particle_grid.cell_map.get(k).copied()
        }).collect();

        let slot_is_active = |s: usize| -> bool {
            let gx = (s % w) as i64 - init_half;
            let gy = (s / w) as i64 - init_half;
            gx.abs() <= act_half && gy.abs() <= act_half
        };

        // ── Init physics; mark every particle in the 5×5 zone active ─────────
        let active_slots: Vec<usize> = (0..n).filter(|&s| slot_is_active(s)).collect();
        self.particle_grid.init_physics_for_cells(keys.iter().copied());
        for &s in &active_slots {
            if let Some(cell_idx) = cell_ind[s] {
                let cell = &mut self.particle_grid.cells[cell_idx];
                cell.is_active.iter_mut().for_each(|a| *a = true);
            }
        }

        let active_cells: HashSet<(i64, i64)> =
            active_slots.iter().map(|&s| keys[s]).collect();
        let active_cell_indices: Vec<usize> =
            active_slots.iter().filter_map(|&s| cell_ind[s]).collect();
        let initialized_cell_indices: Vec<usize> =
            cell_ind.iter().filter_map(|&idx| idx).collect();

        // let t0 = rdtsc();
        // ── Relaxation loop ───────────────────────────────────────────────────
        for _ in 0..self.config.relax_steps {
            // Zero forces.
            for &idx in &initialized_cell_indices {
                if let Some(fs) = self.particle_grid.cell_forces.get_mut(idx) {
                    for f in fs.iter_mut() { *f = DVec2::ZERO; }
                }
            }

            // Accumulate LJ forces via cached cell indices. 
            for &ai_slot in &active_slots {
                let a_ind = cell_ind[ai_slot];
                if a_ind.is_none() { continue; }
                let a_ind = a_ind.unwrap();
                let a_cell = &self.particle_grid.cells[a_ind];
                let a_gx = (ai_slot % w) as i64;
                let a_gy = (ai_slot / w) as i64;

                for ai in 0..a_cell.particles.len() {
                    if !a_cell.is_active[ai] { continue; }
                    let a_pos = a_cell.particles[ai].pos;
                    let a_r   = a_cell.particles[ai].radius;

                    for dx in -cell_r..=cell_r {
                        let nb_gx = a_gx + dx;
                        if nb_gx < 0 || nb_gx >= w as i64 { continue; }
                        for dy in -cell_r..=cell_r {
                            let nb_gy = a_gy + dy;
                            if nb_gy < 0 || nb_gy >= w as i64 { continue; }
                            let nb_slot = nb_gy as usize * w + nb_gx as usize;

                            let b_ind = cell_ind[nb_slot];
                            if b_ind.is_none() { continue; }
                            let b_ind = b_ind.unwrap();
                            let b_cell = &self.particle_grid.cells[b_ind];

                            for bi in 0..b_cell.particles.len() {
                                if nb_slot == ai_slot && bi == ai { continue; }
                                let b_is_active = b_cell.is_active[bi];
                                if b_is_active && (nb_slot < ai_slot || (nb_slot == ai_slot && bi < ai)) {
                                    continue;
                                }
                                let b_pos = b_cell.particles[bi].pos;
                                let b_r   = b_cell.particles[bi].radius;
                                let r_contact = a_r + b_r;
                                let epsilon = self.config.epsilon(a_cell.particles[ai].type_id, b_cell.particles[bi].type_id);
                                let f = -lj_force_vec(b_pos - a_pos, r_contact, epsilon,
                                                        r_contact * lj_cutoff_factor);
                                if f == DVec2::ZERO { continue; }

                                let af = &mut self.particle_grid.cell_forces[a_ind];
                                if ai < af.len() { af[ai] += f; }

                                // Equal-and-opposite reaction is needed for active neighbors
                                // (we skip mirrored active-active pairs), and for movable
                                // inactive neighbors used in activation checks.
                                if b_is_active || !b_cell.is_frozen[bi] {
                                    let bf = &mut self.particle_grid.cell_forces[b_ind];
                                    if bi < bf.len() { bf[bi] -= f; }
                                }
                            }
                        }
                    }
                }
            }

            // Velocity step.
            let mut converged = true;
            for &cp in &active_cell_indices {
                let cell = &mut self.particle_grid.cells[cp];
                let fp = &mut self.particle_grid.cell_forces[cp];

                for ind in 0..cell.particles.len() {
                    if !cell.is_active[ind] { continue; }
                    let f_raw     = if fp.is_empty() { DVec2::ZERO } else { fp[ind] };
                    let f_clamped = f_raw.clamp_length_max(3.0);
                    let f_mag     = f_clamped.length();
                    let f_eff     = if f_mag > static_friction {
                        f_clamped * ((f_mag - static_friction) / f_mag)
                    } else { DVec2::ZERO };

                    let v_raw = damping * cell.velocities[ind] + alpha * f_eff;
                    let v_new = if v_raw.dot(f_eff) >= 0.0 { v_raw } else { DVec2::ZERO };

                    cell.new_pos[ind]    = cell.particles[ind].pos + v_new;
                    cell.velocities[ind] = v_new;

                    // if f_eff.length_squared() >= static_friction * static_friction
                    //     || v_new.length_squared() >= 1e-12
                    // {
                    //     converged = false;
                    // }
                    if f_eff.length_squared() >= 0.05 * 0.05
                    {
                        converged = false;
                    }
                }
            }

            // Commit positions; refresh pointer cache if any cross-cell moves occurred.
            // Cell struct addresses in the HashMap are stable after swap_remove/push
            // (those only mutate the Cell's internal Vecs), but a move outside the
            // preloaded grid would insert a new key and potentially rehash.
            let moved = self.particle_grid.commit_positions(&active_cells);
            if converged { break; }
        }

        // self.timers[TimerEntry::Attach as usize] += rdtsc().wrapping_sub(t0);
        // self.usages[TimerEntry::Attach as usize] += 1;
        // ── Post-relaxation: update rates and regen candidates ────────────────
        let delta   = self.config.delta;
        let regen_r = self.config.max_cutoff() * 2.0;
        let mut rate_particles: Vec<Particle> = Vec::new();
        let mut regen_regions:  Vec<(DVec2, f64)> = Vec::new();
        for &s in &active_slots {
            if let Some(cell_idx) = cell_ind[s] {
                 let cell = &self.particle_grid.cells[cell_idx];
                    rate_particles.extend_from_slice(&cell.particles);
                    for p in &cell.particles { regen_regions.push((p.pos, regen_r)); }
            }
        }

        // Include bonded neighbours outside the active zone. 
        let mut bonded: Vec<Particle> = Vec::new();
        for p in &rate_particles {
            self.particle_grid.query_iter(p.pos.x, p.pos.y, p.radius + max_radius + delta, |nb| {
                if nb.pos != p.pos && p.bonds_to(nb, delta) { bonded.push(nb.clone()); }
            });
        }
        rate_particles.extend(bonded);
        rate_particles.sort_unstable_by(|a, b| {
            a.pos.x.total_cmp(&b.pos.x).then(a.pos.y.total_cmp(&b.pos.y))
        });
        rate_particles.dedup_by(|a, b| a.pos == b.pos);
        for p in &rate_particles { self.update_detach_rate(p); }

        self.regen_candidates_batch(&regen_regions);

        let initialized: HashSet<(i64, i64)> = keys.iter().copied().collect();
        self.particle_grid.clear_physics_for_cells(&initialized);
    }

    /// Return a JSON string with per-type metadata (color, radius) for the renderer.
    pub fn type_metadata_json(&self) -> String {
        let entries: Vec<String> = self
            .config
            .particle_types
            .iter()
            .map(|t| format!(r#"{{"color":"{}","radius":{}}}"#, t.color, t.radius))
            .collect();
        format!("[{}]", entries.join(","))
    }
}
