use std::collections::HashSet;
use std::f32::consts::{PI, TAU};
use glam::Vec2;

use crate::candidates::{CandidateSite, circle_intersections, site_has_overlap};
use crate::config::SimConfig;
use crate::forces::{patchy_force_torque, patchy_pair_energy};
use crate::particle::Particle;
use crate::rates::RateCatalog;
use crate::rng::Rng;
use crate::spatial::{ParticleGrid, SpatialHash};

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

#[derive(Clone)]
struct TestParticleState {
    particle: Particle,
    velocity: Vec2,
    omega: f32,
    frozen: bool,
}

/// Long-lived scratch buffers for `regen_candidates_batch`.
/// Kept between calls so that `clear()` reuses allocated capacity.
#[derive(Default)]
struct RegenScratch {
    to_delete:        Vec<usize>,
    particle_buf:     Vec<Particle>,
    focal_particles:  Vec<Particle>,
    overlap_particles: Vec<Particle>,
    overlap_pos:      Vec<Vec2>,
    overlap_rad:      Vec<f32>,
    overlap_idx:      Vec<usize>,
    bonded_dirs:      Vec<Vec2>,
    arc_candidates:   Vec<Vec2>,
}

pub struct Simulation {
    /// Particles stored directly in the spatial grid; indexed by slot.
    pub particle_grid: ParticleGrid,

    /// candidates[type_id] = current valid attachment sites for that type
    pub candidates: Vec<CandidateSite>,
    pub candidates_spatial: SpatialHash,
    pub attach_rates: RateCatalog, // indexed by candidate index

    pub static_friction: f32,

    pub config: SimConfig,
    pub time: f64,
    pub rng: Rng,

    /// Accumulated TSC cycles per timer bucket.
    pub timers: Vec<u64>,
    pub usages: Vec<u64>,
    cpu_ghz: f64,

    /// Reusable buffer for spatial hash queries (avoids per-query allocation).
    query_buf: Vec<usize>,
    particle_buf: Vec<f32>, // Flat buffer of (x, y, type_id, radius, orientation) for WASM export.
    candidate_buf: Vec<f32>, // Flat buffer of (x, y, type_id, rate) for debug overlay.
    regen_scratch: RegenScratch,

    /// Deterministic test-mode state (only used when `config.testing_mode` is true).
    test_particles: Vec<TestParticleState>,
    test_step_counter: u64,
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
            candidate_buf: Vec::new(),
            regen_scratch: RegenScratch::default(),
            test_particles: Vec::new(),
            test_step_counter: 0,
        };

        if sim.config.testing_mode {
            sim.load_testing_particles_from_config();
            sim.sync_testing_particles_to_grid();
            sim.rebuild_particle_buf();
            println!(
                "[testing_mode] loaded {} particles, patch_types={}, has_patches={}",
                sim.test_particles.len(),
                sim.config.patch_type_count,
                sim.config.has_patches()
            );
            if sim.test_particles.len() < 2 {
                println!(
                    "[testing_mode] warning: fewer than 2 particles loaded; net interaction torques will be zero."
                );
            }
        } else {
            if sim.config.n_types() >= 1 {
                let r0 = sim.config.radius(0);
                sim.add_particle_raw(Particle::new(0.0, 0.0, 0, r0, Vec2::X));
            }

            let regen_r = sim.config.max_cutoff() * 4.0;
            sim.regen_candidates_near(Vec2::ZERO, regen_r);
        }
        sim
    }

    // ── KMC step loop ────────────────────────────────────────────────────────

    pub fn step(&mut self, n: u32) {
        if self.config.testing_mode {
            self.step_testing(n);
            self.rebuild_particle_buf();
            return;
        }

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
                self.attach(candidate_ind);
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
        self.rebuild_candidate_buf();
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

    pub fn is_testing_mode(&self) -> bool {
        self.config.testing_mode
    }

    fn load_testing_particles_from_config(&mut self) {
        self.test_particles.clear();

        if self.config.initial_particles.is_empty() {
            if self.config.n_types() >= 1 {
                let r0 = self.config.radius(0);
                self.test_particles.push(TestParticleState {
                    particle: Particle::new(0.0, 0.0, 0, r0, Vec2::X),
                    velocity: Vec2::ZERO,
                    omega: 0.0,
                    frozen: false,
                });
            }
            return;
        }

        for init in self.config.initial_particles.clone() {
            if init.type_id >= self.config.n_types() {
                continue;
            }
            let radius = self.config.radius(init.type_id);
            let ori_rad = init.orientation_deg * PI / 180.0;
            let orientation = Vec2::new(ori_rad.cos(), ori_rad.sin());
            self.test_particles.push(TestParticleState {
                particle: Particle::new(init.x, init.y, init.type_id, radius, orientation),
                velocity: Vec2::ZERO,
                omega: 0.0,
                frozen: init.frozen,
            });
        }
    }

    fn sync_testing_particles_to_grid(&mut self) {
        let cell_size = self.particle_grid.cell_size;
        self.particle_grid = ParticleGrid::new(cell_size);
        self.candidates.clear();
        self.candidates_spatial = SpatialHash::new(cell_size);
        self.attach_rates = RateCatalog::new();

        for s in &self.test_particles {
            self.particle_grid.insert(s.particle.clone(), 0.0);
        }
    }

    fn compute_testing_diagnostics(&self) -> (Vec<f32>, Vec<Vec2>, Vec<f32>) {
        let n = self.test_particles.len();
        let mut energies = vec![0.0f32; n];
        let mut forces = vec![Vec2::ZERO; n];
        let mut torques = vec![0.0f32; n];

        for i in 0..n {
            let a = &self.test_particles[i].particle;
            for j in (i + 1)..n {
                let b = &self.test_particles[j].particle;
                let r_vec = b.pos - a.pos;
                let r_contact = a.radius + b.radius;

                let pair_energy = if self.config.has_patches() {
                    patchy_pair_energy(
                        r_vec,
                        r_contact,
                        &self.config.particle_types[a.type_id].patches,
                        a.orientation,
                        &self.config.particle_types[b.type_id].patches,
                        b.orientation,
                        &self.config.patch_lut,
                        self.config.patch_type_count,
                        self.config.epsilon(a.type_id, b.type_id),
                        self.config.delta,
                    )
                } else if a.bonds_to(b, self.config.delta) {
                    self.config.epsilon(a.type_id, b.type_id)
                } else {
                    0.0
                };

                energies[i] += pair_energy;
                energies[j] += pair_energy;

                let (f_i, tau_i, tau_j) = patchy_force_torque(
                    r_vec,
                    r_contact,
                    &self.config.particle_types[a.type_id].patches,
                    a.orientation,
                    &self.config.particle_types[b.type_id].patches,
                    b.orientation,
                    &self.config.patch_lut,
                    self.config.patch_type_count,
                    self.config.epsilon(a.type_id, b.type_id),
                    self.config.lj_cutoff_factor,
                );

                forces[i] += f_i;
                forces[j] -= f_i;
                torques[i] += tau_i;
                torques[j] += tau_j;
            }
        }

        (energies, forces, torques)
    }

    fn step_testing(&mut self, n: u32) {
        let alpha = self.config.relax_alpha;
        let damping = self.config.relax_damping;
        let static_friction = self.static_friction;

        for _ in 0..n {
            let (energies, forces, torques) = self.compute_testing_diagnostics();
            self.test_step_counter += 1;
            self.time += 1.0;

            if self.config.print_test_diagnostics {
                let max_abs_tau = torques
                    .iter()
                    .fold(0.0f32, |m, &t| if t.abs() > m { t.abs() } else { m });
                println!(
                    "=== TEST STEP {} === particles={} max|T|={:.6e}",
                    self.test_step_counter,
                    self.test_particles.len(),
                    max_abs_tau
                );
                for i in 0..self.test_particles.len() {
                    let p = &self.test_particles[i].particle;
                    println!(
                        "p[{i}] type={} pos=({:.4},{:.4}) ori={:.4} E={:.6e} F=({:.6e},{:.6e}) T={:.6e} frozen={}",
                        p.type_id,
                        p.pos.x,
                        p.pos.y,
                        p.orientation.y.atan2(p.orientation.x),
                        energies[i],
                        forces[i].x,
                        forces[i].y,
                        torques[i],
                        self.test_particles[i].frozen
                    );
                }
            }

            for i in 0..self.test_particles.len() {
                if self.test_particles[i].frozen {
                    continue;
                }

                // Linear integration (same style as relaxation update).
                let f_raw = forces[i].clamp_length_max(3.0);
                let f_mag = f_raw.length();
                let f_eff = if f_mag > static_friction {
                    f_raw * ((f_mag - static_friction) / f_mag)
                } else {
                    Vec2::ZERO
                };
                let v_raw = damping * self.test_particles[i].velocity + alpha * f_eff;
                let v_new = if v_raw.dot(f_eff) >= 0.0 { v_raw } else { Vec2::ZERO };
                self.test_particles[i].velocity = v_new;
                self.test_particles[i].particle.pos += v_new;

                // Angular integration (same style as relaxation update).
                let tau_raw = torques[i].clamp(-3.0, 3.0);
                let tau_eff = if tau_raw.abs() > static_friction { tau_raw } else { 0.0 };
                let omega_raw = damping * self.test_particles[i].omega + alpha * tau_eff;
                let omega_new = if omega_raw * tau_eff >= 0.0 { omega_raw } else { 0.0 };
                self.test_particles[i].omega = omega_new;
                let (s, c) = omega_new.sin_cos();
                let ori = self.test_particles[i].particle.orientation;
                self.test_particles[i].particle.orientation =
                    Vec2::new(ori.x * c - ori.y * s, ori.x * s + ori.y * c);
            }

            self.sync_testing_particles_to_grid();
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

    fn attach(&mut self, candidate_ind: usize) {
        let site = self.candidates[candidate_ind].clone();
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
            self.del_candidate(candidate_ind);
            return;
        }

        let new_p = Particle::new(pos.x, pos.y, type_c, rc, site.orientation);
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

    fn binding_energy(&self, p: &Particle) -> f32 {
        let cutoff = p.radius + self.config.max_radius() + self.config.delta;
        let config = &self.config;
        let mut energy = 0.0f32;
        self.particle_grid.query_iter(p.pos.x, p.pos.y, cutoff, |neighbor| {
            if neighbor.pos == p.pos { return; }
            if config.has_patches() {
                let r_vec = neighbor.pos - p.pos;
                let r_contact = p.radius + neighbor.radius;
                // patchy_pair_energy returns eps * bump * g, which is negative for
                // attractive patches (eps < 0). Negate so binding_energy > 0 when bound.
                energy -= patchy_pair_energy(
                    r_vec,
                    r_contact,
                    &config.particle_types[p.type_id].patches,
                    p.orientation,
                    &config.particle_types[neighbor.type_id].patches,
                    neighbor.orientation,
                    &config.patch_lut,
                    config.patch_type_count,
                    config.epsilon(p.type_id, neighbor.type_id),
                    config.delta,
                );
            } else if p.bonds_to(neighbor, config.delta) {
                energy += config.epsilon(p.type_id, neighbor.type_id);
            }
        });
        energy
    }

    fn site_binding_energy(&self, site: &CandidateSite) -> f32 {
        let type_c = site.type_id;
        let rc = self.config.radius(type_c);
        let pos = site.pos;
        let query_r = rc + self.config.max_radius() + self.config.delta;
        let config = &self.config;
        let mut energy = 0.0f32;
        self.particle_grid.query_iter(pos.x, pos.y, query_r, |neighbor| {
            if config.has_patches() {
                let r_vec = neighbor.pos - pos;
                let r_contact = rc + neighbor.radius;
                // Negate: patchy_pair_energy < 0 for attractive patches (eps < 0).
                energy -= patchy_pair_energy(
                    r_vec,
                    r_contact,
                    &config.particle_types[type_c].patches,
                    site.orientation,
                    &config.particle_types[neighbor.type_id].patches,
                    neighbor.orientation,
                    &config.patch_lut,
                    config.patch_type_count,
                    config.epsilon(type_c, neighbor.type_id),
                    config.delta,
                );
            } else {
                let contact = neighbor.radius + rc;
                let d = neighbor.pos.distance(pos);
                if d >= contact - config.delta && d <= contact + config.delta {
                    energy += config.epsilon(type_c, neighbor.type_id);
                }
            }
        });
        energy
    }

    fn best_candidate_orientation(&self, pos: Vec2, type_c: usize) -> Vec2 {
        let config = &self.config;
        if !config.has_patches() || config.particle_types[type_c].patches.is_empty() {
            return Vec2::X;
        }

        let rc = config.radius(type_c);
        let query_r = rc + config.max_radius() + config.delta;
        let nearby = self.particle_grid.query(pos.x, pos.y, query_r);
        if nearby.is_empty() {
            return Vec2::X;
        }

        // Build candidate orientations as unit vectors.
        // For each neighbour, the ideal orientation that aligns patch toward it is:
        //   phi_v rotated by -position_rad  =  phi_v * conj(position_cs)
        let mut orientation_candidates: Vec<Vec2> = vec![Vec2::X];
        for nb in &nearby {
            let diff = nb.pos - pos;
            if diff.length_squared() < 1e-24 { continue; }
            let phi_v = diff.normalize();
            for patch in &config.particle_types[type_c].patches {
                let cs = patch.position_cs;
                // rotate phi_v by -position_rad: phi_v * conj(cs)
                let cand = Vec2::new(
                    phi_v.x * cs.x + phi_v.y * cs.y,
                    phi_v.y * cs.x - phi_v.x * cs.y,
                );
                // Only add if not already close to an existing candidate (dot product ~ 1).
                if !orientation_candidates.iter().any(|&u| u.dot(cand) > 1.0 - 1e-9) {
                    orientation_candidates.push(cand);
                }
            }
        }

        let mut best_orientation = orientation_candidates[0];
        let mut best_energy = f32::INFINITY;
        for &orientation in &orientation_candidates {
            let mut energy = 0.0f32;
            for nb in &nearby {
                let r_vec = nb.pos - pos;
                let r_contact = rc + nb.radius;
                energy += patchy_pair_energy(
                    r_vec,
                    r_contact,
                    &config.particle_types[type_c].patches,
                    orientation,
                    &config.particle_types[nb.type_id].patches,
                    nb.orientation,
                    &config.patch_lut,
                    config.patch_type_count,
                    config.epsilon(type_c, nb.type_id),
                    config.delta,
                );
            }
            if energy < best_energy {
                best_energy = energy;
                best_orientation = orientation;
            }
        }

        best_orientation
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

    fn regen_candidates_near(&mut self, center: Vec2, radius: f32) {
        self.regen_candidates_batch(&[(center, radius)]);
    }

    fn regen_candidates_batch(&mut self, regions: &[(Vec2, f32)]) {
        if regions.is_empty() {
            return;
        }

        let s = &mut self.regen_scratch;

        // ── 1. Batch-delete stale candidates across all regions ───────────
        s.to_delete.clear();
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
                    s.to_delete.push(site_idx);
                }
            }
        }
        s.to_delete.sort_unstable();
        s.to_delete.dedup();
        // del_candidates_batch needs ownership; swap out to avoid holding &mut s and &mut self.
        let to_delete = std::mem::take(&mut self.regen_scratch.to_delete);
        self.del_candidates_batch(to_delete);
        // Reclaim the (now-empty) vec back into scratch so capacity is retained.
        // del_candidates_batch consumed and dropped it; re-initialise.
        // (Capacity is lost here but to_delete is typically small — acceptable.)

        let s = &mut self.regen_scratch;

        // ── 2. Build merged focal set (particles near any region) ─────────
        s.particle_buf.clear();
        s.focal_particles.clear();
        for &(center, radius) in regions {
            let search_r = radius + self.config.max_cutoff() * 2.0;
            self.particle_grid.query_into(center.x, center.y, search_r, &mut s.particle_buf);
            s.focal_particles.extend_from_slice(&s.particle_buf);
        }
        s.focal_particles.sort_unstable_by(|a, b| {
            a.pos.x.total_cmp(&b.pos.x).then(a.pos.y.total_cmp(&b.pos.y))
        });
        s.focal_particles.dedup_by(|a, b| a.pos == b.pos);

        // ── 3. Build merged overlap-check arrays ──────────────────────────
        s.overlap_particles.clear();
        for &(center, radius) in regions {
            let search_r = radius + self.config.max_cutoff() * 2.0;
            let overlap_r = search_r + 4.0 * self.config.max_radius();
            self.particle_grid.query_into(center.x, center.y, overlap_r, &mut s.particle_buf);
            s.overlap_particles.extend_from_slice(&s.particle_buf);
        }
        s.overlap_particles.sort_unstable_by(|a, b| {
            a.pos.x.total_cmp(&b.pos.x).then(a.pos.y.total_cmp(&b.pos.y))
        });
        s.overlap_particles.dedup_by(|a, b| a.pos == b.pos);

        s.overlap_pos.clear();
        s.overlap_rad.clear();
        s.overlap_idx.clear();
        for (i, p) in s.overlap_particles.iter().enumerate() {
            s.overlap_pos.push(p.pos);
            s.overlap_rad.push(p.radius);
            s.overlap_idx.push(i);
        }

        // ── 4. For each focal particle A and each candidate type C, generate sites.
        let delta = self.config.delta;
        let max_radius = self.config.max_radius();
        // Snapshot focal set to avoid aliasing with &mut self.regen_scratch below.
        let focal_len = self.regen_scratch.focal_particles.len();
        for fi in 0..focal_len {
            let a_pos;
            let a_r;
            {
                let p = &self.regen_scratch.focal_particles[fi];
                a_pos = p.pos;
                a_r = p.radius;
            }

            let bonded_cutoff = a_r + max_radius + delta;
            self.regen_scratch.bonded_dirs.clear();
            self.particle_grid.query_iter(a_pos.x, a_pos.y, bonded_cutoff, |other| {
                if other.pos == a_pos { return; }
                let d = a_pos.distance(other.pos);
                if d < a_r + other.radius + delta {
                    let diff = other.pos - a_pos;
                    self.regen_scratch.bonded_dirs.push(diff / d);
                }
            });

            for type_c in 0..self.config.n_types() {
                let rc = self.config.radius(type_c);

                // ── 4a. Intersection sites ────────────────────────────────
                let pair_search = a_r + rc + max_radius + rc + delta;
                self.regen_scratch.particle_buf.clear();
                self.particle_grid.query_into(a_pos.x, a_pos.y, pair_search, &mut self.regen_scratch.particle_buf);
                let pb_len = self.regen_scratch.particle_buf.len();
                for bi in 0..pb_len {
                    let (b_pos, b_r) = {
                        let b = &self.regen_scratch.particle_buf[bi];
                        if (b.pos.x, b.pos.y) <= (a_pos.x, a_pos.y) { continue; }
                        (b.pos, b.radius)
                    };

                    let r_ac = a_r + rc;
                    let r_bc = b_r + rc;

                    for maybe_pt in circle_intersections(a_pos, r_ac, b_pos, r_bc) {
                        let Some(site_pos) = maybe_pt else { continue };

                        let da = a_pos.distance(site_pos);
                        let db = b_pos.distance(site_pos);
                        if da < r_ac - delta || da > r_ac + delta { continue; }
                        if db < r_bc - delta || db > r_bc + delta { continue; }

                        if site_has_overlap(site_pos, rc, delta,
                            &self.regen_scratch.overlap_idx,
                            &self.regen_scratch.overlap_pos,
                            &self.regen_scratch.overlap_rad,
                        ) { continue; }
                        let orientation = self.best_candidate_orientation(site_pos, type_c);
                        self.add_candidate(CandidateSite { pos: site_pos, type_id: type_c, orientation });
                    }
                }

                // ── 4b. Arc sites ─────────────────────────────────────────
                let n_ang = self.config.num_isolated_angles;
                let exclusion = TAU / n_ang as f32;
                let r_contact = a_r + rc;
                let max_arc = self.config.max_arc_sites_per_type;

                let cos_excl = exclusion.cos();
                self.regen_scratch.arc_candidates.clear();
                for k in 0..n_ang {
                    let theta = (k as f32) * TAU / (n_ang as f32);
                    let theta_v = Vec2::new(theta.cos(), theta.sin());
                    let blocked = self.regen_scratch.bonded_dirs.iter()
                        .any(|&bd| theta_v.dot(bd) > cos_excl);
                    if blocked { continue; }
                    let site_pos = a_pos + Vec2::new(theta.cos(), theta.sin()) * r_contact;
                    if site_has_overlap(site_pos, rc, delta,
                        &self.regen_scratch.overlap_idx,
                        &self.regen_scratch.overlap_pos,
                        &self.regen_scratch.overlap_rad,
                    ) { continue; }
                    self.regen_scratch.arc_candidates.push(site_pos);
                }

                let arc_len = self.regen_scratch.arc_candidates.len();
                if arc_len <= max_arc {
                    for ai in 0..arc_len {
                        let pos = self.regen_scratch.arc_candidates[ai];
                        let orientation = self.best_candidate_orientation(pos, type_c);
                        self.add_candidate(CandidateSite { pos, type_id: type_c, orientation });
                    }
                } else {
                    let step = arc_len as f32 / max_arc as f32;
                    for i in 0..max_arc {
                        let idx = (i as f32 * step) as usize;
                        let pos = self.regen_scratch.arc_candidates[idx];
                        let orientation = self.best_candidate_orientation(pos, type_c);
                        self.add_candidate(CandidateSite { pos, type_id: type_c, orientation });
                    }
                }
            }
        }
    }

    /// Rebuild the flat f32 particle buffer for the WASM export.
    /// Layout per particle: [x, y, type_id, radius, orientation]
    fn rebuild_particle_buf(&mut self) {
        self.particle_buf.clear();
        for p in self.particle_grid.iter() {
            self.particle_buf.push(p.pos.x);
            self.particle_buf.push(p.pos.y);
            self.particle_buf.push(p.type_id as f32);
            self.particle_buf.push(p.radius);
            self.particle_buf.push(p.orientation.y.atan2(p.orientation.x));
        }
    }

    /// Rebuild the flat f32 candidate buffer for the debug overlay.
    /// Layout per candidate: [x, y, type_id, rate]
    fn rebuild_candidate_buf(&mut self) {
        self.candidate_buf.clear();
        for (i, site) in self.candidates.iter().enumerate() {
            let rate = self.attach_rates.get_rate(i);
            self.candidate_buf.push(site.pos.x);
            self.candidate_buf.push(site.pos.y);
            self.candidate_buf.push(site.type_id as f32);
            self.candidate_buf.push(rate as f32);
        }
    }

    pub fn candidate_buf_ptr(&self) -> *const f32 {
        self.candidate_buf.as_ptr()
    }

    pub fn candidate_count(&self) -> u32 {
        self.candidates.len() as u32
    }

    // ── Public getters ────────────────────────────────────────────────────────

    pub fn particle_count(&self) -> u32 {
        self.particle_grid.len() as u32
    }

    pub fn simulation_time(&self) -> f64 {
        self.time
    }

    pub fn print_timers(&self) {
        let ghz = self.cpu_ghz;
        for i in 0..(TimerEntry::Count as usize) {
            let cycles = self.timers[i];
            let calls = self.usages[i];
            let label = ENUM_STR[i];
            let total_ms = cycles as f64 / ghz / 1e6;
            if calls > 0 {
                let avg_ms = total_ms / calls as f64;
                println!("{label:18} calls={calls:6}  avg={avg_ms:8.3}ms  total={total_ms:.3}ms");
            } else {
                println!("{label:18} calls=     0  avg=       -    total=  0.000ms");
            }
        }
    }

    pub fn particle_buf_ptr(&self) -> *const f32 {
        self.particle_buf.as_ptr()
    }

    /// Fixed-grid relaxation variant.
    ///
    /// Instead of dynamically growing an active set, this activates **all** particles
    /// in a 5×5 cell neighbourhood around the new particle immediately, backed by a
    /// `(5 + 2*cell_r)²` initialised border.  `Cell` indices are cached before
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
        let act_half  = 3i64;              // 5×5 active zone
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
            // Zero forces and torques.
            for &idx in &initialized_cell_indices {
                if let Some(fs) = self.particle_grid.cell_forces.get_mut(idx) {
                    for f in fs.iter_mut() { *f = Vec2::ZERO; }
                }
                if let Some(ts) = self.particle_grid.cell_torques.get_mut(idx) {
                    for t in ts.iter_mut() { *t = 0.0; }
                }
            }

            // Accumulate patchy forces and torques via cached cell indices.
            for &ai_slot in &active_slots {
                let a_ind = cell_ind[ai_slot];
                if a_ind.is_none() { continue; }
                let a_ind = a_ind.unwrap();
                let a_cell = &self.particle_grid.cells[a_ind];
                let a_gx = (ai_slot % w) as i64;
                let a_gy = (ai_slot / w) as i64;

                for ai in 0..a_cell.particles.len() {
                    if !a_cell.is_active[ai] { continue; }
                    let a_pos    = a_cell.particles[ai].pos;
                    let a_r      = a_cell.particles[ai].radius;
                    let a_type   = a_cell.particles[ai].type_id;
                    let a_ori    = a_cell.particles[ai].orientation;

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
                                let b_pos  = b_cell.particles[bi].pos;
                                let b_r    = b_cell.particles[bi].radius;
                                let b_type = b_cell.particles[bi].type_id;
                                let b_ori  = b_cell.particles[bi].orientation;
                                let b_frz  = b_cell.is_frozen[bi];

                                let r_contact = a_r + b_r;
                                // Quick distance cull before full patchy computation.
                                let r2 = (b_pos - a_pos).length_squared();
                                let r_cut = r_contact * lj_cutoff_factor;
                                if r2 >= r_cut * r_cut { continue; }

                                let (f, tau_a, tau_b) = patchy_force_torque(
                                    b_pos - a_pos,
                                    r_contact,
                                    &self.config.particle_types[a_type].patches,
                                    a_ori,
                                    &self.config.particle_types[b_type].patches,
                                    b_ori,
                                    &self.config.patch_lut,
                                    self.config.patch_type_count,
                                    self.config.epsilon(a_type, b_type),
                                    self.config.lj_cutoff_factor,
                                );

                                if f != Vec2::ZERO || tau_a != 0.0 {
                                    let af = &mut self.particle_grid.cell_forces[a_ind];
                                    if ai < af.len() { af[ai] += f; }
                                    let at_ = &mut self.particle_grid.cell_torques[a_ind];
                                    if ai < at_.len() { at_[ai] += tau_a; }
                                }

                                if b_is_active || !b_frz {
                                    let bf = &mut self.particle_grid.cell_forces[b_ind];
                                    if bi < bf.len() { bf[bi] -= f; }
                                    let bt = &mut self.particle_grid.cell_torques[b_ind];
                                    if bi < bt.len() { bt[bi] += tau_b; }
                                }
                            }
                        }
                    }
                }
            }

            // Velocity + angular-velocity step (FIRE-style damped).
            let mut converged = true;
            for &cp in &active_cell_indices {
                let cell = &mut self.particle_grid.cells[cp];
                let fp   = &self.particle_grid.cell_forces[cp];
                let tp   = &self.particle_grid.cell_torques[cp];

                for ind in 0..cell.particles.len() {
                    if !cell.is_active[ind] { continue; }

                    // Linear
                    let f_raw     = fp.get(ind).copied().unwrap_or(Vec2::ZERO);
                    let f_clamped = f_raw.clamp_length_max(3.0);
                    let f_mag     = f_clamped.length();
                    let f_eff     = if f_mag > static_friction {
                        f_clamped * ((f_mag - static_friction) / f_mag)
                    } else { Vec2::ZERO };

                    let v_raw = damping * cell.velocities[ind] + alpha * f_eff;
                    let v_new = if v_raw.dot(f_eff) >= 0.0 { v_raw } else { Vec2::ZERO };

                    cell.new_pos[ind]    = cell.particles[ind].pos + v_new;
                    cell.velocities[ind] = v_new;

                    // Angular
                    let tau_raw     = tp.get(ind).copied().unwrap_or(0.0);
                    let tau_clamped = tau_raw.clamp(-3.0, 3.0);
                    let tau_eff     = if tau_clamped.abs() > static_friction { tau_clamped } else { 0.0 };
                    let omega_raw   = damping * cell.ang_velocities[ind] + alpha * tau_eff;
                    let omega_new   = if omega_raw * tau_eff >= 0.0 { omega_raw } else { 0.0 };

                    let (s, c) = omega_new.sin_cos();
                    let ori = cell.particles[ind].orientation;
                    let new_ori = Vec2::new(ori.x * c - ori.y * s, ori.x * s + ori.y * c);
                    // Renormalize to prevent floating-point drift.
                    let len_sq = new_ori.length_squared();
                    cell.new_orientation[ind] = if (len_sq - 1.0).abs() > 1e-6 {
                        new_ori / len_sq.sqrt()
                    } else {
                        new_ori
                    };
                    cell.ang_velocities[ind]  = omega_new;

                    if f_eff.length_squared() >= 0.05 * 0.05 || tau_eff.abs() >= 0.05 {
                        converged = false;
                    }
                }
            }

            let _moved = self.particle_grid.commit_positions(&active_cells);
            if converged { break; }
        }

        // ── Post-relaxation: update rates and regen candidates ────────────────
        // Query the live grid rather than stale cell_ind: particles can shift
        // into cells that didn't exist when cell_ind was built, so a cell_ind
        // walk would miss their final positions and leave stale candidates.
        let delta   = self.config.delta;
        let regen_r = self.config.max_cutoff() * 2.0;
        let act_radius = act_half as f32 * cell_size + cell_size;
        let mut rate_particles: Vec<Particle> = Vec::new();
        let mut regen_regions:  Vec<(Vec2, f32)> = Vec::new();
        self.particle_grid.query_into(new_p.pos.x, new_p.pos.y, act_radius, &mut rate_particles);
        for p in &rate_particles { regen_regions.push((p.pos, regen_r)); }

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

    /// Return a JSON string with per-type metadata (color, radius, patches) for the renderer.
    pub fn type_metadata_json(&self) -> String {
        let entries: Vec<String> = self
            .config
            .particle_types
            .iter()
            .map(|t| {
                let patches: Vec<String> = t.patches.iter().map(|p| {
                    // r##"..."## avoids the `"#` closing the raw string inside "#ffffff".
                    format!(r##"{{"angle_rad":{:.6},"color":"#ffffff"}}"##, p.position_rad)
                }).collect();
                format!(
                    r#"{{"color":"{}","radius":{},"patches":[{}]}}"#,
                    t.color, t.radius, patches.join(",")
                )
            })
            .collect();
        format!("[{}]", entries.join(","))
    }
}
