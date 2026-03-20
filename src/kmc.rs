use std::collections::HashSet;
use std::f32::consts::PI;
use glam::Vec2;

use crate::candidates::CandidateSite;
use crate::config::SimConfig;
use crate::forces::{patchy_force_torque, patchy_pair_energy};
use crate::particle::Particle;
use crate::rng::Rng;
use crate::spatial::ParticleGrid;
use crate::math_utils::exp3; 

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
    Count,
}

static ENUM_STR: [&str; TimerEntry::Count as usize] = [
    "Select", "Attach", "Detach", "Rebuild", "Relax",
];

#[derive(Clone)]
struct TestParticleState {
    particle: Particle,
    velocity: Vec2,
    omega: f32,
    frozen: bool,
}

pub struct Simulation {
    /// Particles stored directly in the spatial grid; indexed by slot.
    pub particle_grid: ParticleGrid,

    pub static_friction: f32,

    pub config: SimConfig,
    pub time: f64,
    pub rng: Rng,

    /// Accumulated TSC cycles per timer bucket.
    pub timers: Vec<u64>,
    pub usages: Vec<u64>,
    cpu_ghz: f64,

    particle_buf: Vec<f32>, // Flat buffer of (x, y, type_id, radius, orientation) for WASM export.

    /// Deterministic test-mode state (only used when `config.testing_mode` is true).
    test_particles: Vec<TestParticleState>,
    test_step_counter: u64,
}

/*
Simulation steps
- Select parent particle proportional to attach_rate_sum, or select particle to detach.
- If attach:
    - Generate candidates on-demand for parent particle, sample one proportional to rate.
    - Add particle, run relaxation, update neighbor attach/detach rates.
- If detach:
    - Remove particle, update neighbor attach/detach rates.
*/

impl Simulation {
    pub fn new(mut config: SimConfig) -> Self {
        config.init_cache();
        let cell_size = config.max_cutoff().max(0.1)*1.1;

        let mut sim = Self {
            particle_grid: ParticleGrid::new(cell_size),
            static_friction: config.static_friction,
            time: 0.0,
            rng: Rng::new(config.seed),
            config,
            timers: vec![0u64; TimerEntry::Count as usize],
            usages: vec![0u64; TimerEntry::Count as usize],
            cpu_ghz: calibrate_tsc_ghz(),
            particle_buf: Vec::new(),
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
                // Compute initial attach rate for the seed particle.
                let key = sim.particle_grid.cell_key(0.0, 0.0);
                if let Some(&cell_idx) = sim.particle_grid.cell_map.get(&key) {
                    if !sim.particle_grid.cells[cell_idx].particles.is_empty() {
                        sim.update_attach_rate(cell_idx, 0);
                    }
                }
            }
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
            let attach_total = self.particle_grid.total_attach_rate();
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
                let result = self.particle_grid.sample_attach_rate(attach_u);
                self.timers[TimerEntry::Select as usize] += rdtsc().wrapping_sub(t0_sel);
                self.usages[TimerEntry::Select as usize] += 1;

                if let Some((cell, idx)) = result {
                    let t0 = rdtsc();
                    self.attach_from_particle(cell, idx);
                    self.timers[TimerEntry::Attach as usize] += rdtsc().wrapping_sub(t0);
                    self.usages[TimerEntry::Attach as usize] += 1;
                }
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

    // ── Attach ───────────────────────────────────────────────────────────────

    /// Select and attach a particle based on a parent particle's patch candidates.
    fn attach_from_particle(&mut self, cell: usize, idx: usize) {
        let candidates = self.generate_candidates_for_particle(cell, idx);
        if candidates.is_empty() {
            // Rate was stale; zero it out so we don't keep sampling this particle.
            self.particle_grid.cells[cell].aggr_attach_rate -=
                self.particle_grid.cells[cell].attach_rates[idx];
            self.particle_grid.cells[cell].attach_rates[idx] = 0.0;
            return;
        }

        // Weighted sample from candidates proportional to rate.
        let total: f64 = candidates.iter().map(|(_, r)| r).sum();
        let u = self.rng.next_f64() * total;
        let mut accum = 0.0f64;
        let mut chosen = candidates[0].0.clone();
        for (site, rate) in &candidates {
            accum += rate;
            if u < accum {
                chosen = site.clone();
                break;
            }
        }

        self.attach_inner(chosen);
    }

    /// Place a particle at `site`, optionally relax, then update neighbor rates.
    fn attach_inner(&mut self, site: CandidateSite) {
        let type_c = site.type_id;
        let rc = self.config.radius(type_c);
        let pos = site.pos;

        let query_r = rc + self.config.max_radius() + 1e-6;
        let delta = self.config.delta;
        let mut has_overlap = false;
        self.particle_grid.query_iter(pos.x, pos.y, query_r, |p| {
            if has_overlap { return; }
            let overlap_dist = rc + p.radius - delta;
            if overlap_dist > 0.0 && pos.distance_squared(p.pos) < overlap_dist * overlap_dist {
                has_overlap = true;
            }
        });
        if has_overlap {
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
            self.update_neighbors_attach_rates(pos);
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

        self.update_neighbors_attach_rates(pos);
    }

    // ── Internal helpers ─────────────────────────────────────────────────────

    /// Push a particle into the data structures without updating rates.
    fn add_particle_raw(&mut self, p: Particle) {
        self.particle_grid.insert(p, 0.0); // placeholder rate; caller must update
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

    /// Energy of a candidate particle of type `type_c` placed at `pos` with a
    /// fixed `orientation`, summed over all nearby placed particles.
    /// Returns a negative value when the net interaction is attractive.
    fn site_energy_for_orientation(&self, pos: Vec2, type_c: usize, orientation: Vec2) -> f32 {
        let config = &self.config;
        let rc = config.radius(type_c);
        let query_r = rc + config.max_radius() + config.delta;
        let nearby = self.particle_grid.query(pos.x, pos.y, query_r);
        let mut energy = 0.0f32;
        for nb in &nearby {
            let r_vec = nb.pos - pos;
            let r_contact = rc + nb.radius;
            energy += patchy_pair_energy(
                r_vec, r_contact,
                &config.particle_types[type_c].patches, orientation,
                &config.particle_types[nb.type_id].patches, nb.orientation,
                &config.patch_lut, config.patch_type_count,
                config.epsilon(type_c, nb.type_id), config.delta,
            );
        }
        energy
    }

    /// Compute total attachment rate attributable to particle at (cell, idx).
    /// For each patch P_A, for each type C, find the best matching patch on C
    /// and accumulate attach_rate if the interaction is attractive.
    fn compute_attach_rate_sum(&self, cell: usize, idx: usize) -> f64 {
        let (a_pos, a_r, a_type, a_ori) = {
            let p = &self.particle_grid.cells[cell].particles[idx];
            (p.pos, p.radius, p.type_id, p.orientation)
        };

        let config = &self.config;
        let patch_tc = config.patch_type_count;
        let n_a_patches = config.particle_types[a_type].patches.len();

        // Non-patchy particles don't generate candidates in this scheme.
        if n_a_patches == 0 || patch_tc == 0 {
            return 0.0;
        }

        let delta = config.delta;
        let mut total = 0.0f64;

        for pa_idx in 0..n_a_patches {
            let (pa_type_id, pa_cs) = {
                let p_a = &config.particle_types[a_type].patches[pa_idx];
                (p_a.patch_type_id, p_a.position_cs)
            };
            // Lab-frame direction of this patch: a_ori * pa_cs (complex multiply).
            let phi = Vec2::new(
                a_ori.x * pa_cs.x - a_ori.y * pa_cs.y,
                a_ori.x * pa_cs.y + a_ori.y * pa_cs.x,
            );

            for type_c in 0..config.n_types() {
                let (rc, n_c_patches) = {
                    let ct = &config.particle_types[type_c];
                    (ct.radius, ct.patches.len())
                };
                let c_pos = a_pos + phi * (a_r + rc);

                // Overlap check via direct grid query.
                let overlap_query_r = rc + config.max_radius() + delta;
                let mut has_overlap = false;
                self.particle_grid.query_iter(c_pos.x, c_pos.y, overlap_query_r, |other| {
                    if has_overlap { return; }
                    let min_sep = rc + other.radius - delta;
                    if min_sep > 0.0 && c_pos.distance_squared(other.pos) < min_sep * min_sep {
                        has_overlap = true;
                    }
                });
                if has_overlap { continue; }

                // Find the best matching patch on C.
                let mut best_energy = 0.0f32;
                let mut best_c_ori = Vec2::X;
                let mut found = false;

                for pc_idx in 0..n_c_patches {
                    let (pc_type_id, pc_cs) = {
                        let p_c = &config.particle_types[type_c].patches[pc_idx];
                        (p_c.patch_type_id, p_c.position_cs)
                    };
                    if pa_type_id >= patch_tc || pc_type_id >= patch_tc { continue; }

                    let lut = &config.patch_lut[pa_type_id * patch_tc + pc_type_id];
                    // epsilon < 0 means attractive; skip non-attractive pairs.
                    if !lut.enabled || lut.epsilon >= 0.0 { continue; }

                    // Orientation so P_C points toward A (direction -phi).
                    let neg_phi = -phi;
                    let c_ori = Vec2::new(
                        neg_phi.x * pc_cs.x + neg_phi.y * pc_cs.y,
                        neg_phi.y * pc_cs.x - neg_phi.x * pc_cs.y,
                    );

                    let energy = self.site_energy_for_orientation(c_pos, type_c, c_ori);
                    if energy < best_energy {
                        best_energy = energy;
                        best_c_ori = c_ori;
                        found = true;
                    }
                }

                if !found { continue; }
                let _ = best_c_ori; // used in generate_candidates_for_particle
                total += config.attach_rate(type_c, -best_energy) as f64;
            }
        }

        total
    }

    /// Generate on-demand candidates for the particle at (cell, idx).
    /// Returns Vec<(CandidateSite, rate)>.
    fn generate_candidates_for_particle(&self, cell: usize, idx: usize) -> Vec<(CandidateSite, f64)> {
        let (a_pos, a_r, a_type, a_ori) = {
            let p = &self.particle_grid.cells[cell].particles[idx];
            (p.pos, p.radius, p.type_id, p.orientation)
        };

        let config = &self.config;
        let patch_tc = config.patch_type_count;
        let n_a_patches = config.particle_types[a_type].patches.len();

        if n_a_patches == 0 || patch_tc == 0 {
            return Vec::new();
        }

        let delta = config.delta;
        let mut result: Vec<(CandidateSite, f64)> = Vec::new();

        for pa_idx in 0..n_a_patches {
            let (pa_type_id, pa_cs) = {
                let p_a = &config.particle_types[a_type].patches[pa_idx];
                (p_a.patch_type_id, p_a.position_cs)
            };
            let phi = Vec2::new(
                a_ori.x * pa_cs.x - a_ori.y * pa_cs.y,
                a_ori.x * pa_cs.y + a_ori.y * pa_cs.x,
            );

            for type_c in 0..config.n_types() {
                let (rc, n_c_patches) = {
                    let ct = &config.particle_types[type_c];
                    (ct.radius, ct.patches.len())
                };
                let c_pos = a_pos + phi * (a_r + rc);

                let overlap_query_r = rc + config.max_radius() + delta;
                let mut has_overlap = false;
                self.particle_grid.query_iter(c_pos.x, c_pos.y, overlap_query_r, |other| {
                    if has_overlap { return; }
                    let min_sep = rc + other.radius - delta;
                    if min_sep > 0.0 && c_pos.distance_squared(other.pos) < min_sep * min_sep {
                        has_overlap = true;
                    }
                });
                if has_overlap { continue; }

                let mut best_energy = 0.0f32;
                let mut best_c_ori = Vec2::X;
                let mut found = false;

                for pc_idx in 0..n_c_patches {
                    let (pc_type_id, pc_cs) = {
                        let p_c = &config.particle_types[type_c].patches[pc_idx];
                        (p_c.patch_type_id, p_c.position_cs)
                    };
                    if pa_type_id >= patch_tc || pc_type_id >= patch_tc { continue; }

                    let lut = &config.patch_lut[pa_type_id * patch_tc + pc_type_id];
                    if !lut.enabled || lut.epsilon >= 0.0 { continue; }

                    let neg_phi = -phi;
                    let c_ori = Vec2::new(
                        neg_phi.x * pc_cs.x + neg_phi.y * pc_cs.y,
                        neg_phi.y * pc_cs.x - neg_phi.x * pc_cs.y,
                    );

                    let energy = self.site_energy_for_orientation(c_pos, type_c, c_ori);
                    if energy < best_energy {
                        best_energy = energy;
                        best_c_ori = c_ori;
                        found = true;
                    }
                }

                if !found { continue; }
                let rate = config.attach_rate(type_c, -best_energy) as f64;
                result.push((CandidateSite { pos: c_pos, type_id: type_c, orientation: best_c_ori }, rate));
            }
        }

        result
    }

    /// Recompute and store the attach rate sum for the particle at (cell, idx).
    fn update_attach_rate(&mut self, cell: usize, idx: usize) {
        let rate = self.compute_attach_rate_sum(cell, idx);
        let old = self.particle_grid.cells[cell].attach_rates[idx];
        self.particle_grid.cells[cell].attach_rates[idx] = rate;
        self.particle_grid.cells[cell].aggr_attach_rate += rate - old;
    }

    /// Find the particle at `pos` in the grid and update its attach rate.
    fn update_attach_rate_for_particle_at(&mut self, pos: Vec2) {
        let key = self.particle_grid.cell_key(pos.x, pos.y);
        let cell_idx = match self.particle_grid.cell_map.get(&key) {
            Some(&i) => i,
            None => return,
        };
        let part_idx = match self.particle_grid.cells[cell_idx].particles.iter().position(|p| p.pos == pos) {
            Some(i) => i,
            None => return,
        };
        self.update_attach_rate(cell_idx, part_idx);
    }

    /// Update attach rates for all particles within interaction range of `pos`.
    fn update_neighbors_attach_rates(&mut self, pos: Vec2) {
        let radius = self.config.max_cutoff() * 2.0 + self.config.max_radius();
        let neighbor_pos: Vec<Vec2> = self.particle_grid
            .query(pos.x, pos.y, radius)
            .into_iter()
            .map(|p| p.pos)
            .collect();
        for npos in neighbor_pos {
            self.update_attach_rate_for_particle_at(npos);
        }
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

    /// Return a JSON string with per-type metadata (color, radius, patches) for the renderer.
    pub fn type_metadata_json(&self) -> String {
        let entries: Vec<String> = self
            .config
            .particle_types
            .iter()
            .map(|t| {
                let patches: Vec<String> = t.patches.iter().map(|p| {
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

    /// Fixed-grid relaxation variant.
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
        let act_half  = 1i64;              // 5×5 active zone
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

        // ── Post-relaxation: update detach rates and attach rate sums ─────────
        let delta     = self.config.delta;
        let act_radius = act_half as f32 * cell_size + cell_size;
        let mut rate_particles: Vec<Particle> = Vec::new();
        self.particle_grid.query_into(new_p.pos.x, new_p.pos.y, act_radius, &mut rate_particles);

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

        // Update attach rates for all particles in the active zone and neighbours.
        let update_r = act_radius + self.config.max_cutoff() * 2.0;
        let to_update_pos: Vec<Vec2> = self.particle_grid
            .query(new_p.pos.x, new_p.pos.y, update_r)
            .into_iter()
            .map(|p| p.pos)
            .collect();
        for pos in to_update_pos {
            self.update_attach_rate_for_particle_at(pos);
        }

        let initialized: HashSet<(i64, i64)> = keys.iter().copied().collect();
        self.particle_grid.clear_physics_for_cells(&initialized);
    }
}
