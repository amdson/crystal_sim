use std::f64::consts::TAU;
use glam::DVec2;

use crate::candidates::{CandidateSite, circle_intersections, site_has_overlap};
use crate::config::SimConfig;
use crate::forces::lj_force_vec;
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

pub struct Simulation {
    /// Particles stored directly in the spatial grid; indexed by slot.
    pub particles: ParticleGrid,
    pub detach_rates: RateCatalog, // indexed by particle slot

    /// candidates[type_id] = current valid attachment sites for that type
    pub candidates: Vec<CandidateSite>,
    pub candidates_spatial: SpatialHash,
    pub attach_rates: RateCatalog, // indexed by candidate index

    pub static_friction: f64,

    pub config: SimConfig,
    pub time: f64,
    pub rng: Rng,
    /// Flat f32 particle buffer: [x, y, type_id_f32, radius, ...] stride 4
    pub particle_buf: Vec<f32>,
    /// Accumulated TSC cycles per timer bucket.
    pub timers: Vec<u64>,
    pub usages: Vec<u64>,
    cpu_ghz: f64,
    /// Reusable boolean mask for active-set relaxation (avoids per-call allocation).
    in_active: Vec<bool>,
    /// Reusable buffer for spatial hash queries (avoids per-query allocation).
    query_buf: Vec<usize>,
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
        let cell_size = config.max_cutoff().max(0.1) * 2.0; //Large cell size

        let mut sim = Self {
            particles: ParticleGrid::new(cell_size),
            detach_rates: RateCatalog::new(),
            candidates: Vec::new(),
            candidates_spatial: SpatialHash::new(cell_size),
            attach_rates: RateCatalog::new(),
            static_friction: config.static_friction,
            time: 0.0,
            rng: Rng::new(config.seed),
            particle_buf: Vec::new(),
            config,
            timers: vec![0u64; TimerEntry::Count as usize],
            usages: vec![0u64; TimerEntry::Count as usize],
            cpu_ghz: calibrate_tsc_ghz(),
            in_active: Vec::new(),
            query_buf: Vec::new(),
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
            let detach_total = self.detach_rates.total();
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
                let slot = self.detach_rates.select(detach_u);
                self.timers[TimerEntry::Select as usize] += rdtsc().wrapping_sub(t0_sel);
                self.usages[TimerEntry::Select as usize] += 1;

                let t0 = rdtsc();
                self.detach(slot);
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
                let avg_ns = cycles as f64 / ghz / calls as f64;
                println!("{label:18} calls={calls:6}  avg={avg_ns:8.0}ns  total={total_ms:.3}ms");
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
        let nearby = self.particles.query(pos.x, pos.y, query_r);
        let positions: Vec<DVec2> = nearby.iter().map(|&i| self.particles[i].pos).collect();
        let radii: Vec<f64> = nearby.iter().map(|&i| self.particles[i].radius).collect();
        if site_has_overlap(pos, rc, self.config.delta, &(0..nearby.len()).collect::<Vec<_>>(), &positions, &radii) {
            return;
        }

        self.add_particle_raw(Particle::new(pos.x, pos.y, type_c, rc));

        let new_idx = self.particles.len() - 1;

        if self.config.relax_steps > 0 {
            let t0 = rdtsc();
            self.relax_new_particle(new_idx);
            self.timers[TimerEntry::Relax as usize] += rdtsc().wrapping_sub(t0);
            self.usages[TimerEntry::Relax as usize] += 1;
        } else {
            self.update_detach_rate(new_idx);
            let nbs = self.get_neighbors(new_idx);
            for nb in nbs {
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

    fn detach(&mut self, slot: usize) {
        if slot >= self.particles.len() {
            return;
        }

        let pos = self.particles[slot].pos;
        let bonded = self.get_neighbors(slot);
        let last = self.particles.len() - 1;

        // swap_remove on ParticleGrid internally handles cell membership fixup;
        // swap_remove_rate mirrors the same slot movement on the rate catalog.
        self.particles.swap_remove(slot);
        self.detach_rates.swap_remove_rate(slot);

        // Old neighbours lost a bond → update their detach rates.
        // The particle formerly at `last` is now at `slot`.
        for nb in bonded {
            let actual = if nb == last && slot != last { slot } else { nb };
            if actual < self.particles.len() {
                self.update_detach_rate(actual);
            }
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
        self.particles.insert(p);
        self.detach_rates.add_rate(0.0); // placeholder; caller must call update_detach_rate
    }

    /// Compute E_bind for particle at `slot` and refresh its detachment rate.
    pub fn update_detach_rate(&mut self, slot: usize) {
        if slot >= self.particles.len() {
            return;
        }
        let e_bind = self.binding_energy(slot);
        let rate = self.config.detach_rate(e_bind);
        self.detach_rates.set_rate(slot, rate);
    }

    pub fn get_neighbors(&self, slot: usize) -> Vec<usize> {
        if slot >= self.particles.len() {
            return Vec::new();
        }
        let p = &self.particles[slot];
        let cutoff = p.radius + self.config.max_radius() + self.config.delta;
        let nearby = self.particles.query(p.pos.x, p.pos.y, cutoff);
        nearby.into_iter()
            .filter(|&ni| ni != slot && ni < self.particles.len() && p.bonds_to(&self.particles[ni], self.config.delta))
            .collect()
    }

    fn binding_energy(&self, slot: usize) -> f64 {
        let type_i = self.particles[slot].type_id;
        self.get_neighbors(slot)
            .iter()
            .map(|&nb| self.config.epsilon(type_i, self.particles[nb].type_id))
            .sum()
    }

    fn site_binding_energy(&self, site: &CandidateSite) -> f64 {
        let type_c = site.type_id;
        let rc = self.config.radius(type_c);
        let pos = site.pos;
        let query_r = rc + self.config.max_radius() + self.config.delta;
        let nearby = self.particles.query(pos.x, pos.y, query_r);
        nearby
            .iter()
            .filter(|&&i| {
                let d = self.particles[i].pos.distance(pos);
                let contact = self.particles[i].radius + rc;
                d >= contact - self.config.delta && d <= contact + self.config.delta
            })
            .map(|&i| self.config.epsilon(type_c, self.particles[i].type_id))
            .sum()
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
        let mut focal: Vec<usize> = Vec::new();
        for &(center, radius) in regions {
            let search_r = radius + self.config.max_cutoff() * 2.0;
            self.particles.query_into(center.x, center.y, search_r, &mut self.query_buf);
            for &i in &self.query_buf {
                if i < self.particles.len() {
                    focal.push(i);
                }
            }
        }
        focal.sort_unstable();
        focal.dedup();

        // ── 3. Build merged overlap-check arrays ──────────────────────────
        let mut overlap_set: Vec<usize> = Vec::new();
        for &(center, radius) in regions {
            let search_r = radius + self.config.max_cutoff() * 2.0;
            let overlap_r = search_r + 4.0 * self.config.max_radius();
            self.particles.query_into(center.x, center.y, overlap_r, &mut self.query_buf);
            for &i in &self.query_buf {
                if i < self.particles.len() {
                    overlap_set.push(i);
                }
            }
        }
        overlap_set.sort_unstable();
        overlap_set.dedup();
        let overlap_pos: Vec<DVec2> = overlap_set.iter().map(|&i| self.particles[i].pos).collect();
        let overlap_rad: Vec<f64> = overlap_set.iter().map(|&i| self.particles[i].radius).collect();
        let overlap_idx: Vec<usize> = (0..overlap_set.len()).collect();

        // ── 4. For each focal particle A and each candidate type C, generate sites.
        for &a_idx in &focal {
            let a_pos = self.particles[a_idx].pos;
            let a_r = self.particles[a_idx].radius;

            let bonded_cutoff = a_r + self.config.max_radius() + self.config.delta;
            self.particles.query_into(a_pos.x, a_pos.y, bonded_cutoff, &mut self.query_buf);
            let bonded_angles: Vec<f64> = self.query_buf.iter()
                .filter_map(|&ni| {
                    if ni == a_idx || ni >= self.particles.len() { return None; }
                    let other = &self.particles[ni];
                    let d = a_pos.distance(other.pos);
                    if d < a_r + other.radius + self.config.delta {
                        Some(f64::atan2(other.pos.y - a_pos.y, other.pos.x - a_pos.x))
                    } else { None }
                })
                .collect();

            for type_c in 0..self.config.n_types() {
                let rc = self.config.radius(type_c);

                // ── 4a. Intersection sites ────────────────────────────────
                let pair_search = a_r + rc + self.config.max_radius() + rc + self.config.delta;
                self.particles.query_into(a_pos.x, a_pos.y, pair_search, &mut self.query_buf);
                let nearby_indices: Vec<usize> = self.query_buf.clone();

                for &b_idx in &nearby_indices {
                    if b_idx <= a_idx || b_idx >= self.particles.len() {
                        continue;
                    }
                    let b_pos = self.particles[b_idx].pos;
                    let b_r = self.particles[b_idx].radius;

                    let r_ac = a_r + rc;
                    let r_bc = b_r + rc;

                    for maybe_pt in circle_intersections(a_pos, r_ac, b_pos, r_bc) {
                        let Some(site_pos) = maybe_pt else { continue };

                        let da = a_pos.distance(site_pos);
                        let db = b_pos.distance(site_pos);
                        if da < (a_r + rc) - self.config.delta || da > a_r + rc + self.config.delta {
                            continue;
                        }
                        if db < (b_r + rc) - self.config.delta || db > b_r + rc + self.config.delta {
                            continue;
                        }

                        if site_has_overlap(site_pos, rc, self.config.delta, &overlap_idx, &overlap_pos, &overlap_rad) {
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
                    if site_has_overlap(site_pos, rc, self.config.delta, &overlap_idx, &overlap_pos, &overlap_rad) {
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
        let n = self.particles.len();
        self.particle_buf.resize(n * 4, 0.0);
        for slot in 0..n {
            let p = &self.particles[slot];
            let base = slot * 4;
            self.particle_buf[base] = p.pos.x as f32;
            self.particle_buf[base + 1] = p.pos.y as f32;
            self.particle_buf[base + 2] = p.type_id as f32;
            self.particle_buf[base + 3] = p.radius as f32;
        }
    }

    // ── Public getters ────────────────────────────────────────────────────────

    pub fn particle_count(&self) -> u32 {
        self.particles.len() as u32
    }

    pub fn simulation_time(&self) -> f64 {
        self.time
    }

    pub fn particle_buf_ptr(&self) -> *const f32 {
        self.particle_buf.as_ptr()
    }

    pub fn set_temperature(&mut self, t: f64) {
        self.config.temperature = t.max(1e-9);
        for slot in 0..self.particles.len() {
            self.update_detach_rate(slot);
        }
    }

    // ── Post-attachment relaxation ────────────────────────────────────────────

    /// Steepest-descent relaxation of the neighbourhood of a newly attached particle.
    fn relax_new_particle(&mut self, new_idx: usize) {
        if self.config.relax_steps == 0 {
            return;
        }
        let alpha = self.config.relax_alpha;
        let damping = self.config.relax_damping;
        let static_friction = self.static_friction;
        let max_active: usize = 256;
        let lj_cutoff_factor = self.config.lj_cutoff_factor;
        let max_radius_scaled = self.config.max_radius() * lj_cutoff_factor + 1e-6;

        let n_particles = self.particles.len();
        let mut active: Vec<usize> = vec![new_idx];
        if self.in_active.len() < n_particles {
            self.in_active.resize(n_particles, false);
        }
        self.in_active[new_idx] = true;

        let mut initial_positions: Vec<DVec2> = vec![self.particles[new_idx].pos];

        let mut old_positions: Vec<DVec2> = Vec::new();
        let mut forces: Vec<DVec2> = vec![DVec2::ZERO; max_active];
        let mut velocities: Vec<DVec2> = vec![DVec2::ZERO; max_active];
        let mut touched: Vec<(usize, DVec2)> = Vec::new();
        let mut agg: Vec<(usize, DVec2)> = Vec::new();

        for _ in 0..self.config.relax_steps {
            let n_active = active.len();

            // Snapshot positions (Jacobi-style: use old positions for all force calcs)
            old_positions.clear();
            old_positions.extend(active[..n_active].iter().map(|&idx| self.particles[idx].pos));

            // Compute forces on each active particle
            touched.clear();

            for (ai, &idx) in active[..n_active].iter().enumerate() {
                let pos = self.particles[idx].pos;
                let ri = self.particles[idx].radius;

                let query_r = ri + max_radius_scaled;
                self.particles.query_into(pos.x, pos.y, query_r, &mut self.query_buf);

                let mut force = DVec2::ZERO;

                for j_pos in 0..self.query_buf.len() {
                    let j = self.query_buf[j_pos];
                    if j == idx || j >= n_particles {
                        continue;
                    }

                    let r_vec = self.particles[j].pos - pos;
                    let r_contact = ri + self.particles[j].radius;
                    let r_cut = r_contact * lj_cutoff_factor;

                    let f = -lj_force_vec(r_vec, r_contact, 1.0, r_cut);
                    if f != DVec2::ZERO {
                        force += f;
                        if !self.in_active[j] {
                            touched.push((j, -f));
                        }
                    }
                }

                force = force.clamp_length(0.0, 3.0);
                let f_mag = force.length();
                if f_mag > static_friction {
                    force *= (f_mag - static_friction) / f_mag;
                } else {
                    force = DVec2::ZERO;
                }

                forces[ai] = force;
            }

            // Apply damped-velocity Jacobi moves; use move_to to update positions
            // and cell membership atomically.
            let mut converged = true;
            for i in 0..n_active {
                let idx = active[i];

                velocities[i] = damping * velocities[i] + alpha * forces[i];
                if velocities[i].dot(forces[i]) < 0.0 {
                    velocities[i] = DVec2::ZERO;
                }

                let new_pos = old_positions[i] + velocities[i];
                self.particles.move_to(idx, new_pos);

                converged = converged
                    && (forces[i].length_squared() < static_friction * static_friction)
                    && (velocities[i].length_squared() < 1e-12);
            }

            // Activate neighbours whose accumulated reaction force exceeds friction.
            if active.len() < max_active {
                agg.clear();
                for &(nb, f) in &touched {
                    if let Some(entry) = agg.iter_mut().find(|(idx, _)| *idx == nb) {
                        entry.1 += f;
                    } else {
                        agg.push((nb, f));
                    }
                }
                for &(nb_idx, nb_force) in &agg {
                    if active.len() >= max_active {
                        break;
                    }
                    let f_mag = nb_force.length();
                    if f_mag > static_friction {
                        let old_pos = self.particles[nb_idx].pos;
                        initial_positions.push(old_pos);
                        let adjusted = nb_force * ((f_mag - static_friction) / f_mag);
                        let vel = alpha * adjusted;
                        self.particles.move_to(nb_idx, old_pos + vel);
                        self.in_active[nb_idx] = true;
                        active.push(nb_idx);
                        velocities[active.len() - 1] = vel;
                        converged = false;
                    }
                }
            }

            if converged {
                break;
            }
        }

        // ── Post-relaxation: refresh rates and candidate sites ────────────
        let delta = self.config.delta;
        let max_radius = self.config.max_radius();
        let mut rate_indices: Vec<usize> = Vec::new();
        for &idx in &active {
            rate_indices.push(idx);
            let p = &self.particles[idx];
            let cutoff = p.radius + max_radius + delta;
            let nearby = self.particles.query(p.pos.x, p.pos.y, cutoff);
            for ni in nearby {
                if ni != idx && ni < self.particles.len()
                    && p.bonds_to(&self.particles[ni], delta)
                {
                    rate_indices.push(ni);
                }
            }
        }
        rate_indices.sort_unstable();
        rate_indices.dedup();
        for &idx in &rate_indices {
            self.update_detach_rate(idx);
        }

        let regen_r = self.config.max_cutoff() * 2.0;
        let mut regen_regions: Vec<(DVec2, f64)> =
            Vec::with_capacity(active.len() + initial_positions.len());
        for &idx in &active {
            regen_regions.push((self.particles[idx].pos, regen_r));
        }
        for &pos in &initial_positions {
            regen_regions.push((pos, regen_r));
        }
        self.regen_candidates_batch(&regen_regions);

        for &idx in &active {
            self.in_active[idx] = false;
        }
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
