use std::f64::consts::TAU;
use glam::DVec2;

use crate::candidates::{CandidateSite, circle_intersections, site_has_overlap};
use crate::config::SimConfig;
use crate::particle::Particle;
use crate::rates::RateCatalog;
use crate::rng::Rng;
use crate::spatial::SpatialHash;

use std::time::{Duration, Instant};

#[repr(usize)] // Optional: Guarantees the underlying type is usize and ensures sequential values
enum TimerEntry {
    Select,
    Attach,
    Detach,
    Rebuild,
    Count, // Sentinel value to track the number of entries
}

static ENUM_STR: [&str; TimerEntry::Count as usize] = ["Select", "Attach", "Detach", "Rebuild"]; 

pub struct Simulation {
    pub particles: Vec<Particle>, // Flat list of all particles, indexed by particle slot
    pub spatial: SpatialHash, // Spatial hash for particles, indexed by particle slot
    pub detach_rates: RateCatalog, //Detach rates indexed by particle slot

    /// candidates[type_id] = current valid attachment sites for that type
    pub candidates: Vec<CandidateSite>, // Flat list of all candidate sites, indexed by candidate index
    pub candidates_spatial: SpatialHash, // Spatial hash for candidates, indexed by candidate index
    pub attach_rates: RateCatalog, //Attach rates indexed by candidate index

    pub static_friction: f64,

    pub config: SimConfig,
    pub time: f64,
    pub rng: Rng,
    /// Flat f32 particle buffer: [x, y, type_id_f32, radius, ...] stride 4
    pub particle_buf: Vec<f32>,
    pub timers : Vec<Duration>, 
    pub usages: Vec<u32>,
    /// Reusable boolean mask for active-set relaxation (avoids per-call allocation).
    in_active: Vec<bool>,
    /// Reusable buffer for spatial hash queries (avoids per-query allocation).
    query_buf: Vec<usize>,
}

/*
Simulation steps
- Select candidate (add particle t at x, y or / remove particle i)
- If add: 
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
        let cell_size = config.max_cutoff().max(0.1);

        let mut sim = Self {
            particles: Vec::new(), 
            spatial: SpatialHash::new(cell_size), 
            detach_rates: RateCatalog::new(),
            candidates: Vec::new(), 
            candidates_spatial: SpatialHash::new(cell_size),  
            attach_rates: RateCatalog::new(), 
            static_friction: config.static_friction,
            time: 0.0,
            rng: Rng::new(config.seed),  
            particle_buf: Vec::new(), 
            config, 
            timers: (0..TimerEntry::Count as usize).map(|_i| Duration::new(0, 0)).collect(), 
            usages: vec![0; TimerEntry::Count as usize],
            in_active: Vec::new(),
            query_buf: Vec::new(),
        };

        // Seed crystal: one particle of type zero placed at the origin,
        // offset by contact distance so they bond immediately.
        if sim.config.n_types() >= 1 {
            let r0 = sim.config.radius(0); 
            sim.add_particle_raw(Particle::new(0.0, 0.0, 0, r0));  
        }

        // Build initial candidate sites and rates
        let regen_r = sim.config.max_cutoff() * 4.0;
        sim.regen_candidates_near(DVec2::ZERO, regen_r);
        sim
    }

    // ── KMC step loop ────────────────────────────────────────────────────────

    pub fn step(&mut self, n: u32) { 
        for _ in 0..n {
            let now_sel = Instant::now();
            let attach_total = self.attach_rates.total();
            let detach_total = self.detach_rates.total();
            let r_total = attach_total + detach_total;
            if r_total <= 0.0 {
                break;
            }

            let u1 = self.rng.next_f64();
            let u2 = self.rng.next_f64();

            // BKL time advance: Δt = -ln(u2) / R_total
            self.time += -(u2.ln()) / r_total;

            if u1 < attach_total / r_total {
                // Attachment event
                let attach_u = u1 * r_total / attach_total;
                let candidate_ind = self.attach_rates.select(attach_u);
                self.timers[TimerEntry::Select as usize] += now_sel.elapsed();
                self.usages[TimerEntry::Select as usize] += 1;

                let now = Instant::now();
                let site = self.candidates[candidate_ind].clone();
                self.attach(site);

                self.timers[TimerEntry::Attach as usize] += now.elapsed(); 
                self.usages[TimerEntry::Attach as usize] += 1; 
            } else {
                // Detachment event
                let detach_u = (u1 - attach_total / r_total) * r_total / detach_total;
                let slot = self.detach_rates.select(detach_u);
                self.timers[TimerEntry::Select as usize] += now_sel.elapsed();
                self.usages[TimerEntry::Select as usize] += 1;

                let now = Instant::now();
                self.detach(slot);
                self.timers[TimerEntry::Detach as usize] += now.elapsed(); 
                self.usages[TimerEntry::Detach as usize] += 1; 

            }
        }

        let now = Instant::now(); 
        self.rebuild_particle_buf();
        self.timers[TimerEntry::Rebuild as usize] += now.elapsed(); 
        self.usages[TimerEntry::Rebuild as usize] += 1; 
        for i in 0..(TimerEntry::Count as usize) {
            let time = self.timers[i]; 
            let label = ENUM_STR[i]; 
            if self.usages[i] > 0 {
                let avg_time = self.timers[i] / self.usages[i]; 
                println!("{label} avg time: {avg_time:?} total time: {time:?}");
            } else {
                println!("{label} avg time: N/A total time: {time:?}");
            }
        }
    }


    fn del_candidate(&mut self, site_idx: usize) {
        if site_idx >= self.candidates.len() {
            return;
        }
        let site = &self.candidates[site_idx];
        self.candidates_spatial.remove(site_idx, site.pos.x, site.pos.y);
        self.candidates.swap_remove(site_idx); // last element moves into site_idx
        self.attach_rates.swap_remove_rate(site_idx); // same swap-remove in rate catalog to keep indices aligned

        if site_idx < self.candidates.len() {
            // If we swapped, we need to update the spatial hash for the moved candidate
            let moved_site = &self.candidates[site_idx];
            self.candidates_spatial.remove(self.candidates.len(), moved_site.pos.x, moved_site.pos.y); // remove old index
            self.candidates_spatial.insert(site_idx, moved_site.pos.x, moved_site.pos.y); // insert new index
        }
    }

    /// Batch-delete candidate sites. Sorts indices descending so that
    /// swap_remove never invalidates a yet-to-be-deleted index.
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

        // Final hard-core check — site may be stale after previous events
        let query_r = rc + self.config.max_radius() + 1e-6;
        let nearby = self.spatial.query(pos.x, pos.y, query_r);
        // TODO change to filter syntax
        let positions: Vec<DVec2> = nearby.iter().map(|&i| self.particles[i].pos).collect();
        let radii: Vec<f64> = nearby.iter().map(|&i| self.particles[i].radius).collect();
        if site_has_overlap(pos, rc, self.config.delta, &(0..nearby.len()).collect::<Vec<_>>(), &positions, &radii) {
            return; // stale site
        }

        self.add_particle_raw(Particle::new(pos.x, pos.y, type_c, rc));

        let new_idx = self.particles.len() - 1;

        if self.config.relax_steps > 0 {
            // Relaxation handles rate updates and candidate regen
            self.relax_new_particle(new_idx);
        } else {
            self.update_detach_rate(new_idx);
            let nbs = self.get_neighbors(new_idx);
            for nb in nbs {
                self.update_detach_rate(nb);
            }
            let regen_r = (self.config.max_cutoff() + rc) * 2.0;
            self.regen_candidates_near(pos, regen_r);
        }
    }

    // ── Detach ───────────────────────────────────────────────────────────────

    fn detach(&mut self, slot: usize) {
        if slot >= self.particles.len() {
            return;
        }

        let pos = self.particles[slot].pos;
        let bonded = self.get_neighbors(slot);

        // Remove from spatial hash
        self.spatial.remove(slot, pos.x, pos.y);
        self.detach_rates.swap_remove_rate(slot);

        let last = self.particles.len() - 1;
        self.particles.swap_remove(slot);
        if slot != last {
            // The element formerly at `last` is now at `slot`
            let moved = &self.particles[slot];
            self.spatial.remove(last, moved.pos.x, moved.pos.y);
            self.spatial.insert(slot, moved.pos.x, moved.pos.y);
        }

        // Old neighbours lost a bond → update their detach rates
        for nb in bonded {
            // If nb was last, it is now at slot
            let actual = if nb == last && slot != last { slot } else { nb };
            if actual < self.particles.len() {
                self.update_detach_rate(actual);
            }
        }

        let regen_r = self.config.max_cutoff() * 2.0;
        self.regen_candidates_near(pos, regen_r);
    }

    // ── Internal helpers ─────────────────────────────────────────────────────

    /// Push a particle into the data structures without updating rates/candidates.
    fn add_particle_raw(&mut self, p: Particle) {
        let idx = self.particles.len();
        self.spatial.insert(idx, p.pos.x, p.pos.y);
        self.particles.push(p);
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
        // query spatial hash for nearby particles, then filter by actual bond distance
        let p = &self.particles[slot];
        let cutoff = p.radius + self.config.max_radius() + self.config.delta;
        let nearby = self.spatial.query(p.pos.x, p.pos.y, cutoff);
        nearby.into_iter()
            .filter(|&ni| ni != slot && ni < self.particles.len() && p.bonds_to(&self.particles[ni], self.config.delta))
            .collect()
    }

    /// Σ ε(type_i, type_j) over all bonded neighbours j.
    fn binding_energy(&self, slot: usize) -> f64 {
        let type_i = self.particles[slot].type_id;
        self.get_neighbors(slot)
            .iter()
            .map(|&nb| self.config.epsilon(type_i, self.particles[nb].type_id))
            .sum()
    }

    // fn site_binding_energy(&self, site: &CandidateSite) -> f64 {
    //     let type_c = site.type_id;
    //     let rc = self.config.radius(type_c);
    //     let pos = site.pos;
    //     let query_r = rc + self.config.max_radius() + self.config.delta;
    //     let nearby = self.spatial.query(pos.x, pos.y, query_r);
    //     nearby
    //         .iter()
    //         .filter(|&&i| {
    //             let d = self.particles[i].pos.distance(pos);
    //             let contact = self.particles[i].radius + rc;
    //             d >= contact - self.config.delta && d <= contact + self.config.delta
    //         })
    //         .map(|&i| self.config.epsilon(type_c, self.particles[i].type_id))
    //         .sum()
    // }

    fn add_candidate(&mut self, site: CandidateSite) {
        let type_c = site.type_id;
        let pos = site.pos;
        self.candidates.push(site);
        let idx = self.candidates.len() - 1;

        self.candidates_spatial.insert(idx, pos.x, pos.y);

        let rate = self.config.attach_rate(type_c); 
        self.attach_rates.add_rate(rate);
    }
        
    /// Remove stale candidates near `center`, then regenerate from scratch for
    /// all particles within `radius` of `center`.
    fn regen_candidates_near(&mut self, center: DVec2, radius: f64) {
        self.regen_candidates_batch(&[(center, radius)]);
    }

    /// Batch version: delete and regenerate candidate sites for multiple
    /// (center, radius) regions at once.  Merges all deletion, focal, and
    /// overlap queries so overlapping regions don't duplicate work.
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
            self.spatial.query_into(center.x, center.y, search_r, &mut self.query_buf);
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
            self.spatial.query_into(center.x, center.y, overlap_r, &mut self.query_buf);
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

            // Precompute bonded neighbor angles for A (used by arc site filtering below)
            let bonded_cutoff = a_r + self.config.max_radius() + self.config.delta;
            self.spatial.query_into(a_pos.x, a_pos.y, bonded_cutoff, &mut self.query_buf);
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

                // ── 4a. Intersection sites: C touches both A and some neighbor B ──
                let pair_search = a_r + rc + self.config.max_radius() + rc + self.config.delta;
                self.spatial.query_into(a_pos.x, a_pos.y, pair_search, &mut self.query_buf);
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
                        if da < (a_r + rc) - self.config.delta
                            || da > a_r + rc + self.config.delta
                        {
                            continue;
                        }
                        if db < (b_r + rc) - self.config.delta
                            || db > b_r + rc + self.config.delta
                        {
                            continue;
                        }

                        if site_has_overlap(site_pos, rc, self.config.delta, &overlap_idx, &overlap_pos, &overlap_rad) {
                            continue;
                        }
                        self.add_candidate(CandidateSite { pos: site_pos, type_id: type_c });
                    }
                }

                // ── 4b. Arc sites: C touches only A (single-contact bonding) ──
                // Generated for all particles (not just isolated ones).
                // Angles near bonded neighbors are excluded since those regions are
                // already covered by intersection sites above.
                let n_ang = self.config.num_isolated_angles;
                let exclusion = TAU / n_ang as f64; // one angular step
                let r_contact = a_r + rc;
                let max_arc = self.config.max_arc_sites_per_type;

                let mut arc_candidates: Vec<DVec2> = Vec::new();
                for k in 0..n_ang {
                    let theta = (k as f64) * TAU / (n_ang as f64);
                    // Skip angles within `exclusion` of any bonded neighbor direction
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
                    // Select max_arc sites spread evenly across the candidate list
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
        for (i, p) in self.particles.iter().enumerate() {
            let base = i * 4;
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
        // Recompute all detach rates
        for slot in 0..self.particles.len() {
            self.update_detach_rate(slot);
        }
    }

    pub fn set_chemical_potential(&mut self, type_id: usize, mu: f64) {
        if type_id < self.config.n_types() {
            self.config.particle_types[type_id].mu = mu;
            // Update attach rates for all candidates of this type
            let new_rate = self.config.attach_rate(type_id);
            for i in 0..self.candidates.len() {
                if self.candidates[i].type_id == type_id {
                    self.attach_rates.set_rate(i, new_rate);
                }
            }
        }
    }

    // ── Post-attachment relaxation ────────────────────────────────────────────

    /// Steepest-descent relaxation of the neighborhood of a newly attached particle.
    ///
    /// Uses an active-set propagation approach:
    ///   1. Start with only the new particle in the active set.
    ///   2. Each Jacobi step computes forces on all active particles from their
    ///      current neighbors (harmonic spring: F = 2k(r − r_contact) r̂ for any
    ///      pair within r_contact + δ).  Static friction is subtracted from the
    ///      force magnitude (clamped at zero) before applying the displacement.
    ///   3. Reaction forces on non-active neighbors are accumulated; if the total
    ///      exceeds the static-friction threshold the neighbor is activated and
    ///      given an initial displacement in the same step.
    ///   4. Spatial hash is updated at the end of each step.
    ///   5. After all steps, detach rates and candidate sites are refreshed for
    ///      the full active set and its neighborhood.
    fn relax_new_particle(&mut self, new_idx: usize) {
        if self.config.relax_steps == 0 {
            return;
        }

        let k = self.config.spring_k;
        let alpha = self.config.relax_alpha;
        let static_friction = self.static_friction;
        let delta = self.config.delta;
        let max_active: usize = 64; // cap to prevent runaway growth

        let n_particles = self.particles.len();
        let mut active: Vec<usize> = vec![new_idx];
        // Grow the reusable mask if the particle count increased; all new slots are false
        if self.in_active.len() < n_particles {
            self.in_active.resize(n_particles, false);
        }
        self.in_active[new_idx] = true;

        // Save initial positions for post-relaxation candidate regen
        let mut initial_positions: Vec<DVec2> = vec![self.particles[new_idx].pos];

        // Reusable buffers – cleared each step instead of reallocated
        let mut old_positions: Vec<DVec2> = Vec::new();
        let mut forces: Vec<DVec2> = Vec::new();
        let mut touched: Vec<(usize, DVec2)> = Vec::new();
        let mut moved: Vec<(usize, DVec2)> = Vec::new();
        let mut agg: Vec<(usize, DVec2)> = Vec::new();

        for _ in 0..self.config.relax_steps {
            let n_active = active.len();

            // Snapshot positions of currently active particles (Jacobi-style)
            old_positions.clear();
            old_positions.extend(active[..n_active].iter().map(|&idx| self.particles[idx].pos));

            // Compute forces on each active particle
            forces.clear();
            touched.clear();

            for &idx in &active[..n_active] {
                let pos = self.particles[idx].pos;
                let ri = self.particles[idx].radius;
                let type_i = self.particles[idx].type_id;

                let query_r = ri + self.config.max_radius() + delta + 1e-6;
                self.spatial.query_into(pos.x, pos.y, query_r, &mut self.query_buf);

                let mut force = DVec2::ZERO;

                for j_pos in 0..self.query_buf.len() {
                    let j = self.query_buf[j_pos];
                    if j == idx || j >= n_particles {
                        continue;
                    }

                    let r_vec = self.particles[j].pos - pos;
                    let r = r_vec.length();
                    if r < 1e-12 {
                        continue;
                    }
                    let r_hat = r_vec / r;
                    let r_contact = ri + self.particles[j].radius;

                    if r <= r_contact + delta {
                        // Unified harmonic spring: F = 2k(r − r_contact) r̂
                        //   attractive when r > r_contact, repulsive when r < r_contact
                        let f = (2.0 * k * (r - r_contact)) * r_hat;
                        force += f;

                        // For repulsive pairs (ε < 0), add a force derived from
                        // U = −ε·((r_contact + δ − r)/δ)², which costs energy when
                        // the pair is in the bonding shell and goes to zero at r_contact + δ.
                        // F_extra = ε · 2(r_contact + δ − r)/δ² · r̂  (negative → repulsive)
                        let eps = self.config.epsilon(type_i, self.particles[j].type_id);
                        let mut f_nb = -f; // reaction for spring
                        if eps < 0.0 {
                            let f_rep = (2.0 * (r_contact + delta - r)) * r_hat;
                            force += f_rep;
                            f_nb -= f_rep; // reaction for repulsion
                        }

                        // Track total reaction force on non-active neighbors
                        if !self.in_active[j] {
                            touched.push((j, f_nb));
                        }
                    }
                }

                // Apply static friction: reduce magnitude, clamp at zero
                let f_mag = force.length();
                if f_mag > static_friction {
                    force *= (f_mag - static_friction) / f_mag;
                } else {
                    force = DVec2::ZERO;
                }

                forces.push(force);
            }

            // Track all particles that moved this step for spatial hash update
            moved.clear();

            // Apply Jacobi moves to all currently active particles
            for i in 0..n_active {
                let idx = active[i];
                moved.push((idx, old_positions[i]));
                self.particles[idx].pos = old_positions[i] + alpha * forces[i];
            }

            // Aggregate reaction forces on each touched neighbor, activate if
            // total force exceeds static friction, and apply initial displacement.
            if active.len() < max_active {
                // Linear-scan aggregation (small N, no hashing overhead)
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
                        self.particles[nb_idx].pos += alpha * adjusted;
                        self.in_active[nb_idx] = true;
                        active.push(nb_idx);
                        moved.push((nb_idx, old_pos));
                    }
                }
            }

            // Update spatial hash for all particles that moved this step
            for &(idx, ref old_pos) in &moved {
                self.spatial.remove(idx, old_pos.x, old_pos.y);
                let new_pos = self.particles[idx].pos;
                self.spatial.insert(idx, new_pos.x, new_pos.y);
            }
        }

        // ── Post-relaxation: refresh rates and candidate sites ────────────
        // Collect every particle whose bonds may have changed (sorted + dedup, no HashSet)
        let mut rate_indices: Vec<usize> = Vec::new();
        for &idx in &active {
            rate_indices.push(idx);
            rate_indices.extend(self.get_neighbors(idx));
        }
        rate_indices.sort_unstable();
        rate_indices.dedup();
        for &idx in &rate_indices {
            self.update_detach_rate(idx);
        }

        // Regenerate candidates near every active particle (current + initial positions)
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

        // Reset the in_active flags for reuse next call
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

