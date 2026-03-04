use std::f64::consts::TAU;
use glam::DVec2;

use crate::candidates::{CandidateSite, circle_intersections, site_has_overlap};
use crate::config::SimConfig;
use crate::particle::Particle;
use crate::rates::RateCatalog;
use crate::rng::Rng;
use crate::spatial::SpatialHash;

pub struct Simulation {
    pub particles: Vec<Particle>,
    /// neighbors[i] = indices of particles bonded to particle i
    pub neighbors: Vec<Vec<usize>>,
    pub spatial: SpatialHash,
    /// candidates[type_id] = current valid attachment sites for that type
    pub candidates: Vec<Vec<CandidateSite>>,
    pub rates: RateCatalog,
    pub config: SimConfig,
    pub time: f64,
    pub rng: Rng,
    /// Flat f32 particle buffer: [x, y, type_id_f32, radius, ...] stride 4
    pub particle_buf: Vec<f32>,
}

impl Simulation {
    pub fn new(config: SimConfig) -> Self {
        let n_types = config.n_types();
        let cell_size = config.max_cutoff().max(0.1);

        let mut sim = Self {
            particles: Vec::new(),
            neighbors: Vec::new(),
            spatial: SpatialHash::new(cell_size),
            candidates: vec![Vec::new(); n_types],
            rates: RateCatalog::new(n_types),
            time: 0.0,
            rng: Rng::new(config.seed),
            particle_buf: Vec::new(),
            config,
        };

        // Seed crystal: one particle of each type placed at the origin,
        // offset by contact distance so they bond immediately.
        if sim.config.n_types() >= 1 {
            let r0 = sim.config.radius(0);
            sim.add_particle_raw(Particle::new(0.0, 0.0, 0, r0));
        }
        if sim.config.n_types() >= 2 {
            let r0 = sim.config.radius(0);
            let r1 = sim.config.radius(1);
            let contact = r0 + r1;
            sim.add_particle_raw(Particle::new(contact, 0.0, 1, r1));
        }

        // Build initial candidate sites and rates
        let regen_r = sim.config.max_cutoff() * 4.0;
        sim.regen_candidates_near(DVec2::ZERO, regen_r);
        sim.refresh_attach_rates();

        sim
    }

    // ── KMC step loop ────────────────────────────────────────────────────────

    pub fn step(&mut self, n: u32) {
        for _ in 0..n {
            let r_total = self.rates.total();
            if r_total <= 0.0 {
                break;
            }

            let u1 = self.rng.next_f64();
            let u2 = self.rng.next_f64();

            // BKL time advance: Δt = -ln(u2) / R_total
            self.time += -(u2.ln()) / r_total;

            let event = self.rates.select(u1);
            let n_types = self.config.n_types();

            if event < n_types {
                // Attachment event for type `event`
                let type_c = event;
                if self.candidates[type_c].is_empty() {
                    continue;
                }
                let site_idx = self.rng.next_usize(self.candidates[type_c].len());
                let site = self.candidates[type_c][site_idx].clone();
                self.attach(site);
            } else {
                // Detachment event for particle slot
                let slot = event - n_types;
                if slot < self.particles.len() {
                    self.detach(slot);
                }
            }
        }

        self.rebuild_particle_buf();
    }

    // ── Attach ───────────────────────────────────────────────────────────────

    fn attach(&mut self, site: CandidateSite) {
        let type_c = site.type_id;
        let rc = self.config.radius(type_c);
        let pos = site.pos;

        // Final hard-core check — site may be stale after previous events
        let query_r = rc + self.config.max_radius() + 1e-6;
        let nearby = self.spatial.query(pos.x, pos.y, query_r);
        let positions: Vec<DVec2> = nearby.iter().map(|&i| self.particles[i].pos).collect();
        let radii: Vec<f64> = nearby.iter().map(|&i| self.particles[i].radius).collect();
        if site_has_overlap(pos, rc, &(0..nearby.len()).collect::<Vec<_>>(), &positions, &radii) {
            return; // stale site
        }

        self.add_particle_raw(Particle::new(pos.x, pos.y, type_c, rc));

        let new_idx = self.particles.len() - 1;
        let regen_r = (self.config.max_cutoff() + rc) * 2.0;
        self.update_neighbors_for_new(new_idx);
        self.update_detach_rate(new_idx);

        // Neighbours' binding energy increased → lower detach rate
        let nbs = self.neighbors[new_idx].clone();
        for nb in nbs {
            self.update_detach_rate(nb);
        }

        self.regen_candidates_near(pos, regen_r);
        self.refresh_attach_rates();
    }

    // ── Detach ───────────────────────────────────────────────────────────────

    fn detach(&mut self, slot: usize) {
        if slot >= self.particles.len() {
            return;
        }

        let pos = self.particles[slot].pos;
        let bonded = std::mem::take(&mut self.neighbors[slot]);

        // Remove from spatial hash
        self.spatial.remove(slot, pos.x, pos.y);

        // Remove `slot` from every bonded neighbor's list
        for &nb in &bonded {
            self.neighbors[nb].retain(|&x| x != slot);
        }

        let last = self.particles.len() - 1;

        if slot != last {
            // Move the last particle into `slot`
            let last_pos = self.particles[last].pos;
            self.spatial.remove(last, last_pos.x, last_pos.y);
            self.particles.swap(slot, last);
            self.neighbors.swap(slot, last);
            self.spatial.insert(slot, last_pos.x, last_pos.y);

            // Rewrite all neighbor-list references: last → slot
            let moved_nbs = self.neighbors[slot].clone();
            for &nb in &moved_nbs {
                for x in self.neighbors[nb].iter_mut() {
                    if *x == last {
                        *x = slot;
                    }
                }
            }

            // Swap-remove in rate catalog; moved particle keeps its rate
            self.rates.remove_particle(slot);
            // Recompute: moved particle may have lost the bond to `slot`
            self.update_detach_rate(slot);
        } else {
            self.rates.remove_particle(slot);
        }

        self.particles.pop();
        self.neighbors.pop();

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
        self.refresh_attach_rates();
    }

    // ── Internal helpers ─────────────────────────────────────────────────────

    /// Push a particle into the data structures without updating rates/candidates.
    fn add_particle_raw(&mut self, p: Particle) {
        let idx = self.particles.len();
        self.spatial.insert(idx, p.pos.x, p.pos.y);
        self.particles.push(p);
        self.neighbors.push(Vec::new());
        self.rates.add_particle(0.0); // placeholder; caller must call update_detach_rate
    }

    /// Populate neighbors[new_idx] and update symmetric entries.
    fn update_neighbors_for_new(&mut self, new_idx: usize) {
        let p = self.particles[new_idx].clone();
        let cutoff = p.radius + self.config.max_radius() + self.config.delta;
        let nearby = self.spatial.query(p.pos.x, p.pos.y, cutoff);

        for &ni in &nearby {
            if ni == new_idx || ni >= self.particles.len() {
                continue;
            }
            let other = &self.particles[ni];
            if p.bonds_to(other, self.config.delta) {
                self.neighbors[new_idx].push(ni);
                self.neighbors[ni].push(new_idx);
            }
        }
    }

    /// Compute E_bind for particle at `slot` and refresh its detachment rate.
    pub fn update_detach_rate(&mut self, slot: usize) {
        if slot >= self.particles.len() {
            return;
        }
        let e_bind = self.binding_energy(slot);
        let rate = self.config.detach_rate(e_bind);
        self.rates.set_detach_rate(slot, rate);
    }

    /// Σ ε(type_i, type_j) over all bonded neighbours j.
    fn binding_energy(&self, slot: usize) -> f64 {
        let type_i = self.particles[slot].type_id;
        self.neighbors[slot]
            .iter()
            .map(|&nb| self.config.epsilon(type_i, self.particles[nb].type_id))
            .sum()
    }

    /// Remove stale candidates near `center`, then regenerate from scratch for
    /// all particles within `radius` of `center`.
    fn regen_candidates_near(&mut self, center: DVec2, radius: f64) {
        let r_sq = radius * radius;
        let eps_dedup = self.config.max_radius() * 0.05;
        let eps_dedup_sq = eps_dedup * eps_dedup;

        // 1. Drop any existing candidate site within the affected region.
        for type_c in 0..self.config.n_types() {
            self.candidates[type_c].retain(|s| s.pos.distance_squared(center) >= r_sq);
        }

        // 2. Find particles whose neighbourhood overlaps the affected region.
        let search_r = radius + self.config.max_cutoff();
        let focal: Vec<usize> = self
            .spatial
            .query(center.x, center.y, search_r)
            .into_iter()
            .filter(|&i| i < self.particles.len())
            .collect();

        // 3. For each focal particle A and each candidate type C, generate sites.
        for &a_idx in &focal {
            let a_pos = self.particles[a_idx].pos;
            let a_r = self.particles[a_idx].radius;

            for type_c in 0..self.config.n_types() {
                let rc = self.config.radius(type_c);

                // All particles close enough to A that C can touch both
                let pair_search = a_r + rc + self.config.max_radius() + rc + self.config.delta;
                let nearby_indices = self.spatial.query(a_pos.x, a_pos.y, pair_search);

                // Pre-build position/radius slices for overlap check
                let all_pos: Vec<DVec2> = self.particles.iter().map(|p| p.pos).collect();
                let all_rad: Vec<f64> = self.particles.iter().map(|p| p.radius).collect();
                let all_idx: Vec<usize> = (0..self.particles.len()).collect();

                let mut found_pair = false;

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

                        // Must be within bonding shell of both A and B
                        let da = a_pos.distance(site_pos);
                        let db = b_pos.distance(site_pos);
                        if da < (a_r + rc) * (1.0 - 1e-6)
                            || da >= a_r + rc + self.config.delta
                        {
                            continue;
                        }
                        if db < (b_r + rc) * (1.0 - 1e-6)
                            || db >= b_r + rc + self.config.delta
                        {
                            continue;
                        }

                        if site_has_overlap(site_pos, rc, &all_idx, &all_pos, &all_rad) {
                            continue;
                        }

                        // Deduplicate against already-added sites for this type
                        let dup = self.candidates[type_c]
                            .iter()
                            .any(|s| s.pos.distance_squared(site_pos) < eps_dedup_sq);
                        if !dup {
                            self.candidates[type_c]
                                .push(CandidateSite { pos: site_pos, type_id: type_c });
                            found_pair = true;
                        }
                    }
                }

                // Isolated particle: generate arc sites if no bonded neighbours at all
                if !found_pair {
                    let bonded_cutoff = a_r + self.config.max_radius() + self.config.delta;
                    let has_bonded_neighbor = self
                        .spatial
                        .query(a_pos.x, a_pos.y, bonded_cutoff)
                        .iter()
                        .any(|&ni| {
                            if ni == a_idx || ni >= self.particles.len() {
                                return false;
                            }
                            let other = &self.particles[ni];
                            let d = a_pos.distance(other.pos);
                            d < a_r + other.radius + self.config.delta
                        });

                    if !has_bonded_neighbor {
                        let n_ang = self.config.num_isolated_angles;
                        let r_contact = a_r + rc;
                        for k in 0..n_ang {
                            let theta = (k as f64) * TAU / (n_ang as f64);
                            let site_pos =
                                a_pos + DVec2::new(theta.cos(), theta.sin()) * r_contact;

                            if site_has_overlap(site_pos, rc, &all_idx, &all_pos, &all_rad) {
                                continue;
                            }
                            let dup = self.candidates[type_c]
                                .iter()
                                .any(|s| s.pos.distance_squared(site_pos) < eps_dedup_sq);
                            if !dup {
                                self.candidates[type_c]
                                    .push(CandidateSite { pos: site_pos, type_id: type_c });
                            }
                        }
                    }
                }
            }
        }
    }

    /// Recompute aggregate attachment rates from current candidate counts.
    fn refresh_attach_rates(&mut self) {
        for type_c in 0..self.config.n_types() {
            let n_sites = self.candidates[type_c].len() as f64;
            let rate_per_site = self.config.attach_rate(type_c);
            self.rates.set_attach_rate(type_c, n_sites * rate_per_site);
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
        self.refresh_attach_rates();
    }

    pub fn set_chemical_potential(&mut self, type_id: usize, mu: f64) {
        if type_id < self.config.n_types() {
            self.config.particle_types[type_id].mu = mu;
            self.refresh_attach_rates();
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
