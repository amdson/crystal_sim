use std::f64::consts::TAU;
use glam::DVec2;

use crate::candidates::{CandidateSite, circle_intersections, site_has_overlap};
use crate::config::SimConfig;
use crate::particle::Particle;
use crate::rates::RateCatalog;
use crate::rng::Rng;
use crate::spatial::SpatialHash;

use std::time::{Duration, Instant, SystemTime};

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

    pub config: SimConfig,
    pub time: f64,
    pub rng: Rng,
    /// Flat f32 particle buffer: [x, y, type_id_f32, radius, ...] stride 4
    pub particle_buf: Vec<f32>,
    pub timers : Vec<Duration>, 
    pub usages: Vec<u32>, 
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
        let n_types = config.n_types();
        let cell_size = config.max_cutoff().max(0.1);

        let mut sim = Self {
            particles: Vec::new(), 
            spatial: SpatialHash::new(cell_size), 
            detach_rates: RateCatalog::new(),
            candidates: Vec::new(), 
            candidates_spatial: SpatialHash::new(cell_size),  
            attach_rates: RateCatalog::new(), 

            time: 0.0,
            rng: Rng::new(config.seed),  
            particle_buf: Vec::new(), 
            config, 
            timers: (0..TimerEntry::Count as usize).map(|_i| Duration::new(0, 0)).collect(), 
            usages: vec![0; TimerEntry::Count as usize] 
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
        let regen_r = (self.config.max_cutoff() + rc) * 2.0;
        self.update_detach_rate(new_idx); 

        // Neighbours' binding energy increased → lower detach rate
        let nbs = self.get_neighbors(new_idx);
        for nb in nbs {
            self.update_detach_rate(nb);
        }

        // self.relax_new_particle(new_idx); 
        self.regen_candidates_near(pos, regen_r);
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

    fn site_binding_energy(&self, site: &CandidateSite) -> f64 {
        let type_c = site.type_id;
        let rc = self.config.radius(type_c);
        let pos = site.pos;
        let query_r = rc + self.config.max_radius() + self.config.delta;
        let nearby = self.spatial.query(pos.x, pos.y, query_r);
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
        self.candidates.push(site);
        let idx = self.candidates.len() - 1;

        self.candidates_spatial.insert(idx, pos.x, pos.y);

        let rate = self.config.attach_rate(type_c); 
        self.attach_rates.add_rate(rate);
    }
        
    /// Remove stale candidates near `center`, then regenerate from scratch for
    /// all particles within `radius` of `center`.
    fn regen_candidates_near(&mut self, center: DVec2, radius: f64) {
        let r_sq = radius * radius;
        let eps_dedup = self.config.max_radius() * 0.05;
        let eps_dedup_sq = eps_dedup * eps_dedup;

        let nearby_sites = self.candidates_spatial.query(center.x, center.y, radius + self.config.max_cutoff());
        let to_delete: Vec<usize> = nearby_sites
            .iter()
            .filter(|&&site_idx| {
                site_idx < self.candidates.len()
                    && self.candidates[site_idx].pos.distance_squared(center) < r_sq
            })
            .copied()
            .collect();
        self.del_candidates_batch(to_delete);

        // 2. Find particles whose neighbourhood overlaps the affected region.
        let search_r = radius + self.config.max_cutoff() * 2.0;
        let focal: Vec<usize> = self
            .spatial
            .query(center.x, center.y, search_r)
            .into_iter()
            .filter(|&i| i < self.particles.len())
            .collect();

        // 3. Build overlap-check arrays once from a wider query that covers all
        //    particles that could possibly overlap any candidate site generated
        //    from focal particles.  A site is at most ~2*max_radius from a focal
        //    particle, and an overlapper is at most ~2*max_radius from the site.
        let overlap_r = search_r + 4.0 * self.config.max_radius();
        let overlap_set: Vec<usize> = self
            .spatial
            .query(center.x, center.y, overlap_r)
            .into_iter()
            .filter(|&i| i < self.particles.len())
            .collect();
        let overlap_pos: Vec<DVec2> = overlap_set.iter().map(|&i| self.particles[i].pos).collect();
        let overlap_rad: Vec<f64> = overlap_set.iter().map(|&i| self.particles[i].radius).collect();
        let overlap_idx: Vec<usize> = (0..overlap_set.len()).collect();

        // 4. For each focal particle A and each candidate type C, generate sites.
        for &a_idx in &focal {
            let a_pos = self.particles[a_idx].pos;
            let a_r = self.particles[a_idx].radius;

            for type_c in 0..self.config.n_types() {
                let rc = self.config.radius(type_c);

                // All particles close enough to A that C can touch both
                let pair_search = a_r + rc + self.config.max_radius() + rc + self.config.delta;
                let nearby_indices = self.spatial.query(a_pos.x, a_pos.y, pair_search);

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
                        found_pair = true;
                        
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
                            if site_has_overlap(site_pos, rc, self.config.delta, &overlap_idx, &overlap_pos, &overlap_rad) {
                                continue; 
                            }
                            self.add_candidate(CandidateSite { pos: site_pos, type_id: type_c });
                        }
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
        }
    }

    // ── Post-attachment relaxation ────────────────────────────────────────────

    /// Steepest-descent relaxation of the newly attached particle.
    ///
    /// Uses a smooth harmonic potential:
    ///   - bonded pairs:     U = spring_k * (r - r_contact)^2   (min at exact contact)
    ///   - non-bonded nearby: U = spring_k * max(0, r_contact - r)^2  (soft repulsion only)
    ///
    /// After convergence (or max steps), must call regen_candidates_near to update candidate sites/rates.
    // regen_candidates_near should be rewritten to take a list of changed particles instead of a center/radius, so that we can call it on the affected neighborhood instead of a large fixed radius.
    // modify other code calling regen_candidates_near accordingly
    fn relax_new_particle(&mut self, new_idx: usize) {
        if self.config.relax_steps == 0 {
            return;
        }
        let k = self.config.spring_k;
        let alpha = self.config.relax_alpha;
        let static_friction = self.config.max_radius() * 0.05; // below this force, treat as static friction and don't move
        

        for _ in 0..self.config.relax_steps {
            let pos = self.particles[new_idx].pos;
            let ri = self.particles[new_idx].radius;

            // Collect nearby particles (need snapshot to avoid borrow conflict)
            let query_r = ri + self.config.max_radius() + self.config.delta + 1e-6;
            let nearby = self.spatial.query(pos.x, pos.y, query_r);

            let mut force = DVec2::ZERO;

            for &j in &nearby {
                if j == new_idx {
                    continue;
                }
                let r_vec = self.particles[j].pos - pos; // points from new toward j
                let r = r_vec.length();
                if r < 1e-12 {
                    continue;
                }
                let r_hat = r_vec / r;
                let r_contact = ri + self.particles[j].radius;

                if self.neighbors[new_idx].contains(&j) {
                    // Bonded spring: attractive if stretched, repulsive if compressed
                    let stretch = r - r_contact;
                    force += (2.0 * k * stretch) * r_hat;
                } else {
                    // Non-bonded soft-core: repulsive only
                    let overlap = r_contact - r;
                    if overlap > 0.0 {
                        force -= (2.0 * k * overlap) * r_hat;
                    }
                }
            }

            // Move particle and update spatial hash
            let old_pos = self.particles[new_idx].pos;
            self.particles[new_idx].pos += alpha * force;
            let new_pos = self.particles[new_idx].pos;
            self.spatial.remove(new_idx, old_pos.x, old_pos.y);
            self.spatial.insert(new_idx, new_pos.x, new_pos.y);
        }

        // Refresh rates for new particle and all affected neighbors
        self.update_detach_rate(new_idx);
        let nbs = self.get_neighbors(new_idx);
        for nb in nbs {
            self.update_detach_rate(nb);
        }
        // Also update old neighbors that may have lost the bond
        for nb in old_neighbors {
            if !nbs.contains(&nb) {
                self.update_detach_rate(nb);
            }
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

