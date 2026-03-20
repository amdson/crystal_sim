use std::collections::HashMap;
use std::f32::consts::PI;
use serde::Deserialize;
use glam::Vec2;

#[derive(Deserialize, Clone, Debug)]
pub struct PatchTypeDef {
    pub name: String,
    #[serde(default)]
    pub color: Option<String>,
}

#[derive(Deserialize, Clone, Debug, Default)]
pub struct PatchDef {
    /// Unique string tag identifying the patch type (e.g. "A", "B").
    pub patch_type: String,
    /// Angular position of the patch relative to the particle body frame, in degrees.
    pub position_deg: f32,
    /// Computed from `position_deg` during `init_cache`. Not serialised.
    #[serde(skip)]
    pub position_rad: f32,
    /// (cos(position_rad), sin(position_rad)) — precomputed for trig-free patch_dir.
    #[serde(skip)]
    pub position_cs: Vec2,
    /// Interned patch type id for hot-path lookup. Computed during `init_cache`.
    #[serde(skip)]
    pub patch_type_id: usize,
}

#[derive(Deserialize, Clone, Debug)]
pub struct PatchInteraction {
    /// [patch_type_0, patch_type_1] - order determines which sigma applies to which side.
    pub types: [String; 2],
    /// Patch interaction strength. Negative = attractive, positive = repulsive.
    /// Matches the scalar epsilon convention: energy = epsilon * bump * g.
    pub epsilon: f32,
    /// Gaussian half-width in degrees for each side: [sigma_0_deg, sigma_1_deg].
    pub angular_width_deg: [f32; 2],
    /// LJ cutoff expressed as a multiple of r_contact.
    pub cutoff: f32,
    /// Computed from `angular_width_deg` during `init_cache`. Not serialised.
    #[serde(skip)]
    pub angular_width: [f32; 2],
}

#[derive(Clone, Copy, Debug, Default)]
pub struct PatchLutEntry {
    pub epsilon: f32,
    pub sigma_i: f32,
    pub sigma_j: f32,
    pub cutoff: f32,
    pub enabled: bool,
}

#[derive(Deserialize, Clone, Debug)]
pub struct ParticleTypeDef {
    /// Particle radius sigma_i
    pub radius: f32,
    /// CSS hex color for the renderer, e.g. "#4a90d9"
    pub color: String,
    /// Chemical potential mu_t controlling supersaturation
    pub mu: f32,
    /// Directional patches on this particle type.
    #[serde(default)]
    pub patches: Vec<PatchDef>,
}

#[derive(Deserialize, Clone, Debug, Default)]
pub struct InitialParticleDef {
    pub x: f32,
    pub y: f32,
    pub type_id: usize,
    #[serde(default)]
    pub orientation_deg: f32,
    #[serde(default)]
    pub frozen: bool,
}

/// Full simulation configuration passed from JS as a JSON string.
#[derive(Deserialize, Clone, Debug)]
pub struct SimConfig {
    pub particle_types: Vec<ParticleTypeDef>,
    /// NxN symmetric interaction matrix epsilon_ij. Positive = attractive.
    /// Used as fallback when no patch interactions are defined.
    pub epsilon: Vec<Vec<f32>>,
    /// Bonding shell half-width delta: particles bond when |d - contact| <= delta,
    /// and only overlap when d < contact - delta
    pub delta: f32,
    /// Temperature T (k=1, so kT = temperature)
    pub temperature: f32,
    /// Attempt frequency nu
    pub nu: f32,
    /// RNG seed
    pub seed: u64,
    /// How many arc angles to consider around each particle (all are evaluated, best are selected)
    #[serde(default = "default_angles")]
    pub num_isolated_angles: u32,
    /// Max arc sites submitted per (particle, candidate-type) pair after angular filtering
    #[serde(default = "default_max_arc_sites")]
    pub max_arc_sites_per_type: usize,
    /// Steepest-descent steps run after each attachment (0 = disabled)
    #[serde(default = "default_relax_steps")]
    pub relax_steps: u32,
    /// Step size for the steepest-descent relaxation
    #[serde(default = "default_relax_alpha")]
    pub relax_alpha: f32,
    /// Harmonic spring stiffness for the relaxation potential
    #[serde(default = "default_spring_k")]
    pub spring_k: f32,
    /// LJ cutoff as a multiple of r_contact (e.g. 2.0 -> cutoff at 2x contact distance)
    #[serde(default = "default_lj_cutoff_factor")]
    pub lj_cutoff_factor: f32,
    #[serde(default = "default_static_friction")]
    pub static_friction: f32,
    #[serde(default = "default_steps_per_frame")]
    pub steps_per_frame: u32,
    #[serde(default = "default_relax_damping")]
    pub relax_damping: f32,
    /// Optional patch type color definitions: [{name: "A", color: "#e84040"}, ...]
    #[serde(default)]
    pub patch_types: Vec<PatchTypeDef>,
    /// Patch-pair interaction definitions.
    #[serde(default)]
    pub patch_interactions: Vec<PatchInteraction>,
    /// Deterministic test mode: no KMC events, only force/torque integration.
    #[serde(default)]
    pub testing_mode: bool,
    /// Print per-particle diagnostics (energy / force / torque) each test step.
    #[serde(default = "default_print_test_diagnostics")]
    pub print_test_diagnostics: bool,
    /// Explicit initial particle placements for testing mode.
    #[serde(default)]
    pub initial_particles: Vec<InitialParticleDef>,

    // Cached / computed fields (not serialised)
    /// Cached largest particle radius (computed after deserialization)
    #[serde(skip)]
    cached_max_radius: f32,
    /// Cached largest bonding cutoff (computed after deserialization)
    #[serde(skip)]
    cached_max_cutoff: f32,
    /// Number of distinct patch types after interning.
    #[serde(skip)]
    pub patch_type_count: usize,
    /// Dense directed LUT: [i * patch_type_count + j] gives parameters for
    /// interactions where patch type i is on the i-side and type j on the j-side.
    #[serde(skip)]
    pub patch_lut: Vec<PatchLutEntry>,
    /// Patch colors indexed by intern id (same order as patch_type_ids).
    #[serde(skip)]
    pub patch_colors: Vec<Option<String>>,
}

fn default_angles() -> u32 { 16 }
fn default_max_arc_sites() -> usize { 6 }
fn default_steps_per_frame() -> u32 { 150 }
fn default_relax_steps() -> u32 { 30 }
fn default_relax_alpha() -> f32 { 0.01 }
fn default_spring_k() -> f32 { 50.0 }
fn default_lj_cutoff_factor() -> f32 { 2.0 }
fn default_static_friction() -> f32 { 0.01 }
fn default_relax_damping() -> f32 { 0.8 }
fn default_print_test_diagnostics() -> bool { true }

impl SimConfig {
    /// Recompute cached derived values. Must be called after deserialization
    /// or after mutating radii/delta.
    pub fn init_cache(&mut self) {
        let mut patch_type_ids: HashMap<String, usize> = HashMap::new();
        let mut intern_patch_type = |name: &str| -> usize {
            if let Some(&id) = patch_type_ids.get(name) {
                id
            } else {
                let id = patch_type_ids.len();
                patch_type_ids.insert(name.to_string(), id);
                id
            }
        };

        // Intern all patch type names seen on particles.
        for t in &self.particle_types {
            for p in &t.patches {
                intern_patch_type(&p.patch_type);
            }
        }
        // Also intern names that appear only in interaction table.
        for int in &self.patch_interactions {
            intern_patch_type(&int.types[0]);
            intern_patch_type(&int.types[1]);
        }

        // Patch positions in radians and numeric patch ids.
        for t in &mut self.particle_types {
            for p in &mut t.patches {
                p.position_rad = p.position_deg * PI / 180.0;
                p.position_cs = Vec2::new(p.position_rad.cos(), p.position_rad.sin());
                p.patch_type_id = *patch_type_ids
                    .get(&p.patch_type)
                    .expect("interned patch type missing");
            }
        }

        // Patch interaction angular widths in radians.
        for int in &mut self.patch_interactions {
            int.angular_width[0] = int.angular_width_deg[0] * PI / 180.0;
            int.angular_width[1] = int.angular_width_deg[1] * PI / 180.0;
        }

        // Build dense directed interaction LUT for hot-path lookup.
        self.patch_type_count = patch_type_ids.len();
        let lut_len = self.patch_type_count.saturating_mul(self.patch_type_count);
        self.patch_lut.clear();
        self.patch_lut.resize(lut_len, PatchLutEntry::default());

        for int in &self.patch_interactions {
            let a = *patch_type_ids
                .get(&int.types[0])
                .expect("interned patch type missing");
            let b = *patch_type_ids
                .get(&int.types[1])
                .expect("interned patch type missing");
            let n = self.patch_type_count;

            let ab = a * n + b;
            self.patch_lut[ab] = PatchLutEntry {
                epsilon: int.epsilon,
                sigma_i: int.angular_width[0],
                sigma_j: int.angular_width[1],
                cutoff: int.cutoff,
                enabled: true,
            };

            if a != b {
                let ba = b * n + a;
                self.patch_lut[ba] = PatchLutEntry {
                    epsilon: int.epsilon,
                    sigma_i: int.angular_width[1],
                    sigma_j: int.angular_width[0],
                    cutoff: int.cutoff,
                    enabled: true,
                };
            }
        }

        // Build patch color table indexed by intern id.
        self.patch_colors.clear();
        self.patch_colors.resize(self.patch_type_count, None);
        for ptype in &self.patch_types {
            if let Some(&id) = patch_type_ids.get(&ptype.name) {
                self.patch_colors[id] = ptype.color.clone();
            }
        }

        self.cached_max_radius = self.particle_types.iter().map(|t| t.radius).fold(0.0f32, f32::max);
        let n = self.n_types();
        let mut max_c = 0.0f32;
        for i in 0..n {
            for j in 0..n {
                let c = self.cutoff(i, j);
                if c > max_c { max_c = c; }
            }
        }
        self.cached_max_cutoff = max_c;
    }

    pub fn n_types(&self) -> usize {
        self.particle_types.len()
    }

    pub fn radius(&self, type_id: usize) -> f32 {
        self.particle_types[type_id].radius
    }

    pub fn epsilon(&self, i: usize, j: usize) -> f32 {
        self.epsilon[i][j]
    }

    /// Hard-core contact distance between types i and j
    pub fn contact(&self, i: usize, j: usize) -> f32 {
        self.radius(i) + self.radius(j)
    }

    /// Upper bound of bonding shell for types i and j
    pub fn cutoff(&self, i: usize, j: usize) -> f32 {
        self.contact(i, j) + self.delta
    }

    pub fn max_radius(&self) -> f32 {
        self.cached_max_radius
    }

    /// Largest possible bonding cutoff across all type pairs
    pub fn max_cutoff(&self) -> f32 {
        self.cached_max_cutoff
    }

    /// Attachment rate per candidate site for type t
    pub fn attach_rate(&self, type_id: usize, e_bind: f32) -> f64 {
        (self.nu as f64) * (((self.particle_types[type_id].mu as f64) + (e_bind as f64 / 2.0)) / self.temperature as f64).exp()
    }

    /// Detachment rate given a binding energy E_bind
    pub fn detach_rate(&self, e_bind: f32) -> f64 {
        (self.nu as f64) * ((-e_bind as f64 / 2.0) / self.temperature as f64).exp()
    }

    /// Returns true if patchy interactions are configured for any particle type.
    pub fn has_patches(&self) -> bool {
        !self.patch_interactions.is_empty()
            && self.particle_types.iter().any(|t| !t.patches.is_empty())
    }
}
