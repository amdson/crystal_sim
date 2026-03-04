use serde::Deserialize;

#[derive(Deserialize, Clone, Debug)]
pub struct ParticleTypeDef {
    /// Particle radius σ_i
    pub radius: f64,
    /// CSS hex color for the renderer, e.g. "#4a90d9"
    pub color: String,
    /// Chemical potential μ_t controlling supersaturation
    pub mu: f64,
}

/// Full simulation configuration passed from JS as a JSON string.
#[derive(Deserialize, Clone, Debug)]
pub struct SimConfig {
    pub particle_types: Vec<ParticleTypeDef>,
    /// N×N symmetric interaction matrix ε_ij. Positive = attractive.
    pub epsilon: Vec<Vec<f64>>,
    /// Bonding shell width δ: bond forms when contact ≤ r < contact + δ
    pub delta: f64,
    /// Temperature T (k=1, so kT = temperature)
    pub temperature: f64,
    /// Attempt frequency ν
    pub nu: f64,
    /// RNG seed
    pub seed: u64,
    /// How many arc angles to sample around isolated particles
    #[serde(default = "default_angles")]
    pub num_isolated_angles: u32,
    /// Steepest-descent steps run after each attachment (0 = disabled)
    #[serde(default = "default_relax_steps")]
    pub relax_steps: u32,
    /// Step size for the steepest-descent relaxation
    #[serde(default = "default_relax_alpha")]
    pub relax_alpha: f64,
    /// Harmonic spring stiffness for the relaxation potential
    #[serde(default = "default_spring_k")]
    pub spring_k: f64,
}

fn default_angles() -> u32 { 12 }
fn default_relax_steps() -> u32 { 30 }
fn default_relax_alpha() -> f64 { 0.01 }
fn default_spring_k() -> f64 { 50.0 }

impl SimConfig {
    pub fn n_types(&self) -> usize {
        self.particle_types.len()
    }

    pub fn radius(&self, type_id: usize) -> f64 {
        self.particle_types[type_id].radius
    }

    pub fn epsilon(&self, i: usize, j: usize) -> f64 {
        self.epsilon[i][j]
    }

    /// Hard-core contact distance between types i and j
    pub fn contact(&self, i: usize, j: usize) -> f64 {
        self.radius(i) + self.radius(j)
    }

    /// Upper bound of bonding shell for types i and j
    pub fn cutoff(&self, i: usize, j: usize) -> f64 {
        self.contact(i, j) + self.delta
    }

    pub fn max_radius(&self) -> f64 {
        self.particle_types.iter().map(|t| t.radius).fold(0.0f64, f64::max)
    }

    /// Largest possible bonding cutoff across all type pairs
    pub fn max_cutoff(&self) -> f64 {
        let n = self.n_types();
        let mut max = 0.0f64;
        for i in 0..n {
            for j in 0..n {
                let c = self.cutoff(i, j);
                if c > max { max = c; }
            }
        }
        max
    }

    /// Attachment rate per candidate site for type t
    pub fn attach_rate(&self, type_id: usize) -> f64 {
        self.nu * (self.particle_types[type_id].mu / self.temperature).exp()
    }

    /// Detachment rate given a binding energy E_bind
    pub fn detach_rate(&self, e_bind: f64) -> f64 {
        self.nu * ((-e_bind) / self.temperature).exp()
    }
}
