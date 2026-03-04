use glam::DVec2;

#[derive(Clone, Debug)]
pub struct Particle {
    pub pos: DVec2,
    pub type_id: usize,
    pub radius: f64,
}

impl Particle {
    pub fn new(x: f64, y: f64, type_id: usize, radius: f64) -> Self {
        Self { pos: DVec2::new(x, y), type_id, radius }
    }

    pub fn dist(&self, other: &Particle) -> f64 {
        self.pos.distance(other.pos)
    }

    pub fn dist_sq(&self, other: &Particle) -> f64 {
        self.pos.distance_squared(other.pos)
    }

    /// Hard-core overlap: centers are closer than σ_i + σ_j
    pub fn overlaps(&self, other: &Particle) -> bool {
        let contact = self.radius + other.radius;
        self.dist_sq(other) < contact * contact * (1.0 - 1e-9)
    }

    /// Bond: distance is within [contact, contact + δ)
    pub fn bonds_to(&self, other: &Particle, delta: f64) -> bool {
        let r = self.dist(other);
        let contact = self.radius + other.radius;
        r >= contact * (1.0 - 1e-9) && r < contact + delta
    }
}
