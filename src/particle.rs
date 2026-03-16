use glam::Vec2;

#[derive(Clone, Debug)]
pub struct Particle {
    pub pos: Vec2,
    pub type_id: usize,
    pub radius: f32,
    pub orientation: Vec2,  // unit vector; angle = atan2(y, x)
}

impl Particle {
    pub fn new(x: f32, y: f32, type_id: usize, radius: f32, orientation: Vec2) -> Self {
        Self { pos: Vec2::new(x, y), type_id, radius, orientation }
    }

    pub fn dist(&self, other: &Particle) -> f32 {
        self.pos.distance(other.pos)
    }

    pub fn dist_sq(&self, other: &Particle) -> f32 {
        self.pos.distance_squared(other.pos)
    }

    /// Hard-core overlap: centers are closer than σ_i + σ_j - δ
    pub fn overlaps(&self, other: &Particle, delta: f32) -> bool {
        let contact = self.radius + other.radius;
        let overlap_dist = contact - delta;
        overlap_dist > 0.0 && self.dist_sq(other) < overlap_dist * overlap_dist
    }

    /// Bond: distance is within [contact - δ, contact + δ]
    pub fn bonds_to(&self, other: &Particle, delta: f32) -> bool {
        let r = self.dist(other);
        let contact = self.radius + other.radius;
        r >= contact - delta && r <= contact + delta
    }
}
