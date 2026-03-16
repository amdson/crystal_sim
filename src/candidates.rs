use glam::Vec2;

/// A potential attachment location for a particle of `type_id`.
#[derive(Clone, Debug)]
pub struct CandidateSite {
    pub pos: Vec2,
    pub type_id: usize,
    pub orientation: Vec2,  // unit vector
}

/// Compute the (up to 2) intersection points of two circles.
///
/// Circle 1: center `a`, radius `r1`
/// Circle 2: center `b`, radius `r2`
pub fn circle_intersections(a: Vec2, r1: f32, b: Vec2, r2: f32) -> [Option<Vec2>; 2] {
    let d_vec = b - a;
    let d = d_vec.length();

    if d < 1e-12 {
        return [None, None]; // coincident centres
    }
    if d > r1 + r2 + 1e-9 {
        return [None, None]; // circles too far apart
    }
    if d < (r1 - r2).abs() - 1e-9 {
        return [None, None]; // one circle inside the other
    }

    let a_coeff = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
    let h_sq = r1 * r1 - a_coeff * a_coeff;
    if h_sq < 0.0 {
        return [None, None];
    }

    let mid = a + d_vec * (a_coeff / d);

    if h_sq < 1e-18 {
        return [Some(mid), None]; // tangent — single point
    }

    let h = h_sq.sqrt();
    let perp = Vec2::new(-d_vec.y, d_vec.x) * (h / d);
    [Some(mid + perp), Some(mid - perp)]
}

/// True if a proposed site at `pos` with radius `rc` hard-core overlaps any
/// existing particle (center-to-center < σ_i + σ_C - δ).
///
/// `nearby_indices` should come from a spatial hash query with an appropriate radius;
/// `particles_pos` and `particles_radius` are parallel slices indexed by those indices.
pub fn site_has_overlap(
    pos: Vec2,
    rc: f32,
    delta: f32,
    nearby: &[usize],
    positions: &[Vec2],
    radii: &[f32],
) -> bool {
    for &idx in nearby {
        let contact = rc + radii[idx];
        let overlap_dist = contact - delta;
        if overlap_dist > 0.0 && pos.distance_squared(positions[idx]) < overlap_dist * overlap_dist {
            return true;
        }
    }
    false
}
