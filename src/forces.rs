//! Force, potential, and torque functions for particle relaxation.
//!
//! Sign convention:
//! Most low-level LJ helpers return energy gradients (nabla_i U). Negate to get
//! physical force on i. Patchy helpers return physical force directly.

use glam::Vec2;

use crate::config::{PatchDef, PatchLutEntry};

/// Rotate `position_cs` by `orientation` — trig-free patch direction.
/// Both arguments are unit vectors (orientation of particle, precomputed patch direction).
/// Equivalent to: angle = orientation_angle + position_rad; Vec2::new(cos, sin).
#[inline]
pub fn patch_dir(orientation: Vec2, position_cs: Vec2) -> Vec2 {
    Vec2::new(
        orientation.x * position_cs.x - orientation.y * position_cs.y,
        orientation.x * position_cs.y + orientation.y * position_cs.x,
    )
}

/// 2^(-1/6): scale factor so that the LJ minimum is exactly at r_contact.
const SIGMA_FACTOR: f32 = 0.890_898_7;

/// Lennard-Jones potential (no cutoff).
pub fn lj_potential(r: f32, r_contact: f32, epsilon: f32) -> f32 {
    let sigma = r_contact * SIGMA_FACTOR;
    let sr = sigma / r;
    let sr6 = sr * sr * sr * sr * sr * sr;
    let sr12 = sr6 * sr6;
    4.0 * epsilon * (sr12 - sr6)
}

/// LJ gradient nabla_i U. Negate for force on i.
/// Legacy behavior: for epsilon > 0 it uses epsilon=1; for epsilon <= 0 it uses pure repulsion.
pub fn lj_force_vec(r_vec: Vec2, r_contact: f32, epsilon: f32, r_cut: f32) -> Vec2 {
    let r2 = r_vec.length_squared();
    if r2 >= r_cut * r_cut || r2 < 1e-24 {
        return Vec2::ZERO;
    }
    let sigma = r_contact * SIGMA_FACTOR;
    let sigma2 = sigma * sigma;
    let sr2 = sigma2 / r2;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    if epsilon > 0.0 {
        let coeff = 24.0 / r2 * (2.0 * sr12 - sr6);
        coeff * r_vec
    } else {
        let coeff = 24.0 * 2.0 / r2 * sr12;
        coeff * r_vec
    }
}

/// Full LJ gradient nabla_i U with explicit epsilon. Negate for force on i.
pub fn lj_force_vec_full(r_vec: Vec2, r_contact: f32, epsilon: f32, r_cut: f32) -> Vec2 {
    let r2 = r_vec.length_squared();
    if r2 >= r_cut * r_cut || r2 < 1e-24 {
        return Vec2::ZERO;
    }
    let sigma = r_contact * SIGMA_FACTOR;
    let sigma2 = sigma * sigma;
    let sr2 = sigma2 / r2;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    let coeff = 24.0 * epsilon / r2 * (2.0 * sr12 - sr6);
    coeff * r_vec
}

/// Purely repulsive r^-12 gradient. Negate for force on i.
pub fn repulsive_force_vec(r_vec: Vec2, r_contact: f32, r_cut: f32) -> Vec2 {
    let r2 = r_vec.length_squared();
    if r2 >= r_cut * r_cut || r2 < 1e-24 {
        return Vec2::ZERO;
    }
    let sigma = r_contact * SIGMA_FACTOR;
    let sigma2 = sigma * sigma;
    let sr2 = sigma2 / r2;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    let coeff = 24.0 * 2.0 / r2 * sr12;
    coeff * r_vec
}

/// Harmonic spring force on i due to j (actual force, not gradient).
pub fn spring_force_vec(r_vec: Vec2, r_contact: f32, k: f32, shell: f32) -> Vec2 {
    let r = r_vec.length();
    if r < 1e-12 || r > r_contact + shell {
        return Vec2::ZERO;
    }
    let r_hat = r_vec / r;
    (2.0 * k * (r - r_contact)) * r_hat
}

/// Total interaction energy for one pair with optional patch modulation.
pub fn patchy_pair_energy(
    r_vec: Vec2,
    r_contact: f32,
    i_patches: &[PatchDef],
    i_ori: Vec2,
    j_patches: &[PatchDef],
    j_ori: Vec2,
    patch_lut: &[PatchLutEntry],
    patch_type_count: usize,
    epsilon_fallback: f32,
    delta: f32,
) -> f32 {
    let r2 = r_vec.length_squared();
    if r2 < 1e-24 {
        return 0.0;
    }
    let r = r2.sqrt();
    let r_hat = r_vec / r; // Unit vector from i to j.

    // Scalar fallback if either side has no patches.
    if i_patches.is_empty() || j_patches.is_empty() {
        if r < r_contact - delta || r > r_contact + delta {
            return 0.0;
        }
        return epsilon_fallback;
    }

    let mut energy = 0.0;
    for i_patch in i_patches {
        let n_a = patch_dir(i_ori, i_patch.position_cs);
        let cos_ia = r_hat.dot(n_a);   // cos of angle from r_hat to patch n_a
        let i_type = i_patch.patch_type_id;

        for j_patch in j_patches {
            let j_type = j_patch.patch_type_id;
            if i_type >= patch_type_count || j_type >= patch_type_count {
                continue;
            }
            let int = &patch_lut[i_type * patch_type_count + j_type];
            if !int.enabled {
                continue;
            }

            let r_cut = r_contact * int.cutoff;
            if r >= r_cut {
                continue;
            }
            let radial = contact_bump(r, r_contact, r_cut);
            if radial <= 0.0 {
                continue;
            }

            let n_b = patch_dir(j_ori, j_patch.position_cs);
            let cos_jb = -r_hat.dot(n_b);  // cos of angle from -r_hat to patch n_b
            // g = exp(-(1-cos)/sigma^2); matches exp(-theta^2/(2*sigma^2)) for small theta.
            let si2 = int.sigma_i * int.sigma_i;
            let sj2 = int.sigma_j * int.sigma_j;
            let g = ((cos_ia - 1.0) / si2).exp() * ((cos_jb - 1.0) / sj2).exp();
            energy += int.epsilon * radial * g;
        }
    }

    energy
}

/// Physical force on i and torques on both particles for a patchy pair.
/// Returns (F_i, tau_i, tau_j).
pub fn patchy_force_torque(
    r_vec: Vec2,
    r_contact: f32,
    i_patches: &[PatchDef],
    i_ori: Vec2,
    j_patches: &[PatchDef],
    j_ori: Vec2,
    patch_lut: &[PatchLutEntry],
    patch_type_count: usize,
    epsilon_fallback: f32,
    lj_cutoff_factor: f32,
) -> (Vec2, f32, f32) {
    let r2 = r_vec.length_squared();
    if r2 < 1e-24 {
        return (Vec2::ZERO, 0.0, 0.0);
    }
    let r = r2.sqrt();
    let r_hat = r_vec / r;
    let r_cut_default = r_contact * lj_cutoff_factor;

    // Scalar fallback if either side has no patches.
    if i_patches.is_empty() || j_patches.is_empty() {
        return (
            -lj_force_vec_full(r_vec, r_contact, epsilon_fallback, r_cut_default),
            0.0,
            0.0,
        );
    }

    // Isotropic backbone repulsion.
    let mut force = -repulsive_force_vec(r_vec, r_contact, r_cut_default);
    let mut tau_i = 0.0f32;
    let mut tau_j = 0.0f32;

    for i_patch in i_patches {
        let n_a = patch_dir(i_ori, i_patch.position_cs);
        let cos_ia = r_hat.dot(n_a);                      // cos(theta_ia)
        let sin_ia = r_hat.x * n_a.y - r_hat.y * n_a.x;  // sin(theta_ia): r_hat cross n_a
        let i_type = i_patch.patch_type_id;

        for j_patch in j_patches {
            let j_type = j_patch.patch_type_id;
            if i_type >= patch_type_count || j_type >= patch_type_count {
                continue;
            }
            let int = &patch_lut[i_type * patch_type_count + j_type];
            if !int.enabled {
                continue;
            }

            let r_cut = r_contact * int.cutoff;
            if r >= r_cut {
                continue;
            }

            let (radial, d_radial_dr) = contact_bump_with_derivative(r, r_contact, r_cut);
            if radial <= 0.0 {
                continue;
            }

            let n_b   = patch_dir(j_ori, j_patch.position_cs);
            let cos_jb =  -r_hat.dot(n_b);                    // cos(theta_jb): (-r_hat)·n_b
            let sin_jb = -(r_hat.x * n_b.y - r_hat.y * n_b.x); // sin(theta_jb): (-r_hat) cross n_b

            let si2 = int.sigma_i * int.sigma_i;
            let sj2 = int.sigma_j * int.sigma_j;
            // g = exp(-(1-cos)/sigma^2); matches exp(-theta^2/(2*sigma^2)) for small theta.
            let g = ((cos_ia - 1.0) / si2).exp() * ((cos_jb - 1.0) / sj2).exp();
            let eps = int.epsilon;

            // Radial term: F_i = -dU/dr_i; dr/dr_i = -r_hat.
            force += eps * d_radial_dr * g * r_hat;

            // Angular correction to linear force.
            // dg/dr_i = -g * n_a_perp/(si2*r) + g * n_b_perp_rhat/(sj2*r)
            // where n_x_perp = n_x - r_hat*(r_hat·n_x).
            let n_a_perp = n_a - r_hat * cos_ia;             // reuse cos_ia
            let n_b_perp = n_b + r_hat * cos_jb;             // = n_b - r_hat*(r_hat·n_b), reuse cos_jb
            force += eps * radial * g / r * (n_a_perp / si2 - n_b_perp / sj2);

            // Torques: tau = -dU/dtheta = eps*radial*g*sin/sigma^2.
            // sin comes from d(cos)/d(theta) = -sin, and the negation from -dU.
            tau_i += eps * radial * g * sin_ia / si2;
            tau_j += eps * radial * g * sin_jb / sj2;
        }
    }
    (force, tau_i, tau_j)
}

/// Compact, strictly positive bump centered at contact.
/// Support is [2*r_contact - r_cut, r_cut], and zero outside.
/// Normalized to 1 at r = r_contact.
#[inline]
fn contact_bump(r: f32, r_contact: f32, r_cut: f32) -> f32 {
    contact_bump_with_derivative(r, r_contact, r_cut).0
}

/// Compact bump and radial derivative used for patch KMC energy and relaxation.
/// `d_bump_dr` is derivative w.r.t. scalar distance `r`.
#[inline]
fn contact_bump_with_derivative(r: f32, r_contact: f32, r_cut: f32) -> (f32, f32) {
    let width = r_cut - r_contact;
    if width <= 1e-12 {
        return (0.0, 0.0);
    }
    let x = (r - r_contact) / width;
    let x2 = x * x;
    if x2 >= 1.0 {
        return (0.0, 0.0);
    }
    let inv = 1.0 / (1.0 - x2);
    let bump = (1.0 - inv).exp();
    let d_bump_dr = bump * (-(2.0 * x) * inv * inv) / width;
    (bump, d_bump_dr)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    fn patch(patch_type: &str, patch_type_id: usize, position_deg: f32) -> PatchDef {
        let rad = position_deg * PI / 180.0;
        PatchDef {
            patch_type: patch_type.to_string(),
            position_deg,
            position_rad: rad,
            position_cs: Vec2::new(rad.cos(), rad.sin()),
            patch_type_id,
        }
    }

    fn ori(deg: f32) -> Vec2 {
        let r = deg * PI / 180.0;
        Vec2::new(r.cos(), r.sin())
    }

    #[test]
    fn torque_is_nonzero_for_misaligned_ab_pair() {
        let i_patches = vec![patch("A", 0, 20.0)];
        let j_patches = vec![patch("B", 1, 160.0)];

        let mut lut = vec![PatchLutEntry::default(); 4];
        lut[0 * 2 + 1] = PatchLutEntry {
            epsilon: 10.0,
            sigma_i: 40.0 * PI / 180.0,
            sigma_j: 40.0 * PI / 180.0,
            cutoff: 1.35,
            enabled: true,
        };

        // Separation and orientations chosen to mirror testing_mode_demo first pair.
        let r_vec = Vec2::new(2.1, 0.0);
        let r_contact = 2.0;

        let (_f, tau_i, tau_j) = patchy_force_torque(
            r_vec,
            r_contact,
            &i_patches,
            ori(0.0),
            &j_patches,
            ori(35.0),
            &lut,
            2,
            0.0,
            3.0,
        );

        assert!(tau_i.abs() > 1e-6, "expected nonzero tau_i, got {tau_i}");
        assert!(tau_j.abs() > 1e-6, "expected nonzero tau_j, got {tau_j}");
    }

    #[test]
    fn demo_scene_has_nonzero_net_torque() {
        let patches = vec![patch("A", 0, 20.0), patch("B", 1, 160.0)];
        let mut lut = vec![PatchLutEntry::default(); 4];
        let sigma = 40.0 * PI / 180.0;
        lut[0 * 2 + 1] = PatchLutEntry {
            epsilon: 10.0,
            sigma_i: sigma,
            sigma_j: sigma,
            cutoff: 1.35,
            enabled: true,
        };
        lut[1 * 2 + 0] = PatchLutEntry {
            epsilon: 10.0,
            sigma_i: sigma,
            sigma_j: sigma,
            cutoff: 1.35,
            enabled: true,
        };

        let pos = [
            Vec2::new(0.0, 0.0),
            Vec2::new(2.1, 0.0),
            Vec2::new(4.0, 0.7),
            Vec2::new(5.8, 1.8),
        ];
        let oris = [ori(0.0), ori(35.0), ori(70.0), ori(110.0)];

        let mut torques = [0.0f32; 4];
        for i in 0..4 {
            for j in (i + 1)..4 {
                let r_vec = pos[j] - pos[i];
                let (_f, tau_i, tau_j) = patchy_force_torque(
                    r_vec,
                    2.0,
                    &patches,
                    oris[i],
                    &patches,
                    oris[j],
                    &lut,
                    2,
                    0.0,
                    3.0,
                );
                torques[i] += tau_i;
                torques[j] += tau_j;
            }
        }
        let max_abs = torques
            .iter()
            .fold(0.0f32, |m, &t| if t.abs() > m { t.abs() } else { m });
        assert!(
            max_abs > 1e-6,
            "expected nonzero scene torque, got torques={torques:?}"
        );
    }
}
