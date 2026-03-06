//! Force and potential functions for particle relaxation.
//!
//! All functions follow the same vector convention used in the relaxation loop:
//!   `r_vec = pos_j - pos_i`  (points from i toward j)
//! The returned force vector acts on particle i:
//!   - component along +r_vec → attractive (toward j)
//!   - component along −r_vec → repulsive (away from j)

use glam::DVec2;

/// 2^(−1/6): scale factor so that the LJ minimum falls exactly at r_contact.
const SIGMA_FACTOR: f64 = 0.890_898_718_36; // 2^(-1/6)

// ── Lennard-Jones ─────────────────────────────────────────────────────────────

/// Lennard-Jones potential (no cutoff).
///
/// The energy minimum is placed at `r = r_contact` with depth `−epsilon`.
/// When `epsilon < 0` the well becomes a barrier (purely repulsive pair).
///
/// `U(r) = 4ε[(σ/r)^12 − (σ/r)^6]`  where `σ = r_contact · 2^(−1/6)`
pub fn lj_potential(r: f64, r_contact: f64, epsilon: f64) -> f64 {
    let sigma = r_contact * SIGMA_FACTOR;
    let sr = sigma / r;
    let sr6 = sr * sr * sr * sr * sr * sr;
    let sr12 = sr6 * sr6;
    4.0 * epsilon * (sr12 - sr6)
}

/// Lennard-Jones force vector on particle i due to particle j.
///
/// Returns `DVec2::ZERO` when `r > r_cut`.
///
/// Derived from `U_LJ` above:
///   `F(r) = 24ε/r² · [2(σ/r)^12 − (σ/r)^6] · r_vec`
/// (positive along r_vec → repulsive at short range, negative → attractive at long range)
pub fn lj_force_vec(r_vec: DVec2, r_contact: f64, epsilon: f64, r_cut: f64) -> DVec2 {
    let r2 = r_vec.length_squared();
    if r2 >= r_cut * r_cut || r2 < 1e-24 {
        return DVec2::ZERO;
    }
    let sigma = r_contact * SIGMA_FACTOR;
    let sigma2 = sigma * sigma;
    let sr2 = sigma2 / r2;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    // Factor: 24ε/r² · [2(σ/r)^12 − (σ/r)^6]
    // Positive → along r_vec (repulsive), negative → against r_vec (attractive)
    let coeff = 24.0 * epsilon / r2 * (2.0 * sr12 - sr6);
    coeff * r_vec
}

// ── Harmonic spring ───────────────────────────────────────────────────────────

/// Harmonic spring force vector on particle i due to particle j.
///
/// `F = 2k(r − r_contact) · r̂`
///   - repulsive when `r < r_contact`
///   - attractive when `r > r_contact`
///
/// Returns `DVec2::ZERO` when `r > r_contact + shell` or `r < 1e-12`.
pub fn spring_force_vec(r_vec: DVec2, r_contact: f64, k: f64, shell: f64) -> DVec2 {
    let r = r_vec.length();
    if r < 1e-12 || r > r_contact + shell {
        return DVec2::ZERO;
    }
    let r_hat = r_vec / r;
    (2.0 * k * (r - r_contact)) * r_hat
}
