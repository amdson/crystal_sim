use std::f32::consts::PI;

use glam::Vec2;

use crystal_sim::config::{PatchDef, PatchLutEntry};
use crystal_sim::forces::{
    lj_force_vec, lj_force_vec_full, lj_potential, patchy_force_torque, patchy_pair_energy,
    repulsive_force_vec, spring_force_vec,
};

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

fn ab_lut() -> (Vec<PatchLutEntry>, usize) {
    let n = 2usize;
    let mut lut = vec![PatchLutEntry::default(); n * n];
    let sigma = 20.0 * PI / 180.0;
    lut[0 * n + 1] = PatchLutEntry {
        epsilon: 8.0,
        sigma_i: sigma,
        sigma_j: sigma,
        cutoff: 1.35,
        enabled: true,
    };
    lut[1 * n + 0] = PatchLutEntry {
        epsilon: 8.0,
        sigma_i: sigma,
        sigma_j: sigma,
        cutoff: 1.35,
        enabled: true,
    };
    (lut, n)
}

#[test]
fn runs_lj_potential_and_force_variants() {
    let r_contact = 2.0;
    let r_vec = Vec2::new(2.1, 0.1);
    let r_cut = 3.0;

    let u = lj_potential(r_vec.length(), r_contact, 1.5);
    let g_legacy = lj_force_vec(r_vec, r_contact, 1.5, r_cut);
    let g_full = lj_force_vec_full(r_vec, r_contact, 1.5, r_cut);
    let rep = repulsive_force_vec(r_vec, r_contact, r_cut);

    println!("lj_potential={u:.6e}");
    println!("lj_force_vec={g_legacy:?}");
    println!("lj_force_vec_full={g_full:?}");
    println!("repulsive_force_vec={rep:?}");

    assert!(u.is_finite());
    assert!(g_legacy.is_finite() && g_legacy != Vec2::ZERO);
    assert!(g_full.is_finite() && g_full != Vec2::ZERO);
    assert!(rep.is_finite() && rep != Vec2::ZERO);
}

#[test]
fn runs_spring_force() {
    let r_contact = 2.0;
    let shell = 0.4;
    let k = 15.0;
    let in_shell = spring_force_vec(Vec2::new(2.2, 0.0), r_contact, k, shell);
    let out_shell = spring_force_vec(Vec2::new(2.6, 0.0), r_contact, k, shell);

    println!("spring_in_shell={in_shell:?}");
    println!("spring_out_shell={out_shell:?}");

    assert!(in_shell.is_finite() && in_shell != Vec2::ZERO);
    assert_eq!(out_shell, Vec2::ZERO);
}

#[test]
fn runs_patchy_pair_energy_fallback_and_patchy() {
    let (lut, n_patch_types) = ab_lut();
    let r_contact = 2.0;
    let r_vec = Vec2::new(2.05, 0.0);
    let delta = 0.15;

    // Fallback path: one side has no patches.
    let e_fallback = patchy_pair_energy(
        r_vec,
        r_contact,
        &[],
        Vec2::X,
        &[],
        Vec2::X,
        &lut,
        n_patch_types,
        3.0,
        delta,
    );

    // Patchy path: compatible A-B with near alignment.
    let i_patches = vec![patch("A", 0, 0.0)];
    let j_patches = vec![patch("B", 1, 180.0)];
    let e_patch = patchy_pair_energy(
        r_vec,
        r_contact,
        &i_patches,
        ori(0.0),
        &j_patches,
        ori(0.0),
        &lut,
        n_patch_types,
        0.0,
        delta,
    );

    println!("fallback_energy={e_fallback:.6e}");
    println!("patchy_energy={e_patch:.6e}");

    assert!(e_fallback > 0.0);
    assert!(e_patch.is_finite());
}

#[test]
fn runs_patchy_force_torque() {
    let (lut, n_patch_types) = ab_lut();
    let r_contact = 2.0;
    let r_vec = Vec2::new(2.1, 0.0);
    let i_patches = vec![patch("A", 0, 20.0)];
    let j_patches = vec![patch("B", 1, 160.0)];

    let (f, tau_i, tau_j) = patchy_force_torque(
        r_vec,
        r_contact,
        &i_patches,
        ori(0.0),
        &j_patches,
        ori(35.0),
        &lut,
        n_patch_types,
        0.0,
        3.0,
    );

    println!("patchy_force={f:?}");
    println!("tau_i={tau_i:.6e}, tau_j={tau_j:.6e}");

    assert!(f.is_finite());
    assert!(tau_i.is_finite() && tau_j.is_finite());
    assert!(tau_i.abs() > 1e-8 || tau_j.abs() > 1e-8);
}
