//! EditorSim — lightweight physics sandbox for the browser config editor.
//!
//! Runs a global FIRE-style damped gradient descent over all particles,
//! using the same patchy_force_torque kernel as the main simulation's
//! local relaxation (relax_new_particle_fixed_grid in kmc.rs).
//!
//! Unlike the main simulation this has no KMC engine, no spatial hash,
//! and no candidate-site machinery — just a flat particle list and an
//! O(N²) pairwise force loop, which is fine for the small interactive
//! scenes the editor produces (typically ≤50 particles).

use glam::Vec2;
use serde::Deserialize;
use wasm_bindgen::prelude::*;

use crate::config::SimConfig;
use crate::forces::patchy_force_torque;

// ── Internal particle state ────────────────────────────────────────────────

struct EditorParticle {
    pos: Vec2,
    /// Unit vector encoding orientation angle: angle = atan2(y, x).
    orientation: Vec2,
    type_id: usize,
    radius: f32,
    frozen: bool,
}

// ── EditorSim ─────────────────────────────────────────────────────────────

/// Lightweight physics sandbox exposed to the browser config editor.
///
/// JS usage:
/// ```js
/// import init, { EditorSim } from './pkg/crystal_sim.js';
/// await init();
/// const sim = new EditorSim(JSON.stringify(config));
/// sim.add_particle(0, 0, 0, 0, false);
/// sim.relax(30);
/// const buf = new Float32Array(memory.buffer, sim.particle_buffer(), sim.particle_count() * 5);
/// ```
#[wasm_bindgen]
pub struct EditorSim {
    config: SimConfig,
    particles: Vec<EditorParticle>,
    /// Flat stride-5 buffer: [x, y, type_id, radius, orientation_angle] per particle.
    /// Same layout as CrystalSim so renderer2d.js can be shared.
    particle_buf: Vec<f32>,
    // Per-particle scratch space reused each step (avoids repeated allocations).
    velocities: Vec<Vec2>,
    ang_velocities: Vec<f32>,
    forces: Vec<Vec2>,
    torques: Vec<f32>,
}

impl EditorSim {
    /// Rebuild the flat buffer from current particle state.
    fn sync_buf(&mut self) {
        self.particle_buf.clear();
        for p in &self.particles {
            self.particle_buf.push(p.pos.x);
            self.particle_buf.push(p.pos.y);
            self.particle_buf.push(p.type_id as f32);
            self.particle_buf.push(p.radius);
            self.particle_buf.push(p.orientation.y.atan2(p.orientation.x));
        }
    }

    /// Safe fallback epsilon lookup (returns 0 if out of bounds).
    fn epsilon_safe(&self, i: usize, j: usize) -> f32 {
        self.config
            .epsilon
            .get(i)
            .and_then(|row| row.get(j))
            .copied()
            .unwrap_or(0.0)
    }

    /// Clamp type_id to the valid range; returns 0 if config has no types.
    fn clamp_type(&self, type_id: usize) -> usize {
        let n = self.config.particle_types.len();
        if n == 0 { 0 } else { type_id.min(n - 1) }
    }

    /// Grow scratch vecs to match current particle count if needed.
    fn ensure_scratch(&mut self) {
        let n = self.particles.len();
        if self.velocities.len() < n {
            self.velocities.resize(n, Vec2::ZERO);
            self.ang_velocities.resize(n, 0.0);
            self.forces.resize(n, Vec2::ZERO);
            self.torques.resize(n, 0.0);
        }
    }
}

// ── Deserialization helpers ────────────────────────────────────────────────

#[derive(Deserialize)]
struct PlacementDef {
    x: f32,
    y: f32,
    type_id: usize,
    #[serde(default)]
    orientation_deg: f32,
    #[serde(default)]
    frozen: bool,
}

// ── WASM API ───────────────────────────────────────────────────────────────

#[wasm_bindgen]
impl EditorSim {
    /// Construct from a JSON config string. Starts with zero particles.
    #[wasm_bindgen(constructor)]
    pub fn new(config_json: &str) -> Result<EditorSim, JsValue> {
        let mut config: SimConfig = serde_json::from_str(config_json)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        config.init_cache();
        Ok(EditorSim {
            config,
            particles: Vec::new(),
            particle_buf: Vec::with_capacity(256 * 5),
            velocities: Vec::new(),
            ang_velocities: Vec::new(),
            forces: Vec::new(),
            torques: Vec::new(),
        })
    }

    /// Replace all particles from a JSON array.
    /// Schema: `[{"x":f,"y":f,"type_id":u,"orientation_deg":f,"frozen":b}, ...]`
    /// All velocities are zeroed.
    pub fn set_particles(&mut self, placements_json: &str) -> Result<(), JsValue> {
        let placements: Vec<PlacementDef> = serde_json::from_str(placements_json)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;

        self.particles.clear();
        self.velocities.clear();
        self.ang_velocities.clear();
        self.forces.clear();
        self.torques.clear();

        for p in placements {
            let type_id = self.clamp_type(p.type_id);
            let rad = p.orientation_deg * std::f32::consts::PI / 180.0;
            self.particles.push(EditorParticle {
                pos: Vec2::new(p.x, p.y),
                orientation: Vec2::new(rad.cos(), rad.sin()),
                type_id,
                radius: self.config.particle_types[type_id].radius,
                frozen: p.frozen,
            });
            self.velocities.push(Vec2::ZERO);
            self.ang_velocities.push(0.0);
            self.forces.push(Vec2::ZERO);
            self.torques.push(0.0);
        }

        self.sync_buf();
        Ok(())
    }

    /// Add one particle instance. Returns its index.
    pub fn add_particle(
        &mut self,
        x: f32,
        y: f32,
        type_id: u32,
        orientation_deg: f32,
        frozen: bool,
    ) -> u32 {
        let type_id = self.clamp_type(type_id as usize);
        let rad = orientation_deg * std::f32::consts::PI / 180.0;
        let idx = self.particles.len() as u32;
        self.particles.push(EditorParticle {
            pos: Vec2::new(x, y),
            orientation: Vec2::new(rad.cos(), rad.sin()),
            type_id,
            radius: self.config.particle_types[type_id].radius,
            frozen,
        });
        self.velocities.push(Vec2::ZERO);
        self.ang_velocities.push(0.0);
        self.forces.push(Vec2::ZERO);
        self.torques.push(0.0);
        self.sync_buf();
        idx
    }

    /// Swap-remove particle at index `i`. Index stability not guaranteed.
    pub fn remove_particle(&mut self, i: u32) {
        let i = i as usize;
        if i < self.particles.len() {
            self.particles.swap_remove(i);
            self.velocities.swap_remove(i);
            self.ang_velocities.swap_remove(i);
            self.forces.swap_remove(i);
            self.torques.swap_remove(i);
            self.sync_buf();
        }
    }

    /// Zero all particle velocities (call after a manual drag-move so relaxation
    /// resumes from rest rather than from stale momentum).
    pub fn reset_velocities(&mut self) {
        for v in self.velocities.iter_mut() { *v = Vec2::ZERO; }
        for w in self.ang_velocities.iter_mut() { *w = 0.0; }
    }

    /// Replace the config (e.g. after the user edits interaction parameters)
    /// without clearing the current particle layout.
    pub fn update_config(&mut self, config_json: &str) -> Result<(), JsValue> {
        let mut cfg: SimConfig = serde_json::from_str(config_json)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        cfg.init_cache();

        // Re-clamp type_ids and radii in case the number of types changed.
        for p in self.particles.iter_mut() {
            p.type_id = p.type_id.min(cfg.particle_types.len().saturating_sub(1));
            if !cfg.particle_types.is_empty() {
                p.radius = cfg.particle_types[p.type_id].radius;
            }
        }
        self.config = cfg;
        self.sync_buf();
        Ok(())
    }

    /// Run up to `n` FIRE-style relaxation steps.
    ///
    /// This is a global O(N²) version of the local FIRE loop in kmc.rs
    /// (`relax_new_particle_fixed_grid`). The velocity integration math
    /// is identical; the only difference is that all particles participate
    /// rather than just those in a 5×5 cell neighbourhood.
    pub fn relax(&mut self, n: u32) {
        let count = self.particles.len();
        if count == 0 || n == 0 { return; }

        self.ensure_scratch();

        let alpha           = self.config.relax_alpha;
        let damping         = self.config.relax_damping;
        let static_friction = self.config.static_friction;
        let lj_cutoff       = self.config.lj_cutoff_factor;
        let patch_tc        = self.config.patch_type_count;

        for _ in 0..n {
            // ── Zero forces ────────────────────────────────────────────────
            for f in self.forces[..count].iter_mut() { *f = Vec2::ZERO; }
            for t in self.torques[..count].iter_mut() { *t = 0.0; }

            // ── All-pairs patchy force + torque accumulation ────────────────
            for i in 0..count {
                for j in (i + 1)..count {
                    let r_vec = self.particles[j].pos - self.particles[i].pos;
                    let r_c   = self.particles[i].radius + self.particles[j].radius;
                    let r_cut = r_c * lj_cutoff;
                    if r_vec.length_squared() >= r_cut * r_cut { continue; }

                    let ti = self.particles[i].type_id;
                    let tj = self.particles[j].type_id;

                    let (f, tau_i, tau_j) = patchy_force_torque(
                        r_vec,
                        r_c,
                        &self.config.particle_types[ti].patches,
                        self.particles[i].orientation,
                        &self.config.particle_types[tj].patches,
                        self.particles[j].orientation,
                        &self.config.patch_lut,
                        patch_tc,
                        self.epsilon_safe(ti, tj),
                        lj_cutoff,
                    );

                    self.forces[i]  += f;
                    self.forces[j]  -= f;
                    self.torques[i] += tau_i;
                    self.torques[j] += tau_j;
                }
            }

            // ── FIRE-style velocity update (mirrors kmc.rs lines 1063-1099) ─
            let mut converged = true;
            for i in 0..count {
                if self.particles[i].frozen { continue; }

                // Linear
                let f_raw     = self.forces[i].clamp_length_max(3.0);
                let f_mag     = f_raw.length();
                let f_eff     = if f_mag > static_friction {
                    f_raw * ((f_mag - static_friction) / f_mag)
                } else {
                    Vec2::ZERO
                };

                let v_raw = damping * self.velocities[i] + alpha * f_eff;
                let v_new = if v_raw.dot(f_eff) >= 0.0 { v_raw } else { Vec2::ZERO };

                self.particles[i].pos += v_new;
                self.velocities[i]     = v_new;

                // Angular
                let tau_raw    = self.torques[i].clamp(-3.0, 3.0);
                let tau_eff    = if tau_raw.abs() > static_friction { tau_raw } else { 0.0 };
                let omega_raw  = damping * self.ang_velocities[i] + alpha * tau_eff;
                let omega_new  = if omega_raw * tau_eff >= 0.0 { omega_raw } else { 0.0 };

                let (s, c)  = omega_new.sin_cos();
                let ori     = self.particles[i].orientation;
                let new_ori = Vec2::new(ori.x * c - ori.y * s, ori.x * s + ori.y * c);
                let len_sq  = new_ori.length_squared();
                self.particles[i].orientation = if (len_sq - 1.0).abs() > 1e-6 {
                    new_ori / len_sq.sqrt()
                } else {
                    new_ori
                };
                self.ang_velocities[i] = omega_new;

                if f_eff.length_squared() >= 0.05 * 0.05 || tau_eff.abs() >= 0.05 {
                    converged = false;
                }
            }

            if converged { break; }
        }

        self.sync_buf();
    }

    /// Byte offset into WASM linear memory of the stride-5 particle buffer.
    /// Valid until any method that may reallocate the buffer is called.
    ///
    /// JS: `new Float32Array(memory.buffer, sim.particle_buffer(), sim.particle_count() * 5)`
    pub fn particle_buffer(&self) -> u32 {
        self.particle_buf.as_ptr() as u32
    }

    /// Number of particles currently in the sim.
    pub fn particle_count(&self) -> u32 {
        self.particles.len() as u32
    }

    /// JSON array of per-type metadata for the renderer.
    /// Format: `[{"color":"#hex","radius":r,"patches":[{"angle_rad":f,"color":"#fff"},...]}]`
    /// Matches the format returned by `CrystalSim::type_metadata_json`.
    pub fn type_metadata_json(&self) -> String {
        let entries: Vec<String> = self
            .config
            .particle_types
            .iter()
            .map(|t| {
                let patches: Vec<String> = t
                    .patches
                    .iter()
                    .map(|p| {
                        format!(
                            r##"{{"angle_rad":{:.6},"color":"#ffffff"}}"##,
                            p.position_rad
                        )
                    })
                    .collect();
                format!(
                    r#"{{"color":"{}","radius":{},"patches":[{}]}}"#,
                    t.color,
                    t.radius,
                    patches.join(",")
                )
            })
            .collect();
        format!("[{}]", entries.join(","))
    }
}
