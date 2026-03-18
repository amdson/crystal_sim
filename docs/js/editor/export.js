/**
 * export.js — Build a SimConfig-compatible JSON object from editor state
 * and trigger a browser download.
 */

import { state, interactionKey, patchTypeById } from './state.js';

/**
 * Assemble the full config JSON object from current state.
 * The result is valid input for CrystalSim (the main simulation).
 */
export function buildConfigJson() {
  const { patchTypes, interactions, particleTypes, physics, instances } = state;

  // Particle types
  const particle_types = particleTypes.map(pt => ({
    radius: pt.radius,
    color: pt.color,
    mu: pt.mu,
    patches: pt.patches.map(patch => {
      const ptObj = patchTypeById(patch.patchTypeId);
      return {
        patch_type: ptObj?.name ?? 'A',
        position_deg: patch.position_deg,
      };
    }),
  }));

  // Epsilon matrix (NxN, filled with 0 — patch_interactions carries the real params)
  const n = particleTypes.length;
  const epsilon = Array.from({ length: n }, () => Array(n).fill(0));

  // Patch interactions (emit one entry per enabled pair)
  const patch_interactions = [];
  const emitted = new Set();
  for (const [key, val] of interactions) {
    if (emitted.has(key)) continue;
    emitted.add(key);
    if (val.epsilon === 0) continue;  // skip zero-epsilon pairs

    const [nameA, nameB] = key.split(':');
    // Determine canonical order from original key (already sorted alphabetically).
    patch_interactions.push({
      types: [nameA, nameB],
      epsilon: val.epsilon,
      angular_width_deg: val.angular_width_deg,
      cutoff: val.cutoff,
    });
  }

  // Initial particles (from Panel 3 instances)
  const initial_particles = instances.map(inst => ({
    x: inst.x,
    y: inst.y,
    type_id: inst.typeId,
    orientation_deg: inst.orientationDeg,
    frozen: inst.frozen,
  }));

  return {
    particle_types,
    epsilon,
    delta: physics.delta,
    temperature: physics.temperature,
    nu: physics.nu,
    seed: physics.seed,
    num_isolated_angles: physics.num_isolated_angles,
    max_arc_sites_per_type: physics.max_arc_sites_per_type,
    relax_steps: physics.relax_steps,
    relax_alpha: physics.relax_alpha,
    spring_k: physics.spring_k,
    lj_cutoff_factor: physics.lj_cutoff_factor,
    static_friction: physics.static_friction,
    steps_per_frame: physics.steps_per_frame,
    relax_damping: physics.relax_damping,
    patch_interactions,
    ...(initial_particles.length > 0 ? { initial_particles } : {}),
  };
}

/** Trigger a browser download of the config as a .json file. */
export function downloadJson(obj, filename = 'config.json') {
  const blob = new Blob([JSON.stringify(obj, null, 2)], { type: 'application/json' });
  const url  = URL.createObjectURL(blob);
  const a    = document.createElement('a');
  a.href     = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}
