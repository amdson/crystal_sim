/**
 * state.js — Single source of truth for the editor.
 *
 * All panels read from and write to this object. After mutating state,
 * call dispatchChange(section) to notify registered listeners.
 */

export const state = {
  /** Patch type definitions: [{id, name, color}] */
  patchTypes: [],

  /**
   * Pairwise interaction parameters.
   * Key: canonical "A:B" string (sorted alphabetically).
   * Value: {epsilon, angular_width_deg: [f, f], cutoff}
   */
  interactions: new Map(),

  /**
   * Particle type definitions:
   * [{id, name, color, radius, mu, patches: [{patchTypeId, position_deg}]}]
   */
  particleTypes: [],

  /** Physics / simulation parameters */
  physics: {
    temperature: 1.0,
    nu: 0.05,
    delta: 0.12,
    seed: 42,
    lj_cutoff_factor: 3.0,
    relax_alpha: 0.015,
    relax_damping: 0.8,
    static_friction: 0.12,
    relax_steps: 20,
    steps_per_frame: 10,
    num_isolated_angles: 16,
    max_arc_sites_per_type: 6,
    spring_k: 50.0,
  },

  /**
   * Panel 3 particle instances:
   * [{typeId, x, y, orientationDeg, frozen}]
   */
  instances: [],
};

// ── Change notification ──────────────────────────────────────────────────────

const _listeners = { patches: [], particles: [], physics: [], instances: [], all: [] };

export function onStateChange(section, fn) {
  (_listeners[section] ?? _listeners.all).push(fn);
}

export function dispatchChange(section) {
  for (const fn of (_listeners[section] ?? [])) fn();
  for (const fn of _listeners.all) fn();
}

// ── Interaction key helpers ──────────────────────────────────────────────────

/** Canonical map key for a patch pair (order-independent). */
export function interactionKey(nameA, nameB) {
  return [nameA, nameB].sort().join(':');
}

export function getInteraction(nameA, nameB) {
  return state.interactions.get(interactionKey(nameA, nameB)) ?? null;
}

export function setInteraction(nameA, nameB, value) {
  state.interactions.set(interactionKey(nameA, nameB), value);
}

// ── Auto-naming helpers ──────────────────────────────────────────────────────

const PATCH_COLORS = [
  '#e84040', '#4a90d9', '#f5a623', '#7ed321', '#bd10e0',
  '#ff6b6b', '#50e3c2', '#f8e71c', '#9b9b9b', '#ffffff',
];

/** Assign the next unused uppercase letter name for a new patch type. */
export function nextPatchName() {
  const used = new Set(state.patchTypes.map(p => p.name));
  for (let c = 65; c <= 90; c++) {
    const name = String.fromCharCode(c);
    if (!used.has(name)) return name;
  }
  return `P${state.patchTypes.length}`;
}

export function patchColorFor(index) {
  return PATCH_COLORS[index % PATCH_COLORS.length];
}

/** Return the patch type object by name, or undefined. */
export function patchTypeByName(name) {
  return state.patchTypes.find(p => p.name === name);
}

export function patchTypeById(id) {
  return state.patchTypes.find(p => p.id === id);
}

// ── Load an existing config JSON object into state ───────────────────────────

export function loadConfig(json) {
  state.patchTypes = [];
  state.interactions = new Map();
  state.particleTypes = [];
  state.instances = [];

  // Collect patch type names from particle definitions.
  const patchNameOrder = [];
  const patchNameSet = new Set();
  for (const pt of (json.particle_types ?? [])) {
    for (const patch of (pt.patches ?? [])) {
      if (!patchNameSet.has(patch.patch_type)) {
        patchNameSet.add(patch.patch_type);
        patchNameOrder.push(patch.patch_type);
      }
    }
  }
  // Also ensure names from patch_interactions are included.
  for (const pi of (json.patch_interactions ?? [])) {
    for (const t of pi.types) {
      if (!patchNameSet.has(t)) {
        patchNameSet.add(t);
        patchNameOrder.push(t);
      }
    }
  }

  // Build patchTypes array.
  patchNameOrder.forEach((name, idx) => {
    state.patchTypes.push({ id: idx, name, color: patchColorFor(idx) });
  });

  // Build interactions map.
  for (const pi of (json.patch_interactions ?? [])) {
    const key = interactionKey(pi.types[0], pi.types[1]);
    state.interactions.set(key, {
      epsilon: pi.epsilon ?? 0,
      angular_width_deg: pi.angular_width_deg ?? [30, 30],
      cutoff: pi.cutoff ?? 1.35,
    });
  }

  // Build particle types.
  (json.particle_types ?? []).forEach((pt, idx) => {
    const patches = (pt.patches ?? []).map(patch => {
      const ptObj = state.patchTypes.find(p => p.name === patch.patch_type);
      return { patchTypeId: ptObj?.id ?? 0, position_deg: patch.position_deg ?? 0 };
    });
    state.particleTypes.push({
      id: idx,
      name: `Type ${idx}`,
      color: pt.color ?? '#8888cc',
      radius: pt.radius ?? 1.0,
      mu: pt.mu ?? -2.0,
      patches,
    });
  });

  // Load physics params.
  const P = state.physics;
  P.temperature         = json.temperature         ?? P.temperature;
  P.nu                  = json.nu                  ?? P.nu;
  P.delta               = json.delta               ?? P.delta;
  P.seed                = json.seed                ?? P.seed;
  P.lj_cutoff_factor    = json.lj_cutoff_factor    ?? P.lj_cutoff_factor;
  P.relax_alpha         = json.relax_alpha         ?? P.relax_alpha;
  P.relax_damping       = json.relax_damping       ?? P.relax_damping;
  P.static_friction     = json.static_friction     ?? P.static_friction;
  P.relax_steps         = json.relax_steps         ?? P.relax_steps;
  P.steps_per_frame     = json.steps_per_frame     ?? P.steps_per_frame;
  P.num_isolated_angles = json.num_isolated_angles ?? P.num_isolated_angles;
  P.max_arc_sites_per_type = json.max_arc_sites_per_type ?? P.max_arc_sites_per_type;
  P.spring_k            = json.spring_k            ?? P.spring_k;

  // Load initial particles as editor instances.
  for (const ip of (json.initial_particles ?? [])) {
    state.instances.push({
      typeId: ip.type_id ?? 0,
      x: ip.x ?? 0,
      y: ip.y ?? 0,
      orientationDeg: ip.orientation_deg ?? 0,
      frozen: ip.frozen ?? false,
    });
  }
}
