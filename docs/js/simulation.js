/**
 * simulation.js — thin wrapper around the CrystalSim WASM object.
 *
 * Isolates all direct WASM API calls from the rest of the frontend so that
 * the renderer and controls never import wasm-bindgen artefacts directly.
 */

let sim = null;
let wasmMemory = null;
let typeMetadata = [];  // [{color, radius}, ...]

/**
 * Initialise the simulation.
 * @param {object} wasmModule  — the imported wasm-bindgen module object
 * @param {object} config      — plain JS config object (will be JSON-stringified)
 */
export function initSim(wasmModule, config) {
  sim = new wasmModule.CrystalSim(JSON.stringify(config));
  wasmMemory = wasmModule.memory;
  // Derive metadata from the config directly — no round-trip through WASM JSON.
  typeMetadata = (config.particle_types || []).map(t => ({
    color:   t.color,
    radius:  t.radius,
    patches: (t.patches || []).map(p => ({
      angle_rad: p.position_deg * Math.PI / 180,
      color: '#ffffff',
    })),
  }));
}

/** Advance the simulation by n KMC events. */
export function stepSim(n) {
  if (sim) sim.step(n);
}

/** Number of live particles. */
export function getParticleCount() {
  return sim ? sim.particle_count() : 0;
}

/** Current KMC time. */
export function getTime() {
  return sim ? sim.simulation_time() : 0;
}

/**
 * Zero-copy Float32Array view into WASM memory.
 * Layout: [x, y, type_id, radius, orientation]  per particle (stride 5).
 */
export function getParticleBuffer() {
  if (!sim || !wasmMemory) return new Float32Array(0);
  const count = sim.particle_count();
  const ptr   = sim.particle_buffer();
  return new Float32Array(wasmMemory.buffer, ptr, count * 5);
}

/** Per-type metadata array (color strings, radii). */
export function getTypeMetadata() {
  return typeMetadata;
}

export function setTemperature(t) {
  if (sim) sim.set_temperature(t);
}

export function setChemicalPotential(typeId, mu) {
  if (sim) sim.set_chemical_potential(typeId, mu);
}

/** Number of candidate attachment sites. */
export function getCandidateCount() {
  return sim ? sim.candidate_count() : 0;
}

/**
 * Zero-copy Float32Array view into the candidate buffer.
 * Layout: [x, y, type_id, rate]  per candidate (stride 4).
 */
export function getCandidateBuffer() {
  if (!sim || !wasmMemory) return new Float32Array(0);
  const count = sim.candidate_count();
  const ptr   = sim.candidate_buffer();
  return new Float32Array(wasmMemory.buffer, ptr, count * 4);
}
