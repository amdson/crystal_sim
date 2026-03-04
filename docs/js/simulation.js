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
  typeMetadata = JSON.parse(sim.type_metadata_json());
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
 * Layout: [x, y, type_id, radius]  per particle (stride 4).
 */
export function getParticleBuffer() {
  if (!sim || !wasmMemory) return new Float32Array(0);
  const count = sim.particle_count();
  const ptr   = sim.particle_buffer();
  return new Float32Array(wasmMemory.buffer, ptr, count * 4);
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
