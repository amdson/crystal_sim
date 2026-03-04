/**
 * main.js — entry point.
 *
 * Loads the WASM module, wires up simulation + renderer + controls,
 * then runs the requestAnimationFrame render/step loop.
 */

import init, { CrystalSim } from '../../pkg/crystal_sim.js';
import * as simAPI from './simulation.js';
import { Renderer } from './renderer.js';
import { initControls } from './controls.js';

// ── Default configuration ─────────────────────────────────────────────────
// Edit this object (or load from a JSON file) to change the simulation setup.
const DEFAULT_CONFIG = {
  particle_types: [
    { radius: 1.0,  color: '#5b9bd5', mu: -1.5 },
    { radius: 0.75, color: '#ed7d31', mu: -1.5 },
  ],
  // 2×2 interaction matrix — symmetric, units of kT
  epsilon: [
    [2.5, 1.5],
    [1.5, 2.0],
  ],
  delta: 0.3,       // bonding shell width
  temperature: 1.0, // kT
  nu: 1.0,          // attempt frequency
  seed: 42,
  num_isolated_angles: 16,
};

// ── State ─────────────────────────────────────────────────────────────────
const canvas   = document.getElementById('canvas');
const renderer = new Renderer(canvas);
renderer.paused = false;

let controls = null;
let wasmModule = null;

// ── Startup ───────────────────────────────────────────────────────────────
async function start() {
  wasmModule = await init();          // initialise WASM binary
  boot(DEFAULT_CONFIG);
}

function boot(config) {
  // (Re-)initialise simulation
  simAPI.initSim({ CrystalSim, memory: wasmModule.memory }, config);

  // Wire controls on first boot only
  if (!controls) {
    controls = initControls(simAPI, renderer, config, () => boot(config));
  }

  // Start (or restart) the render loop
  requestAnimationFrame(frame);
}

// ── FPS counter ───────────────────────────────────────────────────────────
let lastTime = performance.now();
let frameCount = 0;

function frame(now) {
  requestAnimationFrame(frame);

  // FPS
  frameCount++;
  if (now - lastTime >= 1000) {
    document.getElementById('statFps').textContent = frameCount;
    frameCount = 0;
    lastTime = now;
  }

  // Advance simulation
  if (!renderer.paused && controls) {
    simAPI.stepSim(controls.getStepsPerFrame());
  }

  // Render
  const count = simAPI.getParticleCount();
  const buf   = simAPI.getParticleBuffer();
  const meta  = simAPI.getTypeMetadata();
  renderer.draw(buf, count, meta);

  // Update stats
  document.getElementById('statParticles').textContent = count;
  document.getElementById('statTime').textContent =
    simAPI.getTime().toExponential(3);
}

start().catch(console.error);
