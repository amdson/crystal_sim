/**
 * panel3.js — Physics preview canvas (Panel 3).
 *
 * Owns the live WASM EditorSim and runs the FIRE relaxation loop.
 * Accepts dropped particle type cards from Panel 2 and lets the user
 * drag existing instances to reposition them.
 */

import { state, dispatchChange } from './state.js';
import {
  renderPanel3Frame,
  clientToSim, hitTestParticle,
} from './renderer2d.js';
import { startDrag, isOverElement } from './drag.js';

let _canvas      = null;
let _ctx         = null;
let _editorSim   = null;
let _wasmMemory  = null;
let _paused      = false;
let _rafId       = null;
let _container   = null;

const camera = { x: 0, y: 0, scale: 50 };

// ── Initialization ────────────────────────────────────────────────────────────

export function initPanel3(container, editorSim, wasmMemory) {
  _container  = container;
  _editorSim  = editorSim;
  _wasmMemory = wasmMemory;

  buildUI(container);
  startLoop();
}

function buildUI(container) {
  container.innerHTML = '';

  const title = document.createElement('h2');
  title.className = 'panel-title';
  title.textContent = 'Physics Preview';
  container.appendChild(title);

  // Canvas wrapper
  const wrapper = document.createElement('div');
  wrapper.className = 'panel3-canvas-wrapper';
  _canvas = document.createElement('canvas');
  _canvas.className = 'panel3-canvas';
  wrapper.appendChild(_canvas);
  container.appendChild(wrapper);

  _ctx = _canvas.getContext('2d');

  // Sync canvas size to wrapper
  const resizeCanvas = () => {
    _canvas.width  = wrapper.clientWidth  || _canvas.offsetWidth;
    _canvas.height = wrapper.clientHeight || _canvas.offsetHeight;
  };
  window.addEventListener('resize', resizeCanvas);
  resizeCanvas();

  // Controls footer
  container.appendChild(buildFooter());

  // Pan + zoom + pointer interactions
  setupPointerInteractions();
}

function buildFooter() {
  const footer = document.createElement('div');
  footer.className = 'panel3-footer';

  // Pause/play button
  const pauseBtn = document.createElement('button');
  pauseBtn.className = 'btn';
  pauseBtn.textContent = '⏸ Pause';
  pauseBtn.addEventListener('click', () => {
    _paused = !_paused;
    pauseBtn.textContent = _paused ? '▶ Play' : '⏸ Pause';
  });
  footer.appendChild(pauseBtn);

  // Clear button
  const clearBtn = document.createElement('button');
  clearBtn.className = 'btn';
  clearBtn.textContent = '✕ Clear';
  clearBtn.addEventListener('click', () => {
    state.instances = [];
    _editorSim.set_particles('[]');
    dispatchChange('instances');
  });
  footer.appendChild(clearBtn);

  // Physics param inputs
  const params = [
    { label: 'Steps/frame', key: 'relax_steps',    step: 1,    min: 1 },
    { label: 'Scale',       key: null,              step: null, min: null, isScale: true },
  ];
  for (const p of params) {
    const wrap = document.createElement('label');
    wrap.className = 'footer-input-wrap';
    const lbl = document.createElement('span');
    lbl.textContent = p.label;
    const inp = document.createElement('input');
    inp.type  = 'number';
    if (p.isScale) {
      inp.step  = 5;
      inp.min   = 5;
      inp.value = camera.scale;
      inp.addEventListener('change', () => {
        camera.scale = Math.max(5, parseFloat(inp.value) || 50);
      });
    } else {
      inp.step  = p.step;
      inp.min   = p.min;
      inp.value = state.physics[p.key];
      inp.addEventListener('change', () => {
        state.physics[p.key] = parseFloat(inp.value) || state.physics[p.key];
      });
    }
    wrap.appendChild(lbl);
    wrap.appendChild(inp);
    footer.appendChild(wrap);
  }

  return footer;
}

// ── Animation loop ────────────────────────────────────────────────────────────

function startLoop() {
  if (_rafId) cancelAnimationFrame(_rafId);

  function frame() {
    _rafId = requestAnimationFrame(frame);
    if (!_editorSim || !_wasmMemory) return;

    if (!_paused) {
      _editorSim.relax(state.physics.relax_steps || 20);
    }

    const count = _editorSim.particle_count();
    const ptr   = _editorSim.particle_buffer();
    const buf   = new Float32Array(_wasmMemory.buffer, ptr, count * 5);

    const typeMeta = JSON.parse(_editorSim.type_metadata_json());

    renderPanel3Frame(
      _ctx,
      _canvas.width, _canvas.height,
      buf, count,
      typeMeta,
      camera,
    );
  }

  frame();
}

// ── Receive dropped particle from Panel 2 ────────────────────────────────────

/**
 * Called from main.js when a "particleType" drag is dropped onto Panel 3.
 */
export function receiveParticleDrop(typeId, clientX, clientY) {
  if (!isOverElement(clientX, clientY, _canvas)) return;
  if (!_editorSim) return;

  const sim = clientToSim(clientX, clientY, _canvas, camera);
  _editorSim.add_particle(sim.x, sim.y, typeId, 0, false);

  state.instances.push({ typeId, x: sim.x, y: sim.y, orientationDeg: 0, frozen: false });
  dispatchChange('instances');
}

/** Returns the panel 3 canvas element (used for drop-zone testing in main.js). */
export function getPanel3Canvas() { return _canvas; }

// ── Pointer interactions (pan, zoom, drag instances) ─────────────────────────

function setupPointerInteractions() {
  const c = _canvas;
  if (!c) return;

  // Zoom
  c.addEventListener('wheel', e => {
    e.preventDefault();
    const factor = e.deltaY < 0 ? 1.1 : 0.9;
    camera.scale = Math.min(Math.max(camera.scale * factor, 5), 300);
  }, { passive: false });

  let _panStart = null;
  let _dragInstanceIdx = -1;

  c.addEventListener('pointerdown', e => {
    if (!_editorSim) return;
    e.preventDefault();

    const count  = _editorSim.particle_count();
    const ptr    = _editorSim.particle_buffer();
    const buf    = new Float32Array(_wasmMemory.buffer, ptr, count * 5);
    const simPos = clientToSim(e.clientX, e.clientY, c, camera);
    const hit    = hitTestParticle(simPos.x, simPos.y, buf, count);

    if (hit.index !== -1) {
      // Drag instance
      _dragInstanceIdx = hit.index;
      _editorSim.reset_velocities();
      c.setPointerCapture(e.pointerId);

      const onMove = ev => {
        const s = clientToSim(ev.clientX, ev.clientY, c, camera);
        // Update WASM position by rebuilding all instances from WASM buffer
        // (simpler than tracking state separately)
        const n    = _editorSim.particle_count();
        const p2   = _editorSim.particle_buffer();
        const b2   = new Float32Array(_wasmMemory.buffer, p2, n * 5);
        const placements = [];
        for (let i = 0; i < n; i++) {
          placements.push({
            x: i === _dragInstanceIdx ? s.x : b2[i * 5],
            y: i === _dragInstanceIdx ? s.y : b2[i * 5 + 1],
            type_id: Math.round(b2[i * 5 + 2]),
            orientation_deg: b2[i * 5 + 4] * 180 / Math.PI,
            frozen: false,
          });
        }
        _editorSim.set_particles(JSON.stringify(placements));
        _editorSim.reset_velocities();
      };

      const onUp = () => {
        _dragInstanceIdx = -1;
        c.removeEventListener('pointermove', onMove);
        c.removeEventListener('pointerup', onUp);
      };

      c.addEventListener('pointermove', onMove);
      c.addEventListener('pointerup', onUp);
    } else {
      // Pan
      _panStart = { x: e.clientX, y: e.clientY, cx: camera.x, cy: camera.y };
      c.setPointerCapture(e.pointerId);

      const onMove = ev => {
        if (!_panStart) return;
        const dx = ev.clientX - _panStart.x;
        const dy = ev.clientY - _panStart.y;
        camera.x = _panStart.cx - dx / camera.scale;
        camera.y = _panStart.cy + dy / camera.scale;
      };
      const onUp = () => {
        _panStart = null;
        c.removeEventListener('pointermove', onMove);
        c.removeEventListener('pointerup', onUp);
      };
      c.addEventListener('pointermove', onMove);
      c.addEventListener('pointerup', onUp);
    }
  });
}

// ── Config hot-reload ─────────────────────────────────────────────────────────

/**
 * Push a new config JSON string to the WASM sim without clearing particles.
 * Called from main.js whenever panel 1 or 2 changes.
 */
export function pushConfig(configJson) {
  if (!_editorSim) return;
  try {
    _editorSim.update_config(configJson);
  } catch (e) {
    console.warn('EditorSim.update_config failed:', e);
  }
}
