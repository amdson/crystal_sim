/**
 * main.js — Editor entry point.
 *
 * Loads the WASM module, initialises the three panels, wires inter-panel
 * drag-and-drop events, and connects the Export / Load Config buttons.
 */

import init, { EditorSim } from '../../pkg/crystal_sim.js';

import {
  state, dispatchChange, loadConfig,
} from './state.js';
import { buildConfigJson, downloadJson } from './export.js';
import { initPanel1, refreshPanel1 } from './panel1.js';
import {
  initPanel2, refreshPanel2, redrawAllCards,
  receivePatchDrop, updatePatchDropPreview,
} from './panel2.js';
import { initPanel3, receiveParticleDrop, getPanel3Canvas, pushConfig } from './panel3.js';
import { isDragging, getActiveDrag, isOverElement } from './drag.js';

let _wasmModule = null;
let _editorSim  = null;

// ── Default starter config ────────────────────────────────────────────────────

const DEFAULT_CONFIG = {
  particle_types: [],
  epsilon: [],
  delta: 0.12,
  temperature: 1.0,
  nu: 0.05,
  seed: 42,
  patch_interactions: [],
};

// ── Bootstrap ─────────────────────────────────────────────────────────────────

async function main() {
  _wasmModule = await init();

  // Create a minimal EditorSim from the default config
  _editorSim = new EditorSim(JSON.stringify(DEFAULT_CONFIG));

  // Get WASM linear memory reference
  // The memory object is on the wasm import object returned by init()
  const wasmMemory = _wasmModule.memory;

  // Initialise panels
  initPanel1(document.getElementById('panel1'), onConfigChange);
  initPanel2(document.getElementById('panel2'), onConfigChange);
  initPanel3(document.getElementById('panel3'), _editorSim, wasmMemory);

  // Wire export button
  document.getElementById('btn-export').addEventListener('click', () => {
    downloadJson(buildConfigJson());
  });

  // Wire load config button
  document.getElementById('btn-load').addEventListener('click', () => {
    document.getElementById('file-input').click();
  });
  document.getElementById('file-input').addEventListener('change', async e => {
    const file = e.target.files[0];
    if (!file) return;
    const text = await file.text();
    try {
      const json = JSON.parse(text);
      loadConfig(json);
      // Rebuild EditorSim with new config
      const cfgJson = JSON.stringify(buildConfigJson());
      _editorSim = new EditorSim(cfgJson);
      // Re-init Panel 3 with the new sim instance
      initPanel3(document.getElementById('panel3'), _editorSim, wasmMemory);
      // Place initial particles if any
      if (state.instances.length > 0) {
        _editorSim.set_particles(JSON.stringify(
          state.instances.map(i => ({
            x: i.x, y: i.y, type_id: i.typeId,
            orientation_deg: i.orientationDeg, frozen: i.frozen,
          }))
        ));
      }
      refreshPanel1();
      refreshPanel2();
      redrawAllCards();
    } catch (err) {
      alert('Failed to parse config: ' + err.message);
    }
    e.target.value = '';
  });

  // Wire global pointer events for cross-panel drag-and-drop
  wireGlobalDragEvents(wasmMemory);

  // Wire JSON preview textarea
  const preview = document.getElementById('json-preview');
  if (preview) {
    dispatchChange('all'); // will be caught by onStateChange below
  }
}

// ── Config change handler ─────────────────────────────────────────────────────

function onConfigChange() {
  const cfgJson = JSON.stringify(buildConfigJson());
  pushConfig(cfgJson);

  // Update JSON preview if visible
  const preview = document.getElementById('json-preview');
  if (preview && !preview.closest('.json-section')?.classList.contains('collapsed')) {
    preview.value = JSON.stringify(buildConfigJson(), null, 2);
  }

  // After patch changes, redraw particle cards so new patch colors appear
  redrawAllCards();
}

// ── Cross-panel drag wiring ───────────────────────────────────────────────────

function wireGlobalDragEvents(wasmMemory) {
  // We need to intercept global pointermove/pointerup to handle drops that
  // cross panel boundaries.

  document.addEventListener('pointermove', e => {
    const drag = getActiveDrag();
    if (!drag) return;

    if (drag.type === 'patchType') {
      // Show drop preview on whichever particle card is under the cursor
      const pt = drag.payload;
      updatePatchDropPreview(e.clientX, e.clientY, pt.color);
    }
  });

  document.addEventListener('pointerup', e => {
    const drag = getActiveDrag();
    if (!drag) return;

    if (drag.type === 'patchType') {
      // Find which particle canvas, if any, received the drop
      const canvases = document.querySelectorAll('.particle-canvas');
      for (const canvas of canvases) {
        if (isOverElement(e.clientX, e.clientY, canvas)) {
          const typeIndex = parseInt(canvas.dataset.typeIndex, 10);
          receivePatchDrop(typeIndex, drag.payload.patchTypeId, e.clientX, e.clientY);
          break;
        }
      }
      // Clear all previews
      updatePatchDropPreview(-9999, -9999, '#fff');
    }

    if (drag.type === 'particleType') {
      const canvas = getPanel3Canvas();
      if (canvas && isOverElement(e.clientX, e.clientY, canvas)) {
        receiveParticleDrop(drag.payload.typeId, e.clientX, e.clientY);
        // Rebuild WASM sim with updated config so new particle type is recognised
        const cfgJson = JSON.stringify(buildConfigJson());
        _editorSim.update_config(cfgJson);
      }
    }
  });
}

// ── JSON preview toggle ───────────────────────────────────────────────────────

document.addEventListener('DOMContentLoaded', () => {
  const toggle = document.getElementById('json-toggle');
  const section = document.querySelector('.json-section');
  const preview = document.getElementById('json-preview');

  if (toggle && section) {
    toggle.addEventListener('click', () => {
      section.classList.toggle('collapsed');
      if (!section.classList.contains('collapsed') && preview) {
        preview.value = JSON.stringify(buildConfigJson(), null, 2);
      }
    });
  }

  main().catch(err => {
    console.error('Editor init failed:', err);
    document.body.innerHTML = `<pre style="color:red;padding:20px">Failed to load editor:\n${err}</pre>`;
  });
});
