/**
 * panel2.js — Particle type designer.
 *
 * Shows a scrollable list of "particle type cards". Each card has:
 *   - A canvas preview of the particle circle with its patches
 *   - Color picker, radius input, mu input
 *   - Drag-and-drop: receive patches from Panel 1, rotate/delete patches on canvas
 *   - Drag source for Panel 3 particle placement
 */

import { state, dispatchChange } from './state.js';
import {
  renderParticleCard, hitTestPatch, pixelToAngleDeg,
} from './renderer2d.js';
import { startDrag, isOverElement } from './drag.js';

const CARD_CANVAS_SIZE = 120;

let _container      = null;
let _onConfigChange = null;

export function initPanel2(container, onConfigChange) {
  _container      = container;
  _onConfigChange = onConfigChange;
  render();
}

export function refreshPanel2() {
  render();
}

// ── Top-level render ──────────────────────────────────────────────────────────

function render() {
  if (!_container) return;
  _container.innerHTML = '';

  const title = document.createElement('h2');
  title.className = 'panel-title';
  title.textContent = 'Particles';
  _container.appendChild(title);

  const addBtn = document.createElement('button');
  addBtn.className = 'btn btn-add';
  addBtn.textContent = '+ Particle Type';
  addBtn.addEventListener('click', addParticleType);
  _container.appendChild(addBtn);

  for (let i = 0; i < state.particleTypes.length; i++) {
    _container.appendChild(buildCard(i));
  }
}

// ── Particle type card ────────────────────────────────────────────────────────

function buildCard(typeIndex) {
  const pt = state.particleTypes[typeIndex];

  const card = document.createElement('div');
  card.className = 'particle-card';
  card.dataset.typeIndex = typeIndex;

  // ── Card header ──────────────────────────────────────────────────────────
  const header = document.createElement('div');
  header.className = 'card-header';

  const colorPicker = document.createElement('input');
  colorPicker.type  = 'color';
  colorPicker.value = pt.color;
  colorPicker.title = 'Particle color';
  colorPicker.addEventListener('input', () => {
    pt.color = colorPicker.value;
    redrawCardCanvas(typeIndex);
    dispatchChange('particles');
    _onConfigChange?.();
  });

  const nameLabel = document.createElement('span');
  nameLabel.className = 'card-name';
  nameLabel.textContent = pt.name;

  const delBtn = document.createElement('button');
  delBtn.className = 'card-delete';
  delBtn.textContent = '×';
  delBtn.title = 'Delete particle type';
  delBtn.addEventListener('click', () => {
    state.particleTypes.splice(typeIndex, 1);
    // Re-index ids
    state.particleTypes.forEach((p, i) => p.id = i);
    // Remove instances of this type
    state.instances = state.instances.filter(inst => inst.typeId !== typeIndex);
    dispatchChange('particles');
    _onConfigChange?.();
    render();
  });

  // Drag handle for dropping onto Panel 3
  header.style.cursor = 'grab';
  header.addEventListener('pointerdown', e => {
    e.preventDefault();
    startDrag({
      sourceEl: header,
      pointerId: e.pointerId,
      type: 'particleType',
      payload: { typeId: typeIndex },
      ghostColor: pt.color,
      ghostLabel: pt.name.slice(0, 3),
    });
  });

  header.appendChild(colorPicker);
  header.appendChild(nameLabel);
  header.appendChild(delBtn);
  card.appendChild(header);

  // ── Canvas preview ────────────────────────────────────────────────────────
  const canvas = document.createElement('canvas');
  canvas.width  = CARD_CANVAS_SIZE;
  canvas.height = CARD_CANVAS_SIZE;
  canvas.className = 'particle-canvas';
  canvas.dataset.typeIndex = typeIndex;
  card.appendChild(canvas);

  drawCardCanvas(canvas, typeIndex);
  wireCardCanvas(canvas, typeIndex, card);

  // ── Input row ────────────────────────────────────────────────────────────
  const inputs = document.createElement('div');
  inputs.className = 'card-inputs';

  inputs.appendChild(makeNumberInput('r', 'Radius', pt.radius, 0.1, 0.1, v => {
    pt.radius = v;
    redrawCardCanvas(typeIndex);
    dispatchChange('particles');
    _onConfigChange?.();
  }));

  inputs.appendChild(makeNumberInput('μ', 'Chemical potential (mu)', pt.mu, 0.5, -Infinity, v => {
    pt.mu = v;
    dispatchChange('particles');
    _onConfigChange?.();
  }));

  card.appendChild(inputs);
  return card;
}

function makeNumberInput(label, title, value, step, min, onChange) {
  const wrap = document.createElement('label');
  wrap.className = 'input-wrap';
  wrap.title = title;
  const lbl = document.createElement('span');
  lbl.textContent = label;
  const inp = document.createElement('input');
  inp.type = 'number';
  inp.step = step;
  if (min !== -Infinity) inp.min = min;
  inp.value = value;
  inp.addEventListener('change', () => onChange(parseFloat(inp.value) || 0));
  wrap.appendChild(lbl);
  wrap.appendChild(inp);
  return wrap;
}

// ── Canvas drawing ────────────────────────────────────────────────────────────

function drawCardCanvas(canvas, typeIndex, opts = {}) {
  const ctx = canvas.getContext('2d');
  const pt  = state.particleTypes[typeIndex];
  renderParticleCard(
    ctx,
    canvas.width,
    canvas.height,
    pt,
    state.patchTypes,
    opts.hoveredPatchIdx   ?? null,
    opts.draggingPatchIdx  ?? null,
    opts.previewAngleDeg   ?? null,
    opts.previewPatchColor ?? null,
  );
}

function redrawCardCanvas(typeIndex) {
  const canvas = document.querySelector(`.particle-canvas[data-type-index="${typeIndex}"]`);
  if (canvas) drawCardCanvas(canvas, typeIndex);
}

/** Redraw all card canvases (called after patch types change). */
export function redrawAllCards() {
  for (let i = 0; i < state.particleTypes.length; i++) {
    redrawCardCanvas(i);
  }
}

// ── Canvas pointer interactions ───────────────────────────────────────────────

function wireCardCanvas(canvas, typeIndex, card) {
  let hoveredPatchIdx = null;

  // ── Hover: highlight nearest patch ────────────────────────────────────────
  canvas.addEventListener('pointermove', e => {
    if (e.buttons !== 0) return; // don't change highlight while dragging
    const rect = canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    const pt = state.particleTypes[typeIndex];
    const idx = hitTestPatch(mx, my, canvas.width, canvas.height, pt.patches, 12);
    if (idx !== hoveredPatchIdx) {
      hoveredPatchIdx = idx;
      drawCardCanvas(canvas, typeIndex, { hoveredPatchIdx });
    }
  });

  canvas.addEventListener('pointerleave', () => {
    if (hoveredPatchIdx !== null) {
      hoveredPatchIdx = null;
      drawCardCanvas(canvas, typeIndex);
    }
  });

  // ── Pointerdown: start patch drag or canvas drag ───────────────────────────
  canvas.addEventListener('pointerdown', e => {
    e.preventDefault();
    const rect = canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    const pt = state.particleTypes[typeIndex];
    const patchIdx = hitTestPatch(mx, my, canvas.width, canvas.height, pt.patches, 14);

    if (patchIdx !== -1) {
      // Start patch-on-particle drag
      const patch = pt.patches[patchIdx];
      const patchType = state.patchTypes.find(p => p.id === patch.patchTypeId);
      startPatchDrag(canvas, typeIndex, patchIdx, patchType, e.pointerId, card);
    }
    // No drag started for canvas background — pointer events pass through to parent if needed
  });

  // ── Drop target: receive patch from Panel 1 ────────────────────────────────
  // We handle this in the global pointerup listener wired in main.js,
  // but we add a mouseover highlight here.
  canvas.addEventListener('dragover', e => e.preventDefault()); // not used but safe
}

// ── Patch drag on particle ────────────────────────────────────────────────────

function startPatchDrag(canvas, typeIndex, patchIdx, patchType, pointerId, card) {
  const pt = state.particleTypes[typeIndex];

  startDrag({
    sourceEl: canvas,
    pointerId,
    type: 'patchOnParticle',
    payload: { typeIndex, patchIdx },
    ghostColor: patchType?.color ?? '#ffffff',
    ghostLabel: patchType?.name ?? '?',

    onMove: (clientX, clientY) => {
      // Update the preview angle while dragging over the same canvas
      if (isOverElement(clientX, clientY, canvas)) {
        const rect = canvas.getBoundingClientRect();
        const mx = clientX - rect.left;
        const my = clientY - rect.top;
        const deg = pixelToAngleDeg(mx, my, canvas.width, canvas.height);
        drawCardCanvas(canvas, typeIndex, {
          draggingPatchIdx: patchIdx,
          previewAngleDeg: deg,
          previewPatchColor: patchType?.color ?? '#fff',
        });
      } else {
        // Show as being dragged off — indicate deletion
        drawCardCanvas(canvas, typeIndex, { draggingPatchIdx: patchIdx });
      }
    },

    onDrop: (clientX, clientY) => {
      if (isOverElement(clientX, clientY, canvas)) {
        // Commit new angle
        const rect = canvas.getBoundingClientRect();
        const mx = clientX - rect.left;
        const my = clientY - rect.top;
        pt.patches[patchIdx].position_deg =
          pixelToAngleDeg(mx, my, canvas.width, canvas.height);
      } else {
        // Dragged off — delete patch
        pt.patches.splice(patchIdx, 1);
      }
      drawCardCanvas(canvas, typeIndex);
      dispatchChange('particles');
      _onConfigChange?.();
    },

    onCancel: () => drawCardCanvas(canvas, typeIndex),
  });
}

// ── Receive patch drop from Panel 1 ──────────────────────────────────────────

/**
 * Called from main.js when a "patchType" drag is dropped over a particle canvas.
 *
 * @param {number} typeIndex      Target particle type index
 * @param {number} patchTypeId    The patch type being dropped
 * @param {number} clientX/Y      Drop position in client coordinates
 */
export function receivePatchDrop(typeIndex, patchTypeId, clientX, clientY) {
  const canvas = document.querySelector(`.particle-canvas[data-type-index="${typeIndex}"]`);
  if (!canvas) return;

  const rect   = canvas.getBoundingClientRect();
  const mx     = clientX - rect.left;
  const my     = clientY - rect.top;
  const deg    = pixelToAngleDeg(mx, my, canvas.width, canvas.height);

  state.particleTypes[typeIndex].patches.push({ patchTypeId, position_deg: deg });
  drawCardCanvas(canvas, typeIndex);
  dispatchChange('particles');
  _onConfigChange?.();
}

/**
 * Called from main.js on pointermove during a "patchType" drag
 * to show a drop preview on the card canvas under the cursor.
 */
export function updatePatchDropPreview(clientX, clientY, patchColor) {
  // Clear all previews first
  for (let i = 0; i < state.particleTypes.length; i++) {
    const canvas = document.querySelector(`.particle-canvas[data-type-index="${i}"]`);
    if (!canvas) continue;
    if (isOverElement(clientX, clientY, canvas)) {
      const rect = canvas.getBoundingClientRect();
      const mx = clientX - rect.left;
      const my = clientY - rect.top;
      const deg = pixelToAngleDeg(mx, my, canvas.width, canvas.height);
      drawCardCanvas(canvas, i, { previewAngleDeg: deg, previewPatchColor: patchColor });
    } else {
      drawCardCanvas(canvas, i);
    }
  }
}

// ── Add particle type ─────────────────────────────────────────────────────────

function addParticleType() {
  const id = state.particleTypes.length;
  const colors = ['#4a90d9', '#e84040', '#f5a623', '#7ed321', '#bd10e0'];
  state.particleTypes.push({
    id,
    name: `Type ${id}`,
    color: colors[id % colors.length],
    radius: 1.0,
    mu: -2.0,
    patches: [],
  });
  dispatchChange('particles');
  _onConfigChange?.();
  render();
}
