/**
 * panel1.js — Patch palette and KxK interaction table.
 *
 * Responsibilities:
 *   - Render the row of patch type tokens
 *   - "New Patch" button
 *   - KxK symmetric interaction table with inline editing
 *   - Notify the rest of the app (via onConfigChange callback) whenever the
 *     patch definitions or interaction parameters change
 */

import {
  state, dispatchChange,
  nextPatchName, patchColorFor,
  interactionKey, setInteraction, getInteraction,
} from './state.js';
import { startDrag } from './drag.js';

let _container = null;
let _onConfigChange = null;

/**
 * Initialize Panel 1.
 *
 * @param {Element}  container     The #panel1 DOM element
 * @param {Function} onConfigChange Called whenever patches/interactions change
 */
export function initPanel1(container, onConfigChange) {
  _container = container;
  _onConfigChange = onConfigChange;
  render();
}

export function refreshPanel1() {
  render();
}

// ── Rendering ────────────────────────────────────────────────────────────────

function render() {
  if (!_container) return;
  _container.innerHTML = '';

  // Section title
  const title = document.createElement('h2');
  title.className = 'panel-title';
  title.textContent = 'Patches';
  _container.appendChild(title);

  // Palette row
  _container.appendChild(buildPalette());

  // Interaction table
  if (state.patchTypes.length > 0) {
    const tableSection = document.createElement('div');
    tableSection.className = 'interaction-section';
    const tableTitle = document.createElement('h3');
    tableTitle.className = 'section-subtitle';
    tableTitle.textContent = 'Patch Interactions';
    tableSection.appendChild(tableTitle);
    tableSection.appendChild(buildInteractionTable());
    _container.appendChild(tableSection);
  }
}

// ── Palette row ───────────────────────────────────────────────────────────────

function buildPalette() {
  const row = document.createElement('div');
  row.className = 'patch-palette';

  for (const pt of state.patchTypes) {
    row.appendChild(buildPatchToken(pt));
  }

  // New Patch button
  const addBtn = document.createElement('button');
  addBtn.className = 'btn btn-add';
  addBtn.textContent = '+ Patch';
  addBtn.addEventListener('click', addNewPatch);
  row.appendChild(addBtn);

  return row;
}

function buildPatchToken(pt) {
  const token = document.createElement('div');
  token.className = 'patch-token';
  token.dataset.patchId = pt.id;

  // Colored swatch circle
  const swatch = document.createElement('div');
  swatch.className = 'patch-swatch';
  swatch.style.background = pt.color;

  // Editable name label
  const label = document.createElement('span');
  label.className = 'patch-name';
  label.textContent = pt.name;
  label.contentEditable = 'true';
  label.spellcheck = false;
  label.addEventListener('blur', () => {
    const newName = label.textContent.trim();
    if (newName && newName !== pt.name) {
      pt.name = newName;
      dispatchChange('patches');
      _onConfigChange?.();
    }
  });
  label.addEventListener('keydown', e => {
    if (e.key === 'Enter') { e.preventDefault(); label.blur(); }
  });

  // Delete button
  const del = document.createElement('button');
  del.className = 'patch-delete';
  del.textContent = '×';
  del.title = 'Delete patch type';
  del.addEventListener('click', () => {
    // Only allow deletion if no particle type uses this patch.
    const inUse = state.particleTypes.some(t =>
      t.patches.some(p => p.patchTypeId === pt.id)
    );
    if (inUse) {
      alert(`Patch "${pt.name}" is used by one or more particle types. Remove it from all particles first.`);
      return;
    }
    state.patchTypes = state.patchTypes.filter(p => p.id !== pt.id);
    // Also remove interactions involving this patch.
    for (const [key] of state.interactions) {
      const parts = key.split(':');
      if (parts.includes(pt.name)) state.interactions.delete(key);
    }
    dispatchChange('patches');
    _onConfigChange?.();
    render();
  });

  token.appendChild(swatch);
  token.appendChild(label);
  token.appendChild(del);

  // Make patch token draggable (for dropping onto particle cards)
  token.addEventListener('pointerdown', e => {
    e.preventDefault();
    startDrag({
      sourceEl: token,
      pointerId: e.pointerId,
      type: 'patchType',
      payload: { patchTypeId: pt.id, patchName: pt.name, color: pt.color },
      ghostColor: pt.color,
      ghostLabel: pt.name,
    });
  });

  return token;
}

function addNewPatch() {
  const name  = nextPatchName();
  const color = patchColorFor(state.patchTypes.length);
  const id    = state.patchTypes.length === 0
    ? 0
    : Math.max(...state.patchTypes.map(p => p.id)) + 1;

  state.patchTypes.push({ id, name, color });

  // Initialize zero-epsilon interactions with all existing patch types.
  for (const existing of state.patchTypes) {
    if (existing.id === id) continue;
    const key = interactionKey(name, existing.name);
    if (!state.interactions.has(key)) {
      state.interactions.set(key, {
        epsilon: 0,
        angular_width_deg: [30, 30],
        cutoff: 1.35,
      });
    }
  }
  // Self-interaction
  const selfKey = interactionKey(name, name);
  if (!state.interactions.has(selfKey)) {
    state.interactions.set(selfKey, {
      epsilon: 0,
      angular_width_deg: [30, 30],
      cutoff: 1.35,
    });
  }

  dispatchChange('patches');
  _onConfigChange?.();
  render();
}

// ── KxK Interaction table ────────────────────────────────────────────────────

function buildInteractionTable() {
  const pts = state.patchTypes;
  const K = pts.length;

  const wrapper = document.createElement('div');
  wrapper.className = 'interaction-table-wrapper';

  const table = document.createElement('table');
  table.className = 'interaction-table';

  // Header row
  const thead = document.createElement('thead');
  const headRow = document.createElement('tr');
  headRow.appendChild(document.createElement('th')); // corner
  for (const pt of pts) {
    const th = document.createElement('th');
    th.style.color = pt.color;
    th.textContent = pt.name;
    headRow.appendChild(th);
  }
  thead.appendChild(headRow);
  table.appendChild(thead);

  // Body — upper triangle is editable, lower mirrors it
  const tbody = document.createElement('tbody');
  for (let i = 0; i < K; i++) {
    const tr = document.createElement('tr');

    // Row header
    const rowHead = document.createElement('th');
    rowHead.style.color = pts[i].color;
    rowHead.textContent = pts[i].name;
    tr.appendChild(rowHead);

    for (let j = 0; j < K; j++) {
      const td = document.createElement('td');
      const key = interactionKey(pts[i].name, pts[j].name);
      const int = state.interactions.get(key) ?? { epsilon: 0, angular_width_deg: [30, 30], cutoff: 1.35 };

      if (j < i) {
        // Lower triangle: mirror display
        td.className = 'interaction-cell mirror';
        td.textContent = int.epsilon.toFixed(1);
        td.style.opacity = '0.5';
      } else {
        // Upper triangle (including diagonal): editable
        td.className = 'interaction-cell editable';
        td.dataset.key = key;
        td.title = 'Click to edit';
        td.textContent = int.epsilon.toFixed(1);
        td.style.background = epsilonColor(int.epsilon);

        td.addEventListener('click', () => openInteractionPopover(td, pts[i].name, pts[j].name, key));
      }
      tr.appendChild(td);
    }
    tbody.appendChild(tr);
  }
  table.appendChild(tbody);
  wrapper.appendChild(table);
  return wrapper;
}

function epsilonColor(eps) {
  if (eps < 0)  return `rgba(74,144,217,${Math.min(0.6, Math.abs(eps) / 10)})`;
  if (eps > 0)  return `rgba(232,64,64,${Math.min(0.6, eps / 10)})`;
  return 'transparent';
}

// ── Inline popover for editing an interaction ─────────────────────────────────

let _activePopover = null;

function openInteractionPopover(cell, nameA, nameB, key) {
  closePopover();

  const int = state.interactions.get(key) ?? { epsilon: 0, angular_width_deg: [30, 30], cutoff: 1.35 };

  const pop = document.createElement('div');
  pop.className = 'interaction-popover';

  const title = document.createElement('div');
  title.className = 'popover-title';
  title.textContent = `${nameA} ↔ ${nameB}`;
  pop.appendChild(title);

  const fields = [
    { label: 'ε (epsilon)',     key: 'epsilon',              value: int.epsilon,                step: 0.5 },
    { label: 'σ₀ width (deg)', key: 'angular_width_deg[0]', value: int.angular_width_deg[0],   step: 5   },
    { label: 'σ₁ width (deg)', key: 'angular_width_deg[1]', value: int.angular_width_deg[1],   step: 5   },
    { label: 'Cutoff (×r_c)',  key: 'cutoff',               value: int.cutoff,                 step: 0.05 },
  ];

  const inputs = {};
  for (const f of fields) {
    const row = document.createElement('div');
    row.className = 'popover-row';
    const lbl = document.createElement('label');
    lbl.textContent = f.label;
    const inp = document.createElement('input');
    inp.type = 'number';
    inp.step = f.step;
    inp.value = f.value;
    inputs[f.key] = inp;
    row.appendChild(lbl);
    row.appendChild(inp);
    pop.appendChild(row);
  }

  const applyBtn = document.createElement('button');
  applyBtn.className = 'btn btn-apply';
  applyBtn.textContent = 'Apply';
  applyBtn.addEventListener('click', () => {
    const updated = {
      epsilon: parseFloat(inputs['epsilon'].value) || 0,
      angular_width_deg: [
        parseFloat(inputs['angular_width_deg[0]'].value) || 30,
        parseFloat(inputs['angular_width_deg[1]'].value) || 30,
      ],
      cutoff: parseFloat(inputs['cutoff'].value) || 1.35,
    };
    state.interactions.set(key, updated);
    closePopover();
    dispatchChange('patches');
    _onConfigChange?.();
    render();
  });
  pop.appendChild(applyBtn);

  // Position below the cell
  const rect = cell.getBoundingClientRect();
  pop.style.position = 'fixed';
  pop.style.left = rect.left + 'px';
  pop.style.top  = (rect.bottom + 4) + 'px';
  pop.style.zIndex = '1000';

  document.body.appendChild(pop);
  _activePopover = pop;

  // Close on outside click
  setTimeout(() => {
    document.addEventListener('pointerdown', closeOnOutside, { once: true });
  }, 10);
}

function closeOnOutside(e) {
  if (_activePopover && !_activePopover.contains(e.target)) closePopover();
}

function closePopover() {
  if (_activePopover) { _activePopover.remove(); _activePopover = null; }
}
