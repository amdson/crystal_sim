/**
 * drag.js — Pointer-event drag orchestration.
 *
 * All drag interactions use pointerdown/pointermove/pointerup + setPointerCapture
 * so the browser doesn't interfere with custom rendering.
 *
 * Drag types:
 *   "patchType"       — dragging a patch token from Panel 1 onto a Panel 2 card
 *   "patchOnParticle" — dragging a patch within a Panel 2 card canvas
 *   "particleType"    — dragging a particle type card onto Panel 3
 *   "instanceMove"    — dragging a Panel 3 particle instance
 */

// ── Ghost element ────────────────────────────────────────────────────────────

let _ghost = null;

function createGhost(color, label = '') {
  if (_ghost) _ghost.remove();
  _ghost = document.createElement('div');
  _ghost.className = 'drag-ghost';
  _ghost.style.cssText = `
    position: fixed; pointer-events: none; z-index: 9999;
    width: 28px; height: 28px; border-radius: 50%;
    background: ${color}; border: 2px solid rgba(255,255,255,0.6);
    display: flex; align-items: center; justify-content: center;
    font-size: 11px; font-weight: bold; color: #fff;
    text-shadow: 0 0 3px #000;
    transform: translate(-50%, -50%);
  `;
  _ghost.textContent = label;
  document.body.appendChild(_ghost);
  return _ghost;
}

function moveGhost(clientX, clientY) {
  if (_ghost) {
    _ghost.style.left = clientX + 'px';
    _ghost.style.top  = clientY + 'px';
  }
}

function removeGhost() {
  if (_ghost) { _ghost.remove(); _ghost = null; }
}

// ── Active drag state ────────────────────────────────────────────────────────

let _active = null;

/**
 * Start a drag operation.
 *
 * @param {Element}  sourceEl   The element the drag started on (for pointer capture)
 * @param {number}   pointerId  Pointer ID from the pointerdown event
 * @param {string}   type       Drag type discriminant
 * @param {object}   payload    Type-specific data
 * @param {string}   ghostColor Ghost circle background color
 * @param {string}   ghostLabel Ghost label text (short)
 * @param {Function} onMove     (clientX, clientY) → void  — called on each move
 * @param {Function} onDrop     (clientX, clientY) → void  — called on pointer-up
 * @param {Function} onCancel   () → void                  — called on pointer-cancel
 */
export function startDrag({
  sourceEl, pointerId, type, payload,
  ghostColor = '#8888cc', ghostLabel = '',
  onMove = null, onDrop = null, onCancel = null,
}) {
  if (_active) cancelDrag();

  sourceEl.setPointerCapture(pointerId);

  createGhost(ghostColor, ghostLabel);

  _active = { type, payload, sourceEl, pointerId, onMove, onDrop, onCancel };

  const handleMove = e => {
    if (e.pointerId !== pointerId) return;
    moveGhost(e.clientX, e.clientY);
    onMove?.(e.clientX, e.clientY);
  };

  const handleUp = e => {
    if (e.pointerId !== pointerId) return;
    cleanup();
    onDrop?.(e.clientX, e.clientY);
    // Delay clearing _active so the pointerup event can finish bubbling to
    // document-level listeners (e.g. main.js drop handlers) before it goes null.
    queueMicrotask(() => { _active = null; });
  };

  const handleCancel = e => {
    if (e.pointerId !== pointerId) return;
    cleanup();
    onCancel?.();
    queueMicrotask(() => { _active = null; });
  };

  function cleanup() {
    removeGhost();
    sourceEl.removeEventListener('pointermove', handleMove);
    sourceEl.removeEventListener('pointerup',   handleUp);
    sourceEl.removeEventListener('pointercancel', handleCancel);
  }

  sourceEl.addEventListener('pointermove',   handleMove);
  sourceEl.addEventListener('pointerup',     handleUp);
  sourceEl.addEventListener('pointercancel', handleCancel);
}

function cancelDrag() {
  if (!_active) return;
  _active.onCancel?.();
  removeGhost();
  _active = null;
}

/** Returns true if a drag of the given type is in progress. */
export function isDragging(type = null) {
  if (!_active) return false;
  return type === null || _active.type === type;
}

export function getActiveDrag() {
  return _active;
}

// ── Drop zone helpers ─────────────────────────────────────────────────────────

/**
 * Check whether a client coordinate is inside an element's bounding rect.
 */
export function isOverElement(clientX, clientY, el) {
  const r = el.getBoundingClientRect();
  return clientX >= r.left && clientX <= r.right &&
         clientY >= r.top  && clientY <= r.bottom;
}
