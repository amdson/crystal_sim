/**
 * renderer2d.js — Shared canvas 2D drawing helpers.
 *
 * Used by both Panel 2 (particle type card previews) and Panel 3
 * (physics preview canvas). Keeps the circle+patch rendering logic
 * in one place instead of duplicating it across panels.
 *
 * Coordinate system: simulation units, y-axis pointing up (y is negated
 * when converting to canvas pixels, matching the existing renderer.js).
 */

// ── Core drawing primitives ──────────────────────────────────────────────────

/**
 * Draw a single particle (circle + patches) at canvas pixel coords (px, py).
 *
 * @param {CanvasRenderingContext2D} ctx
 * @param {number} px      Canvas x in pixels
 * @param {number} py      Canvas y in pixels
 * @param {number} pr      Radius in pixels
 * @param {number} oriRad  Orientation angle in radians
 * @param {object} typeDef {color, patches:[{angle_rad, color}]}
 * @param {boolean} [highlight] Draw a selection ring if true
 */
export function drawParticle(ctx, px, py, pr, oriRad, typeDef, highlight = false) {
  // Circle fill
  ctx.beginPath();
  ctx.arc(px, py, Math.max(pr, 1), 0, Math.PI * 2);
  ctx.fillStyle = typeDef.color ?? '#8888cc';
  ctx.fill();
  ctx.strokeStyle = highlight ? 'rgba(255,255,100,0.8)' : 'rgba(255,255,255,0.25)';
  ctx.lineWidth = highlight ? 2 : 0.5;
  ctx.stroke();

  // Patches — colored squares at the boundary
  const patches   = typeDef.patches ?? [];
  const patchHalf = Math.max(pr * 0.18, 2);
  for (const patch of patches) {
    const labAngle = oriRad + patch.angle_rad;
    const sx = px + Math.cos(labAngle) * pr;
    const sy = py - Math.sin(labAngle) * pr;  // y-axis flipped
    ctx.save();
    ctx.translate(sx, sy);
    ctx.rotate(-labAngle);
    ctx.fillStyle = patch.color ?? '#ffffff';
    ctx.fillRect(-patchHalf, -patchHalf, patchHalf * 2, patchHalf * 2);
    ctx.restore();
  }
}

/**
 * Draw a single patch square at a specific canvas position/angle.
 * Used for drag ghosts and hover previews.
 */
export function drawPatchSquare(ctx, px, py, size, color, angle = 0) {
  ctx.save();
  ctx.translate(px, py);
  ctx.rotate(angle);
  ctx.fillStyle = color;
  ctx.fillRect(-size / 2, -size / 2, size, size);
  ctx.restore();
}

/**
 * Draw a patch highlight ring (used when hovering over a patch in Panel 2).
 */
export function drawPatchHighlight(ctx, px, py, size) {
  ctx.beginPath();
  ctx.arc(px, py, size * 1.5, 0, Math.PI * 2);
  ctx.strokeStyle = 'rgba(255,255,100,0.9)';
  ctx.lineWidth = 1.5;
  ctx.stroke();
}

/**
 * Draw bond lines between particles in a flat stride-5 buffer.
 * Matches the heuristic in the existing renderer.js.
 */
export function drawBonds(ctx, buf, count, scale, cx, cy) {
  if (count < 2) return;
  ctx.strokeStyle = 'rgba(255,255,255,0.15)';
  ctx.lineWidth   = 1;
  ctx.beginPath();

  const S = 5;
  for (let i = 0; i < count; i++) {
    const bi = i * S;
    const xi = cx + buf[bi]     * scale;
    const yi = cy - buf[bi + 1] * scale;
    const ri = buf[bi + 3]      * scale;

    for (let j = i + 1; j < count; j++) {
      const bj = j * S;
      const xj = cx + buf[bj]     * scale;
      const yj = cy - buf[bj + 1] * scale;
      const rj = buf[bj + 3]      * scale;

      const dx = xi - xj, dy = yi - yj;
      const bondCutoff = (ri + rj) * 1.25;
      if (dx * dx + dy * dy < bondCutoff * bondCutoff) {
        ctx.moveTo(xi, yi);
        ctx.lineTo(xj, yj);
      }
    }
  }
  ctx.stroke();
}

// ── Full frame render for Panel 3 ────────────────────────────────────────────

/**
 * Render one frame of the physics preview (Panel 3).
 *
 * @param {CanvasRenderingContext2D} ctx
 * @param {number}       W          Canvas width in pixels
 * @param {number}       H          Canvas height in pixels
 * @param {Float32Array} buf        Stride-5 particle buffer from WASM
 * @param {number}       count      Number of particles
 * @param {Array}        typeMeta   Per-type metadata (same format as CrystalSim)
 * @param {object}       camera     {x, y, scale}  sim coords at canvas centre
 * @param {number|null}  hoveredIdx Particle index currently under pointer, or null
 */
export function renderPanel3Frame(ctx, W, H, buf, count, typeMeta, camera, hoveredIdx = null) {
  ctx.clearRect(0, 0, W, H);

  const cx = W / 2 - camera.x * camera.scale;
  const cy = H / 2 + camera.y * camera.scale;

  drawBonds(ctx, buf, count, camera.scale, cx, cy);

  const S = 5;
  for (let i = 0; i < count; i++) {
    const base   = i * S;
    const px     = cx + buf[base]     * camera.scale;
    const py     = cy - buf[base + 1] * camera.scale;
    const typeId = Math.round(buf[base + 2]);
    const pr     = buf[base + 3]      * camera.scale;
    const oriRad = buf[base + 4];
    const meta   = typeMeta[typeId] ?? { color: '#8888cc', patches: [] };
    drawParticle(ctx, px, py, pr, oriRad, meta, i === hoveredIdx);
  }
}

// ── Particle card preview for Panel 2 ────────────────────────────────────────

/**
 * Render a particle type card canvas.
 *
 * @param {CanvasRenderingContext2D} ctx
 * @param {number}  W            Canvas width
 * @param {number}  H            Canvas height
 * @param {object}  particleType {color, radius, patches:[{patchTypeId, position_deg}]}
 * @param {Array}   patchTypes   Global patch type list [{id, name, color}]
 * @param {number|null} hoveredPatchIdx   Which patch is hovered (for highlight)
 * @param {number|null} draggingPatchIdx  Which patch is being dragged
 * @param {number|null} previewAngleDeg   Hover-drop preview angle (null if not dragging)
 * @param {string|null} previewPatchColor Color for preview patch
 */
export function renderParticleCard(
  ctx, W, H,
  particleType, patchTypes,
  hoveredPatchIdx = null,
  draggingPatchIdx = null,
  previewAngleDeg = null,
  previewPatchColor = null,
) {
  ctx.clearRect(0, 0, W, H);

  const cx = W / 2;
  const cy = H / 2;
  const scale = Math.min(W, H) * 0.38;  // radius in pixels
  const pr = scale;

  // Build typeDef in the format drawParticle expects
  const typeDef = {
    color: particleType.color,
    patches: particleType.patches
      .filter((_, idx) => idx !== draggingPatchIdx)
      .map(patch => {
        const pt = patchTypes.find(p => p.id === patch.patchTypeId);
        return {
          angle_rad: patch.position_deg * Math.PI / 180,
          color: pt?.color ?? '#ffffff',
        };
      }),
  };

  drawParticle(ctx, cx, cy, pr, 0, typeDef);

  // Draw hovered patch highlight
  if (hoveredPatchIdx !== null && hoveredPatchIdx >= 0 && hoveredPatchIdx < particleType.patches.length) {
    const patch = particleType.patches[hoveredPatchIdx];
    const ang = patch.position_deg * Math.PI / 180;
    const sx = cx + Math.cos(ang) * pr;
    const sy = cy - Math.sin(ang) * pr;
    drawPatchHighlight(ctx, sx, sy, Math.max(pr * 0.18, 2));
  }

  // Draw drop preview patch
  if (previewAngleDeg !== null && previewPatchColor !== null) {
    const ang = previewAngleDeg * Math.PI / 180;
    const sx = cx + Math.cos(ang) * pr;
    const sy = cy - Math.sin(ang) * pr;
    const half = Math.max(pr * 0.18, 2);
    ctx.globalAlpha = 0.6;
    drawPatchSquare(ctx, sx, sy, half * 2, previewPatchColor, -ang);
    ctx.globalAlpha = 1.0;
  }
}

// ── Utilities ────────────────────────────────────────────────────────────────

/**
 * Hit-test: find index of the patch closest to (mx, my) on a particle card canvas.
 * Returns patch index or -1 if none is within the threshold.
 *
 * @param {number} mx    Mouse x in canvas pixels
 * @param {number} my    Mouse y in canvas pixels
 * @param {number} W     Canvas width
 * @param {number} H     Canvas height
 * @param {Array}  patches [{position_deg}]
 * @param {number} [threshold]  Hit radius in pixels (default 10)
 */
export function hitTestPatch(mx, my, W, H, patches, threshold = 10) {
  const cx = W / 2;
  const cy = H / 2;
  const pr = Math.min(W, H) * 0.38;

  let best = -1, bestDist = threshold * threshold;
  for (let i = 0; i < patches.length; i++) {
    const ang = patches[i].position_deg * Math.PI / 180;
    const sx = cx + Math.cos(ang) * pr;
    const sy = cy - Math.sin(ang) * pr;
    const d2 = (mx - sx) ** 2 + (my - sy) ** 2;
    if (d2 < bestDist) { bestDist = d2; best = i; }
  }
  return best;
}

/**
 * Convert canvas pixel (mx, my) relative to the card center to a position_deg
 * angle on the particle circle.
 */
export function pixelToAngleDeg(mx, my, W, H) {
  const cx = W / 2;
  const cy = H / 2;
  const dx = mx - cx;
  const dy = cy - my;  // y-axis flip
  return Math.atan2(dy, dx) * 180 / Math.PI;
}

/** Convert client coords to sim-space coords given a camera. */
export function clientToSim(clientX, clientY, canvas, camera) {
  const rect = canvas.getBoundingClientRect();
  const px = clientX - rect.left;
  const py = clientY - rect.top;
  const W = canvas.width;
  const H = canvas.height;
  const cx = W / 2 - camera.x * camera.scale;
  const cy = H / 2 + camera.y * camera.scale;
  return {
    x:  (px - cx) / camera.scale,
    y: -(py - cy) / camera.scale,  // y-axis flip
  };
}

/**
 * Hit-test Panel 3 canvas: find index of particle closest to sim coords (sx, sy).
 * Returns {index, dist} or {index: -1}.
 */
export function hitTestParticle(sx, sy, buf, count) {
  let best = -1, bestDist = Infinity;
  const S = 5;
  for (let i = 0; i < count; i++) {
    const bx = buf[i * S];
    const by = buf[i * S + 1];
    const r  = buf[i * S + 3];
    const d  = Math.sqrt((sx - bx) ** 2 + (sy - by) ** 2);
    if (d < r * 1.5 && d < bestDist) { bestDist = d; best = i; }
  }
  return { index: best, dist: bestDist };
}

/** Hex color #rrggbb → rgba() string */
export function hexToRgba(hex, alpha) {
  const h = hex.replace('#', '');
  const r = parseInt(h.slice(0, 2), 16);
  const g = parseInt(h.slice(2, 4), 16);
  const b = parseInt(h.slice(4, 6), 16);
  return `rgba(${r},${g},${b},${alpha.toFixed(2)})`;
}
