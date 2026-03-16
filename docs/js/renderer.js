/**
 * renderer.js — Canvas 2D render loop.
 *
 * Reads the flat particle buffer from WASM memory on every frame and draws
 * circles (and optionally bond lines) without any copies.
 *
 * Coordinate mapping: simulation units → canvas pixels via a camera object
 * { x, y, scale } where (x,y) is the simulation point at the canvas centre.
 */

export class Renderer {
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx    = canvas.getContext('2d');
    this.camera = { x: 0, y: 0, scale: 30 }; // 30 px per sim unit initially
    this.showBonds = true;
    this.showCandidates = false;

    this._setupResize();
    this._setupPanZoom();
  }

  /** Draw one frame given the particle buffer, type metadata, and optional candidate buffer. */
  draw(buf, count, typeMeta, candBuf, candCount) {
    const { ctx, canvas, camera } = this;
    const W = canvas.width;
    const H = canvas.height;

    ctx.clearRect(0, 0, W, H);

    const cx = W / 2 - camera.x * camera.scale;
    const cy = H / 2 + camera.y * camera.scale;  // y-axis flipped

    // ── Candidate sites ─────────────────────────────────────────────────────
    if (this.showCandidates && candBuf && candCount > 0) {
      const CSTRIDE = 4;
      // Find max rate for normalisation.
      let maxRate = 0;
      for (let i = 0; i < candCount; i++) {
        const r = candBuf[i * CSTRIDE + 3];
        if (r > maxRate) maxRate = r;
      }
      if (maxRate === 0) maxRate = 1;

      for (let i = 0; i < candCount; i++) {
        const base   = i * CSTRIDE;
        const sx     = cx + candBuf[base]     * camera.scale;
        const sy     = cy - candBuf[base + 1] * camera.scale;
        const typeId = Math.round(candBuf[base + 2]);
        const rate   = candBuf[base + 3];
        const t      = rate / maxRate;   // 0..1
        const dotR   = 1.5 + t * 4.5;   // 1.5..6 px

        const meta  = typeMeta[typeId] ?? { color: '#8888cc' };
        ctx.beginPath();
        ctx.arc(sx, sy, dotR, 0, Math.PI * 2);
        ctx.fillStyle = meta.color.replace(/^#/, '') === meta.color
          ? `rgba(255,200,60,${0.25 + t * 0.55})`
          : hexToRgba(meta.color, 0.25 + t * 0.55);
        ctx.fill();
      }
    }

    if (count === 0) return;

    const STRIDE = 5;

    // ── Bond lines ─────────────────────────────────────────────────────────
    if (this.showBonds) {
      ctx.strokeStyle = 'rgba(255,255,255,0.12)';
      ctx.lineWidth   = 1;
      ctx.beginPath();

      for (let i = 0; i < count; i++) {
        const bi  = i * STRIDE;
        const xi  = cx + buf[bi]     * camera.scale;
        const yi  = cy - buf[bi + 1] * camera.scale;
        const ri  = buf[bi + 3]      * camera.scale;

        for (let j = i + 1; j < count; j++) {
          const bj  = j * STRIDE;
          const xj  = cx + buf[bj]     * camera.scale;
          const yj  = cy - buf[bj + 1] * camera.scale;
          const rj  = buf[bj + 3]      * camera.scale;

          const dx = xi - xj;
          const dy = yi - yj;
          const d2 = dx * dx + dy * dy;
          const bondCutoff = (ri + rj) * 1.25; // ~contact + some slack
          if (d2 < bondCutoff * bondCutoff) {
            ctx.moveTo(xi, yi);
            ctx.lineTo(xj, yj);
          }
        }
      }
      ctx.stroke();
    }

    // ── Particle circles + patches ───────────────────────────────────────────
    for (let i = 0; i < count; i++) {
      const base        = i * STRIDE;
      const px          = cx + buf[base]     * camera.scale;
      const py          = cy - buf[base + 1] * camera.scale;
      const typeId      = Math.round(buf[base + 2]);
      const pr          = buf[base + 3]      * camera.scale;
      const orientation = buf[base + 4];

      const meta  = typeMeta[typeId] ?? { color: '#8888cc', patches: [] };
      const color = meta.color;

      // Circle
      ctx.beginPath();
      ctx.arc(px, py, Math.max(pr, 1), 0, Math.PI * 2);
      ctx.fillStyle   = color;
      ctx.fill();
      ctx.strokeStyle = 'rgba(255,255,255,0.25)';
      ctx.lineWidth   = 0.5;
      ctx.stroke();

      // Patches — squares at the particle boundary
      const patches   = meta.patches ?? [];
      const patchHalf = Math.max(pr * 0.18, 1.5);
      for (const patch of patches) {
        const labAngle = orientation + patch.angle_rad;
        const sx = px + Math.cos(labAngle)  * pr;
        const sy = py - Math.sin(labAngle)  * pr;  // y-axis flipped
        ctx.save();
        ctx.translate(sx, sy);
        ctx.rotate(-labAngle);
        ctx.fillStyle = patch.color ?? '#ffffff';
        ctx.fillRect(-patchHalf, -patchHalf, patchHalf * 2, patchHalf * 2);
        ctx.restore();
      }
    }
  }

  // ── Resize handling ───────────────────────────────────────────────────────

  _setupResize() {
    const resize = () => {
      this.canvas.width  = this.canvas.offsetWidth;
      this.canvas.height = this.canvas.offsetHeight;
    };
    window.addEventListener('resize', resize);
    resize();
  }

  // ── Pan + zoom ────────────────────────────────────────────────────────────

  _setupPanZoom() {
    const c = this.canvas;
    let dragging = false;
    let lastX = 0, lastY = 0;

    c.addEventListener('mousedown', e => {
      dragging = true;
      lastX = e.clientX;
      lastY = e.clientY;
    });
    window.addEventListener('mouseup', () => { dragging = false; });
    window.addEventListener('mousemove', e => {
      if (!dragging) return;
      const dx = e.clientX - lastX;
      const dy = e.clientY - lastY;
      this.camera.x -= dx / this.camera.scale;
      this.camera.y += dy / this.camera.scale;
      lastX = e.clientX;
      lastY = e.clientY;
    });

    c.addEventListener('wheel', e => {
      e.preventDefault();
      const factor = e.deltaY < 0 ? 1.1 : 0.9;
      this.camera.scale = Math.min(Math.max(this.camera.scale * factor, 2), 500);
    }, { passive: false });
  }
}

/** Convert a #rrggbb hex color to an rgba() string with the given alpha. */
function hexToRgba(hex, alpha) {
  const h = hex.replace('#', '');
  const r = parseInt(h.slice(0, 2), 16);
  const g = parseInt(h.slice(2, 4), 16);
  const b = parseInt(h.slice(4, 6), 16);
  return `rgba(${r},${g},${b},${alpha.toFixed(2)})`;
}
