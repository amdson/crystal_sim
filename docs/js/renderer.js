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

    this._setupResize();
    this._setupPanZoom();
  }

  /** Draw one frame given the particle buffer and type metadata. */
  draw(buf, count, typeMeta) {
    const { ctx, canvas, camera } = this;
    const W = canvas.width;
    const H = canvas.height;

    ctx.clearRect(0, 0, W, H);

    if (count === 0) return;

    const cx = W / 2 - camera.x * camera.scale;
    const cy = H / 2 + camera.y * camera.scale;  // y-axis flipped

    // ── Bond lines ─────────────────────────────────────────────────────────
    if (this.showBonds) {
      ctx.strokeStyle = 'rgba(255,255,255,0.12)';
      ctx.lineWidth   = 1;
      ctx.beginPath();

      for (let i = 0; i < count; i++) {
        const bi  = i * 4;
        const xi  = cx + buf[bi]     * camera.scale;
        const yi  = cy - buf[bi + 1] * camera.scale;
        const ri  = buf[bi + 3]      * camera.scale;

        for (let j = i + 1; j < count; j++) {
          const bj  = j * 4;
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

    // ── Particle circles ────────────────────────────────────────────────────
    for (let i = 0; i < count; i++) {
      const base   = i * 4;
      const px     = cx + buf[base]     * camera.scale;
      const py     = cy - buf[base + 1] * camera.scale;
      const typeId = Math.round(buf[base + 2]);
      const pr     = buf[base + 3]      * camera.scale;

      const meta  = typeMeta[typeId] ?? { color: '#8888cc' };
      const color = meta.color;

      ctx.beginPath();
      ctx.arc(px, py, Math.max(pr, 1), 0, Math.PI * 2);
      ctx.fillStyle   = color;
      ctx.fill();
      ctx.strokeStyle = 'rgba(255,255,255,0.25)';
      ctx.lineWidth   = 0.5;
      ctx.stroke();
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
