/**
 * controls.js — Wires DOM inputs to WASM parameter updates.
 *
 * Call `initControls(simAPI, renderer, config, onReset)` once after WASM is ready.
 *  simAPI   — the simulation.js module exports
 *  renderer — Renderer instance
 *  config   — the same config object used to init the simulation
 *  onReset  — callback that re-initialises the simulation from scratch
 */
export function initControls(simAPI, renderer, config, onReset) {
  // ── Steps per frame ───────────────────────────────────────────────────────
  const stepsSlider = document.getElementById('stepsPerFrame');
  const stepsVal    = document.getElementById('stepsPerFrameVal');
  stepsSlider.addEventListener('input', () => {
    stepsVal.textContent = stepsSlider.value;
  });

  // ── Temperature ───────────────────────────────────────────────────────────
  const tempSlider = document.getElementById('temperature');
  const tempVal    = document.getElementById('temperatureVal');
  tempSlider.value = config.temperature;
  tempVal.textContent = Number(config.temperature).toFixed(2);
  tempSlider.addEventListener('input', () => {
    const t = parseFloat(tempSlider.value);
    tempVal.textContent = t.toFixed(2);
    simAPI.setTemperature(t);
  });

  // ── Chemical potentials (one slider per type) ─────────────────────────────
  const muSection = document.getElementById('muSection');
  config.particle_types.forEach((pt, i) => {
    const label = document.createElement('label');
    label.style.marginTop = '6px';

    const colorDot = `<span style="display:inline-block;width:10px;height:10px;border-radius:50%;background:${pt.color};margin-right:4px;"></span>`;
    const valSpan  = document.createElement('span');
    valSpan.textContent = pt.mu.toFixed(2);

    const slider = document.createElement('input');
    slider.type  = 'range';
    slider.min   = -6;
    slider.max   = 6;
    slider.step  = 0.1;
    slider.value = pt.mu;

    slider.addEventListener('input', () => {
      const mu = parseFloat(slider.value);
      valSpan.textContent = mu.toFixed(2);
      simAPI.setChemicalPotential(i, mu);
    });

    label.innerHTML = `${colorDot}Type ${i} μ`;
    label.appendChild(valSpan);
    label.appendChild(slider);
    muSection.appendChild(label);
  });

  // ── Show bonds ────────────────────────────────────────────────────────────
  const bondsCheck = document.getElementById('showBonds');
  bondsCheck.addEventListener('change', () => {
    renderer.showBonds = bondsCheck.checked;
  });

  // ── Pause / Resume ────────────────────────────────────────────────────────
  const pauseBtn = document.getElementById('pauseBtn');
  pauseBtn.addEventListener('click', () => {
    renderer.paused = !renderer.paused;
    pauseBtn.textContent = renderer.paused ? 'Resume' : 'Pause';
  });

  // ── Reset ─────────────────────────────────────────────────────────────────
  document.getElementById('resetBtn').addEventListener('click', onReset);

  // ── Expose steps-per-frame getter ─────────────────────────────────────────
  return {
    getStepsPerFrame: () => parseInt(stepsSlider.value, 10),
  };
}
