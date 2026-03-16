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
  let currentConfig = config;
  // ── Steps per frame ───────────────────────────────────────────────────────
  const stepsSlider = document.getElementById('stepsPerFrame');
  const stepsVal    = document.getElementById('stepsPerFrameVal');
  stepsSlider.addEventListener('input', () => {
    stepsVal.textContent = stepsSlider.value;
  });

  function setStepsFromConfig(cfg) {
    const spf = cfg.steps_per_frame ?? 200;
    stepsSlider.value = spf;
    stepsVal.textContent = spf;
  }
  setStepsFromConfig(config);

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

  // ── Show candidates ───────────────────────────────────────────────────────
  const candidatesCheck = document.getElementById('showCandidates');
  candidatesCheck.addEventListener('change', () => {
    renderer.showCandidates = candidatesCheck.checked;
  });

  // ── Pause / Run ───────────────────────────────────────────────────────────
  const pauseBtn = document.getElementById('pauseBtn');
  renderer.paused = true;        // start paused so Step button is the default
  pauseBtn.textContent = 'Run';
  pauseBtn.addEventListener('click', () => {
    renderer.paused = !renderer.paused;
    pauseBtn.textContent = renderer.paused ? 'Run' : 'Pause';
  });

  // ── Step (advance exactly one frame's worth of steps while paused or not) ─
  document.getElementById('stepBtn').addEventListener('click', () => {
    simAPI.stepSim(parseInt(stepsSlider.value, 10));
  });

  // ── Reset ─────────────────────────────────────────────────────────────────
  document.getElementById('resetBtn').addEventListener('click', () => {
    renderer.paused = true;
    pauseBtn.textContent = 'Run';
    onReset(currentConfig);
  });

  // ── Config loader ─────────────────────────────────────────────────────────
  const configSelect  = document.getElementById('configSelect');
  const loadConfigBtn = document.getElementById('loadConfigBtn');
  const configJson    = document.getElementById('configJson');
  const applyJsonBtn  = document.getElementById('applyJsonBtn');
  const jsonError     = document.getElementById('jsonError');

  function applyConfig(cfg) {
    currentConfig = cfg;
    configJson.value = JSON.stringify(cfg, null, 2);
    jsonError.textContent = '';
    setStepsFromConfig(cfg);
    renderer.paused = true;
    pauseBtn.textContent = 'Run';
    // Rebuild mu sliders for new particle types
    const muSection = document.getElementById('muSection');
    muSection.innerHTML = '<h3>Chemical potentials (μ)</h3>';
    cfg.particle_types.forEach((pt, i) => {
      const label = document.createElement('label');
      label.style.marginTop = '6px';
      const colorDot = `<span style="display:inline-block;width:10px;height:10px;border-radius:50%;background:${pt.color};margin-right:4px;"></span>`;
      const valSpan  = document.createElement('span');
      valSpan.textContent = (pt.mu ?? 0).toFixed(2);
      const slider = document.createElement('input');
      slider.type = 'range'; slider.min = -6; slider.max = 6; slider.step = 0.1;
      slider.value = pt.mu ?? 0;
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
    onReset(cfg);
  }

  // Populate textarea with current config on page load
  configJson.value = JSON.stringify(currentConfig, null, 2);

  loadConfigBtn.addEventListener('click', async () => {
    try {
      const cfg = await fetch(configSelect.value).then(r => {
        if (!r.ok) throw new Error(`HTTP ${r.status}`);
        return r.json();
      });
      applyConfig(cfg);
    } catch (e) {
      jsonError.textContent = `Failed to load: ${e.message}`;
    }
  });

  applyJsonBtn.addEventListener('click', () => {
    try {
      const cfg = JSON.parse(configJson.value);
      applyConfig(cfg);
    } catch (e) {
      jsonError.textContent = `JSON parse error: ${e.message}`;
    }
  });

  // ── Expose steps-per-frame getter ─────────────────────────────────────────
  return {
    getStepsPerFrame: () => parseInt(stepsSlider.value, 10),
    onNewConfig: (cfg) => { currentConfig = cfg; configJson.value = JSON.stringify(cfg, null, 2); setStepsFromConfig(cfg); },
  };
}
