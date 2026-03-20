mod candidates;
pub mod config;
pub mod editor_sim;
pub mod forces;
pub mod kmc;
pub mod particle;
mod rates;
mod rng;
mod spatial;

pub use editor_sim::EditorSim;

use kmc::Simulation;
use config::SimConfig;
use wasm_bindgen::prelude::*;

/// WASM-exposed crystal growth simulator.
///
/// JS usage:
/// ```js
/// import init, { CrystalSim } from './pkg/crystal_sim.js';
/// await init();
/// const sim = new CrystalSim(JSON.stringify(config));
/// sim.step(500);
/// const ptr = sim.particle_buffer();
/// const buf = new Float32Array(sim.memory().buffer, ptr, sim.particle_count() * 5);
/// ```
#[wasm_bindgen]
pub struct CrystalSim {
    inner: Simulation,
}

#[wasm_bindgen]
impl CrystalSim {
    /// Create a new simulation from a JSON config string.
    #[wasm_bindgen(constructor)]
    pub fn new(config_json: &str) -> Result<CrystalSim, JsValue> {
        let mut config: SimConfig = serde_json::from_str(config_json)
            .map_err(|e| JsValue::from_str(&e.to_string()))?;
        config.init_cache();
        Ok(CrystalSim { inner: Simulation::new(config) })
    }

    /// Advance the simulation by `n` KMC events.
    pub fn step(&mut self, n: u32) {
        self.inner.step(n);
    }

    /// Byte offset into WASM linear memory of the particle buffer.
    /// Buffer layout: [x, y, type_id, radius, orientation]  per particle (f32, stride 5)
    ///
    /// JS: `new Float32Array(wasm.memory.buffer, sim.particle_buffer(), sim.particle_count() * 5)`
    pub fn particle_buffer(&self) -> u32 {
        self.inner.particle_buf_ptr() as u32
    }

    /// Number of particles currently in the simulation.
    pub fn particle_count(&self) -> u32 {
        self.inner.particle_count()
    }

    /// Simulated time in KMC time units.
    pub fn simulation_time(&self) -> f64 {
        self.inner.simulation_time()
    }

    /// Update temperature at runtime.
    pub fn set_temperature(&mut self, t: f64) {
        // self.inner.set_temperature(t);
        //pass
        return;
    }

    // /// Update chemical potential for a particle type at runtime.
    // pub fn set_chemical_potential(&mut self, type_id: u32, mu: f64) {
    //     self.inner.set_chemical_potential(type_id as usize, mu);
    // }

    /// JSON array of per-type metadata: [{"color":"#hex","radius":r}, ...]
    pub fn type_metadata_json(&self) -> String {
        self.inner.type_metadata_json()
    }

}
