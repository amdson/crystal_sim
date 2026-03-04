# Crystal Growth Simulator — Project Structure

## Overview

Rust/WASM simulation core + vanilla JS/HTML frontend. No bundler required; built with `wasm-pack` and served via a static HTTP server.

---

## Directory Tree

```
crystal_sim/
│
├── Cargo.toml                   # Rust workspace manifest
├── Cargo.lock
├── crystal_sim_roadmap.md       # Design spec (this project's source of truth)
├── project_structure.md         # This file
│
├── src/                         # Rust simulation core (compiled to WASM)
│   ├── lib.rs                   # WASM entry point; exposes public #[wasm_bindgen] exports
│   │
│   ├── config.rs                # Config deserialization: particle types, ε matrix, control params
│   │
│   ├── particle.rs              # Particle struct: position, type_id, radius, neighbor list
│   │
│   ├── spatial.rs               # Grid-based spatial hash — O(1) neighbor queries within cutoff
│   │
│   ├── surface.rs               # Surface set: tracks particles with exposed bonds
│   │
│   ├── candidates.rs            # Candidate site generation (circle-circle intersections + isolated)
│   │
│   ├── rates.rs                 # Rate catalog: per-event rates, cumulative sum, incremental update
│   │
│   ├── kmc.rs                   # KMC engine: BKL event selection, time advance, event execution
│   │
│   └── rng.rs                   # WASM-compatible PRNG (no OS entropy; seeded deterministically)
│
├── pkg/                         # wasm-pack output (generated, do not edit)
│   ├── crystal_sim_bg.wasm      # Compiled WASM binary
│   ├── crystal_sim_bg.js        # WASM glue (imports)
│   ├── crystal_sim.js           # JS bindings (wasm-bindgen generated)
│   └── crystal_sim.d.ts         # TypeScript declarations (optional, for IDE hints)
│
└── docs/                        # Frontend — served by the HTTP dev server
    ├── index.html               # Main page: canvas + control panel layout
    │
    ├── css/
    │   └── styles.css           # Layout, panel, slider, button styling
    │
    └── js/
        ├── main.js              # Entry point: loads WASM module, wires everything up
        ├── simulation.js        # Wrapper around WASM exports (init, step, buffer access)
        ├── renderer.js          # Canvas 2D render loop (requestAnimationFrame, draw particles/bonds)
        └── controls.js          # UI event handlers: sliders, buttons → set_temperature, set_mu, etc.
```

---

## File Responsibilities

### Rust (`src/`)

| File | Responsibility |
|---|---|
| `lib.rs` | Declares the crate as `cdylib`. Exports `init`, `step`, `particle_buffer`, `particle_count`, `simulation_time`, `set_temperature`, `set_chemical_potential` via `#[wasm_bindgen]`. Owns the global simulation state (`static mut SIM`). |
| `config.rs` | Defines `SimConfig` struct (particle types, radii, N×N ε matrix, δ, T, μ per type, ν). Deserializes from a byte slice passed by JS (likely JSON or a compact binary format like bincode). |
| `particle.rs` | `Particle` struct with `pos: DVec2`, `type_id: u32`, `radius: f64`, and a neighbor index list. Helper: `bond_energy(&self, other: &Particle, eps: f64) -> f64`. |
| `spatial.rs` | `SpatialHash` struct. Cell size = max cutoff distance. Methods: `insert`, `remove`, `neighbors_within(pos, radius)`. |
| `surface.rs` | `SurfaceSet` wrapping a `HashSet<usize>` of particle indices. Tracks which particles are eligible for detachment or neighbor-site generation. Updated on every attach/detach event. |
| `candidates.rs` | `CandidateSite` struct: `pos: DVec2`, `type_id: u32`, `delta_e: f64`. `generate_sites(surface, particles, spatial, config) -> Vec<CandidateSite>`. Circle-circle intersection geometry + isolated-particle arc sampling. |
| `rates.rs` | `RateCatalog` holding a flat `Vec<f64>` of event rates and a Fenwick tree (or prefix-sum array) for O(log N) cumulative-rate queries. Methods: `total()`, `select(u)`, `update(idx, new_rate)`. |
| `kmc.rs` | `Simulation` struct composing all of the above. Methods: `step(n)` runs n BKL iterations; each iteration: sample event, execute (attach or detach), locally update spatial hash, surface set, candidate cache, and rate catalog. |
| `rng.rs` | Lightweight xoshiro256** or similar PRNG. No `getrandom`/OS dependency (WASM-safe). Seeded from a user-supplied u64 or a fixed default. |

### JS (`docs/js/`)

| File | Responsibility |
|---|---|
| `main.js` | `import init from '../../pkg/crystal_sim.js'`. Calls `await init()`, builds a config JSON blob, calls `wasm.init(configBytes)`, then starts the render loop. |
| `simulation.js` | Thin wrapper: exposes `stepSim(n)`, `getParticles()` (returns a `Float32Array` view into WASM memory), `getParticleCount()`, `getTime()`. Isolates WASM API from render/control code. |
| `renderer.js` | `Renderer` class. `draw(particleBuffer, count)`: clears canvas, iterates buffer (stride 4: x, y, type, radius), draws circles with per-type color, optionally draws bond lines. |
| `controls.js` | Reads DOM inputs (temperature slider, μ sliders, steps-per-frame input, pause/resume button). On change: calls `wasm.set_temperature(val)` or `wasm.set_chemical_potential(typeId, val)`. |

---

## Build & Run

```bash
# 1. Build WASM
wasm-pack build --target web --release

# 2. Serve
python3 -m http.server 8080

# 3. Open
# http://localhost:8080/docs/index.html
```

---

## Key Dependencies (`Cargo.toml`)

```toml
[dependencies]
wasm-bindgen = "0.2"
glam = { version = "0.29", features = ["libm"] }   # DVec2; libm required for WASM
serde = { version = "1", features = ["derive"] }
serde_json = "1"

[lib]
crate-type = ["cdylib"]
```

---

## Notes

- `pkg/` is gitignored (generated artifact).
- The particle buffer layout is a flat `f32` array with stride 4: `[x, y, type_id_as_f32, radius, ...]`. JS reads this via `new Float32Array(wasm.memory.buffer, wasm.particle_buffer(), count * 4)` — zero-copy.
- Candidate sites are recomputed locally (not globally) on each attach/detach to maintain performance at N ≈ 5000.
- The Fenwick tree in `rates.rs` keeps event selection at O(log N) instead of O(N) linear scan.
