# Crystal Growth Simulator — Design Specification

## Physical Model

**Off-lattice Kinetic Monte Carlo (KMC) in 2D.** Particles occupy continuous positions in ℝ². The simulation advances by stochastic discrete events (attachment/detachment) rather than by integrating equations of motion.

## Particles

- Spherical (circular in 2D), each with a **type** and a **radius** σ_i.
- Position (x, y) ∈ ℝ², no orientation (isotropic). Orientation will be added in a future extension.
- N particle types supported, where N is user-defined.

## Interactions

- **Pairwise, isotropic, constant within a bonding shell.**
- For particle types i and j:
  - Hard-core repulsion for center-to-center distance r < σ_i + σ_j.
  - Bond energy −ε_ij if σ_i + σ_j ≤ r < σ_i + σ_j + δ (bonding shell of width δ).
  - No interaction for r ≥ σ_i + σ_j + δ.
- ε_ij is specified via a symmetric N×N interaction matrix. ε_ij > 0 means attractive.

## KMC Dynamics

### Events

1. **Attachment:** A particle of type t appears at a candidate site on the crystal surface. Rate: r_attach = ν · exp(μ_t / kT), where μ_t is the chemical potential (controls supersaturation) for type t and ν is an attempt frequency.
2. **Detachment:** A particle currently on the crystal surface detaches. Rate: r_detach = ν · exp(−E_bind / kT), where E_bind = Σ_j ε(type_i, type_j) over all bonded neighbors j.

### Event Selection

BKL (Bortz-Kazhdan-Lebowitz) algorithm:
1. Compute cumulative rate R_total = Σ_k r_k over all possible events k.
2. Draw uniform u ∈ (0, 1), select event k such that Σ_{i<k} r_i < u · R_total ≤ Σ_{i≤k} r_i.
3. Execute event k.
4. Advance time by Δt = −ln(u') / R_total, where u' is a second uniform draw.
5. Update affected rates (local update, not full recomputation).

### Candidate Site Generation

For each pair of occupied surface particles (A, B) with |x_A − x_B| < σ_A + σ_C + σ_B + σ_C (for candidate type C):
- Compute the two intersection points of circles of radius σ_A + σ_C and σ_B + σ_C centered at A and B.
- Filter out points that overlap with any existing particle (hard-core violation).
- Each valid intersection point is a candidate attachment site.

Additionally, for isolated surface particles with exposed surface, generate candidate sites at contact distance at uniformly spaced angles (to allow growth initiation from single particles or sparse surfaces).

## Control Parameters

| Parameter | Symbol | Role |
|---|---|---|
| Temperature | T | Controls detachment rate via Boltzmann factor |
| Chemical potential (per type) | μ_t | Controls supersaturation / attachment favorability |
| Bonding shell width | δ | Tolerance for bond formation beyond hard-core contact |
| Interaction matrix | ε_ij | Pairwise bond energies |
| Particle radii | σ_i | Per-type particle size |

## Architecture

### Language & Build

- **Simulation core:** Rust, compiled to WebAssembly via wasm-pack/wasm-bindgen.
- **Math library:** glam (lightweight 2D vector math via `DVec2`). Requires `libm` feature for WASM (software `sin`/`cos`/`sqrt`).
- **Rendering & UI:** JavaScript/HTML, running in the browser.
- **Build pipeline:** wasm-pack produces a .wasm binary + JS glue module. The JS side imports the WASM module and calls into it.

### Module Boundary

The Rust WASM module exposes:
- `init(config: &[u8])` — initialize simulation from a serialized configuration (particle types, interaction matrix, control parameters).
- `step(n: u32)` — advance the simulation by n KMC events.
- `particle_buffer() -> *const f32` — pointer into WASM linear memory containing particle data as a flat array: [x0, y0, type0, radius0, x1, y1, type1, radius1, ...]. The JS side reads this directly from the SharedArrayBuffer without copying.
- `particle_count() -> u32` — current number of particles.
- `simulation_time() -> f64` — current simulation time.
- `set_temperature(t: f64)`, `set_chemical_potential(type_id: u32, mu: f64)` — runtime parameter updates.

### Rendering

- **JS reads the particle buffer from WASM linear memory once per frame.**
- **Renderer:** Canvas 2D initially (simple, sufficient for colored circles + bond lines). Can be swapped to PixiJS for WebGL-accelerated rendering if Canvas 2D becomes a bottleneck at high N.
- Render loop runs on `requestAnimationFrame`, decoupled from simulation steps.
- Bond lines drawn between particles within bonding distance (optional, togglable).

### Frame Loop

```
requestAnimationFrame:
  1. Call wasm.step(events_per_frame)     // Rust advances KMC
  2. Read particle_buffer from WASM memory // Zero-copy via typed array view
  3. Clear canvas
  4. For each particle: draw circle at (x, y) with color(type) and radius
  5. (Optional) Draw bond lines
```

### Local Development

Build and serve locally without any bundler:

```bash
wasm-pack build --target web --release
python3 -m http.server 8080     # serves from repo root; browser loads .wasm via HTTP
# open http://localhost:8080/docs/index.html
```

A bare `file://` URL will not work — browsers refuse to load .wasm without an HTTP server.

### Data Flow

```
  ┌─────────────────────────┐       shared WASM memory        ┌──────────────────┐
  │     Rust (WASM)         │  ──── particle buffer ────────▶  │   JS / Canvas    │
  │                         │                                  │                  │
  │  - KMC engine           │  ◀──── control parameter ─────── │  - Render loop   │
  │  - Spatial hash         │        updates via exports       │  - UI controls   │
  │  - Candidate sites      │                                  │  - Event handlers│
  │  - Rate catalog         │                                  │                  │
  └─────────────────────────┘                                  └──────────────────┘
```

## Data Structures

- **Particle list:** Array of {position, type, neighbor list}.
- **Spatial index:** Grid-based spatial hash for O(1) neighbor queries within cutoff distance.
- **Surface set:** Set of particles with fewer than maximum possible neighbors (exposed to solution).
- **Candidate site cache:** List of valid attachment sites with precomputed ΔE values. Updated locally when a particle is added or removed.
- **Rate catalog:** Array of rates for all current events, maintained incrementally.

## Performance Target

- Thousands of particles (N ≈ 5000+).
- Interactive in a browser (~60fps rendering, with KMC steps decoupled from render loop).
- Hundreds to thousands of KMC events per frame at N = 5000.

## Future Extensions

- **Oriented bonds:** Extend particle state to include orientation θ_i ∈ [0, 2π). Energy becomes ε(r_ij, θ_i, θ_j, r̂_ij, type_i, type_j). Candidate site generation must additionally search/enumerate over orientation of new particle.
- **3D:** Circle-circle intersections become sphere-sphere-sphere intersections; visualization moves to WebGL.