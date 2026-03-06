# Crystal Growth Simulator — Architecture

## Objective

A 2D kinetic Monte Carlo (KMC) crystal growth simulator. Starting from a seed
particle at the origin, the simulation stochastically attaches and detaches
particles according to thermodynamic rates, growing a crystal structure over
time. It supports multiple particle types with configurable radii, interaction
energies, and chemical potentials.

The simulator runs natively via an **eframe/egui** desktop GUI and can also be
compiled to **WebAssembly** for browser-based use.

---

## Physical Model

| Concept | Formula |
|---|---|
| **Bonding shell** | Two particles bond when their center distance is within `contact ± δ`, where `contact = r_i + r_j` |
| **Hard-core overlap** | Overlap when `d < contact − δ` |
| **Binding energy** | `E_bind(i) = Σ_j ε(type_i, type_j)` over bonded neighbours `j` |
| **Attachment rate** | `ν · exp(μ_t / kT)` per candidate site of type `t` |
| **Detachment rate** | `ν · exp(−E_bind / kT)` per particle |
| **Time advance** | BKL algorithm: `Δt = −ln(u) / R_total` |

---

## Module Map

```
src/
├── config.rs       Configuration (SimConfig), deserialized from JSON
├── particle.rs     Particle struct: position, type, radius, bond/overlap tests
├── candidates.rs   CandidateSite struct, circle-intersection geometry, overlap checks
├── spatial.rs      SpatialHash — HashMap-based grid for O(1) neighbor lookups
├── rates.rs        RateCatalog — Fenwick tree for O(log N) rate selection & updates
├── rng.rs          Rng — splitmix64 PRNG (no OS entropy, WASM-safe)
├── kmc.rs          Simulation — the KMC engine tying everything together
├── lib.rs          CrystalSim — wasm-bindgen wrapper for browser use
└── main.rs         Desktop app — eframe/egui GUI with pan, zoom, live controls
```

---

## Key Abstractions

### `SimConfig` (config.rs)
JSON-deserializable configuration. Holds particle type definitions (radius,
color, chemical potential μ), the interaction matrix ε, bonding shell half-width
δ, temperature, attempt frequency ν, and relaxation parameters. Caches derived
values (`max_radius`, `max_cutoff`) via `init_cache()`.

### `Particle` (particle.rs)
A positioned circle with a type ID. Provides `bonds_to()` and `overlaps()`
predicates using the ±δ bonding shell semantics.

### `CandidateSite` (candidates.rs)
A potential attachment position and target type. Generated geometrically: for
each pair of existing particles, compute circle-circle intersections at the
contact radius to find positions where a new particle would bond to both.
Isolated particles get uniformly-spaced arc sites instead.

### `SpatialHash` (spatial.rs)
Unbounded 2D grid hash mapping `(i64, i64)` cell keys to vectors of particle
indices. Cell size equals `max_cutoff` so a one-ring query always covers the
full interaction range. Supports `insert`, `remove`, and `query(cx, cy, r)`.

### `RateCatalog` (rates.rs)
A Fenwick (binary indexed) tree over a flat rate array. Provides O(log N)
`total()`, `select(u)` (binary-lifting search for the event whose prefix sum
reaches `u · R_total`), `set_rate()`, `add_rate()`, and `swap_remove_rate()`.
The simulation maintains two separate catalogs — one for attachment rates
(indexed by candidate site) and one for detachment rates (indexed by particle
slot).

### `Simulation` (kmc.rs)
The core engine. Each KMC step:
1. Computes `R_total = attach_total + detach_total`
2. Draws two uniform randoms: one for event selection, one for time advance
3. Selects attach or detach, then the specific event via Fenwick binary lifting
4. Executes the event (mutating particles, spatial hashes, and rate catalogs)
5. Regenerates candidate sites in the affected neighborhood

Uses swap-remove everywhere to keep index arrays dense and avoid O(N) shifts.

### `CrystalSim` (lib.rs)
Thin `wasm_bindgen` wrapper exposing `step()`, `particle_buffer()`,
`particle_count()`, temperature/μ setters, and type metadata to JavaScript.

### Desktop GUI (main.rs)
eframe/egui app with a side panel for controls (run/pause, steps-per-frame,
temperature, per-type μ sliders, bond visibility toggle) and a central canvas
that renders particles as filled circles with optional bond lines. Supports
mouse drag panning and scroll-wheel zoom.

---

## Remaining Issues & Suggestions

### Bugs

6. **Bond rendering calls `get_neighbors()` per particle every frame** — this
   does a spatial query + distance filter for each particle, every frame. For
   the desktop GUI this is O(N·k) work purely for rendering. Consider caching a
   neighbor list that is rebuilt only when the simulation advances, or skipping
   bond rendering above a particle-count threshold.
