# Current Physics Model (As Implemented)

This document describes the **current** simulation physics in `crystal_sim` as implemented in code, not an idealized target model.

## 1. High-level model

The simulator combines two layers:

1. **KMC event dynamics** (stochastic):
- Attach a particle at a candidate site, or
- Detach an existing particle.

2. **Local mechanical relaxation** (deterministic, optional):
- After each attachment, run damped translational + rotational relaxation for nearby particles.

Core flow lives in `src/kmc.rs` (`Simulation::step`, `attach`, `detach`, `relax_new_particle_fixed_grid`).

## 2. Particle state

Each particle has:
- Position: `(x, y)`
- Type id
- Radius
- Orientation angle (radians)

See `src/particle.rs`.

## 3. Energetics used for KMC rates

### 3.1 Binding energy

For a particle `p`, the binding energy `E_bind` is computed from neighbors in a cutoff shell (`binding_energy` in `src/kmc.rs`):

- If patch mode is active (`config.has_patches()`):
  - Pair energy is `forces::patchy_pair_energy(...)`.
- Otherwise:
  - Pair contributes `epsilon[type_i][type_j]` only if distance is in the bond shell.

### 3.2 KMC rates

Rates are Arrhenius-like (in `src/config.rs`):

- Attach rate at candidate site:
  - `k_attach = nu * exp((mu_type + E_bind / 2) / T)`
- Detach rate for a particle:
  - `k_detach = nu * exp((-E_bind / 2) / T)`

Event selection uses a Fenwick-tree catalog for attachment candidates (`src/rates.rs`) and weighted sampling over particle detach rates in the grid (`src/spatial.rs`).

## 4. Pair interaction physics

## 4.1 Non-patch fallback

If patch interactions are not configured for a pair:

- Energy for KMC bonding:
  - Uses scalar `epsilon` matrix with bond-shell distance gating.
- Relaxation force:
  - Uses isotropic Lennard-Jones gradient/force fallback through `patchy_force_torque` fallback path.
  - No torque contribution.

## 4.2 Patchy interactions

Patch data is configured in `src/config.rs`:
- Each patch has type + body-frame angle.
- Patch-patch interactions define:
  - `epsilon`
  - Angular widths (`sigma_a`, `sigma_b`)
  - Pair cutoff multiplier.

For a compatible patch pair `(a on i, b on j)`:

- Angular weights:
  - `g_ia = exp(-theta_ia^2 / (2 * sigma_a^2))`
  - `g_jb = exp(-theta_jb^2 / (2 * sigma_b^2))`
  - `g = g_ia * g_jb`

- KMC pair energy contribution:
  - `U_pair = epsilon * bump(r; r_contact, cutoff) * g`
  - `bump(...)` is a compact, strictly positive radial function centered at `r_contact`.

- Relaxation pair mechanics (`forces::patchy_force_torque`):
  - Backbone isotropic repulsion term (`r^-12`-style) to prevent overlap.
  - Patch radial term from `d/dr [epsilon * bump(r) * g]`.
  - Angular correction to linear force from orientation alignment gradients.
  - Torques on both particles:
    - `tau_i ~ epsilon * bump(r) * g * theta_ia / sigma_a^2`
    - `tau_j ~ epsilon * bump(r) * g * theta_jb / sigma_b^2`

Implementation references:
- `src/forces.rs` (`patchy_pair_energy`, `patchy_force_torque`, LJ helpers)

## 4.3 Sign convention and epsilon case-by-case behavior

Sign convention in code:
- `lj_force_vec`, `lj_force_vec_full`, and `repulsive_force_vec` return `grad_i(U)`, not physical force.
- Physical force is `F_i = -grad_i(U)`.
- `patchy_force_torque` already returns physical force on `i` plus torques (`tau_i`, `tau_j`).

### 4.3.1 Meaning of epsilon in this codebase

- Scalar species matrix `config.epsilon[i][j]`:
  - Used only in fallback paths (no patch interaction for that pair path).
  - In KMC fallback (`patchy_pair_energy` when one side has no patches), contributes a constant inside bond shell: `+epsilon[i][j]`.
  - In force fallback (`patchy_force_torque` when one side has no patches), used in full LJ force.

- Patch interaction `patch_interactions[k].epsilon`:
  - In KMC energy (`patchy_pair_energy`), scales a strictly positive radial bump times angular weight.
  - In relaxation (`patchy_force_torque`), scales bump-based force/torque terms.
  - Positive epsilon increases KMC binding contribution; negative epsilon decreases it.

Important: current implementation still treats positive epsilon as binding-favorable in KMC energy. If older docs/configs assume "negative = attractive", those docs are out of sync with current code.

### 4.3.2 Case matrix: what each specification does

Case A: patch mode effectively off (`!has_patches()`)
- Trigger: `patch_interactions` empty, or no particle type has patches.
- KMC pair energy: shell-gated constant `+epsilon[i][j]` (when `|r-r_contact| <= delta`).
- Relaxation force: isotropic LJ fallback with `epsilon[i][j]`; no torque.
- Effect of epsilon:
  - `epsilon > 0`: favors binding (higher `E_bind`, higher attach rate, lower detach rate).
  - `epsilon = 0`: neutral for that pair.
  - `epsilon < 0`: disfavors binding (lower `E_bind`, lower attach rate, higher detach rate).

Case B: patch mode on, but one particle has zero patches
- Trigger: one side `patches.is_empty()`.
- KMC pair energy: same scalar shell fallback as Case A.
- Relaxation: same isotropic LJ fallback as Case A; no torque.
- Practical implication: "mixed" patch/non-patch interactions ignore patch table and act like isotropic species epsilon.

Case C: both particles have patches, but no matching patch-type pair in the patch LUT
- Trigger: the directed LUT entry for that patch-type pair is disabled (no configured interaction).
- KMC pair energy: zero from patchy term.
- Relaxation: only isotropic backbone repulsion (`r^-12` style) is applied.
- Practical implication: particles can still avoid overlap, but get no patch attraction and no orientational torque.

Case D: matching patch pair with `patch_epsilon > 0`
- KMC pair energy: `+epsilon * bump(r) * g(theta_i, theta_j)` inside interaction cutoff.
- Relaxation:
  - radial force from bump derivative scaled by `g`,
  - angular correction to linear force,
  - nonzero torques that align patches toward each other.
- Strongest binding occurs for:
  - distances near `r_contact` (bump peak),
  - small angular mismatch (large Gaussian weight `g`),
  - not exceeding per-interaction cutoff (`r < r_contact * cutoff`).

Case E: matching patch pair with `patch_epsilon = 0`
- KMC pair energy: zero from that interaction.
- Relaxation: zero patch bump/torque from that interaction; only backbone repulsion remains.
- Practical implication: interaction exists structurally in table, but is physically inert.

Case F: matching patch pair with `patch_epsilon < 0`
- KMC pair energy: negative contribution (`epsilon * positive_bump * g`).
- Relaxation: radial and torque patch terms flip sign relative to Case D.
- Practical implication: this behaves like "anti-binding" for KMC energetic weighting; use intentionally.

### 4.3.3 Other patch specification fields that change behavior

- `angular_width_deg` (converted to radians in cache):
  - Smaller width -> sharper directionality, rapid dropoff with angular error.
  - Larger width -> broader acceptance, weaker orientation selectivity.
  - Widths apply per side of the interaction in the order declared by `types`.

- `types` ordering:
  - Lookup is symmetric (`A,B` and `B,A` are both mapped), but sigma assignment is side-aware.
  - `sigma_for_pair` ensures sigma associated with the i-side patch type is applied to `theta_ia`.

- `cutoff` (multiplied by `r_contact`):
  - Hard distance gate for this patch interaction.
  - If `r >= r_contact * cutoff`, that patch pair contributes exactly zero even if perfectly aligned.

## 5. Candidate site generation (attachment geometry)

Candidate regeneration (`regen_candidates_batch` in `src/kmc.rs`) builds candidate sites near changed regions:

1. Remove stale candidates in affected regions.
2. Build local focal particle set.
3. Build overlap check arrays.
4. Generate candidate positions by:
- Circle-circle intersections around focal pairs.
- Arc sampling around focal particles with blocked-angle filtering.

Overlap is rejected by hard-core check (`site_has_overlap` in `src/candidates.rs`).

Important current behavior:
- Candidate orientation is computed per candidate (`best_candidate_orientation`) by sampling plausible orientations and selecting the one with best local patchy energy.
- That orientation is stored on `CandidateSite` and used in candidate energy/rate calculation.
- If patch mode is inactive (or candidate type has no patches), orientation defaults to `0.0`.

## 6. Mechanical relaxation step

Relaxation runs after attach if `relax_steps > 0` (`relax_new_particle_fixed_grid` in `src/kmc.rs`).

Current mechanics:

1. Define fixed active/initialized cell windows around the new particle.
2. Mark particles in active-zone cells as active.
3. For each relax iteration:
- Zero per-particle forces and torques.
- Accumulate pair interactions for active particles against nearby cells.
- Use symmetry skip for active-active pairs (compute each pair once).
- Apply equal/opposite force updates and torque updates.
- Integrate translation and rotation with damping/friction thresholds.
- Commit updated positions/orientations back to grid.
- Stop early if converged.

Translational update (per active particle):
- Clamp force magnitude.
- Apply static-friction-style thresholding.
- `v_new = damping * v_old + alpha * f_eff` (with sign consistency guard).
- `x_new = x_old + v_new`.

Rotational update (per active particle):
- Clamp torque.
- Threshold small torques.
- `omega_new = damping * omega_old + alpha * tau_eff` (with sign consistency guard).
- `ori_new = ori_old + omega_new`.

## 7. Detach/attach side effects after topology change

After attach/detach (+ optional relax), the sim updates:

- Detach rates for affected particles (`update_detach_rate`).
- Local attachment candidate catalog (`regen_candidates_batch`).

This keeps KMC event weights consistent with the latest local structure.

## 8. Spatial representation

`ParticleGrid` (`src/spatial.rs`) is a cell-hashed, append-only cell store with per-cell arrays for:
- particles
- detach rates
- linear velocities
- target positions
- angular velocities
- target orientations
- active/frozen flags
- force and torque buffers.

It supports neighborhood queries, weighted detach sampling, and cross-cell commit of moved particles.

## 9. Notes relevant for refactor

These are "as implemented" details you will likely want to revisit during refactoring:

- Physics is hybrid: stochastic KMC events + local damped relaxation.
- Patch mechanics are active in relaxation/energy if patch config exists.
- KMC rates depend on local energy via the formulas above, not explicit barriers per patch pair.
- Relaxation window is fixed-grid around new particles (not globally adaptive).

