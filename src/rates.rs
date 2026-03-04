/// Flat rate catalog for BKL event selection.
///
/// Layout (indices):
///   [0 .. n_types)            — aggregate attachment rates for each particle type
///   [n_types .. n_types+N)    — per-particle detachment rates
///
/// Attachment rate for type t  = n_candidate_sites[t] * ν * exp(μ_t / kT)
/// Detachment rate for particle i = ν * exp(-E_bind_i / kT)
///
/// Event selection is an O(N) linear scan — sufficient for N ≤ ~10 000 events.
/// A Fenwick tree can be substituted later without changing the public API.
pub struct RateCatalog {
    rates: Vec<f64>,
    n_types: usize,
}

impl RateCatalog {
    pub fn new(n_types: usize) -> Self {
        Self { rates: vec![0.0; n_types], n_types }
    }

    // ── Attachment ──────────────────────────────────────────────────────────

    pub fn set_attach_rate(&mut self, type_id: usize, rate: f64) {
        self.rates[type_id] = rate.max(0.0);
    }

    pub fn attach_rate(&self, type_id: usize) -> f64 {
        self.rates[type_id]
    }

    // ── Detachment ──────────────────────────────────────────────────────────

    /// Append a new particle slot and return its slot index (0-based among particles).
    pub fn add_particle(&mut self, detach_rate: f64) -> usize {
        self.rates.push(detach_rate.max(0.0));
        self.n_particles() - 1
    }

    pub fn set_detach_rate(&mut self, particle_slot: usize, rate: f64) {
        self.rates[self.n_types + particle_slot] = rate.max(0.0);
    }

    pub fn detach_rate(&self, particle_slot: usize) -> f64 {
        self.rates[self.n_types + particle_slot]
    }

    /// Swap-remove particle slot `i`. Returns the old slot of the particle
    /// moved into `i` (previously the last slot), or `None` if `i` was last.
    pub fn remove_particle(&mut self, particle_slot: usize) -> Option<usize> {
        let n = self.n_particles();
        assert!(particle_slot < n, "particle slot out of range");
        let idx = self.n_types + particle_slot;
        let last_idx = self.n_types + n - 1;
        if particle_slot < n - 1 {
            self.rates.swap(idx, last_idx);
        }
        self.rates.pop();
        if particle_slot < n - 1 { Some(n - 1) } else { None }
    }

    // ── Global queries ──────────────────────────────────────────────────────

    pub fn total(&self) -> f64 {
        self.rates.iter().sum()
    }

    /// BKL selection: given `u` uniform in [0,1), return the event index such
    /// that it is the first index where the prefix sum reaches u * R_total.
    /// Index < n_types → attachment event for that type.
    /// Index >= n_types → detachment for particle (index - n_types).
    pub fn select(&self, u: f64) -> usize {
        let target = u * self.total();
        let mut cumsum = 0.0;
        for (i, &r) in self.rates.iter().enumerate() {
            cumsum += r;
            if cumsum >= target {
                return i;
            }
        }
        self.rates.len().saturating_sub(1)
    }

    pub fn n_types(&self) -> usize {
        self.n_types
    }

    pub fn n_particles(&self) -> usize {
        self.rates.len().saturating_sub(self.n_types)
    }
}
