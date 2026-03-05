/// Flat rate catalog for BKL event selection.
/// Uses a Fenwick tree (binary indexed tree) to maintain prefix sums for O(log n) updates and selection.


pub struct RateCatalog {
    rates: Vec<f64>,
    tree: Vec<f64>,
}

impl RateCatalog {
    pub fn new() -> Self {
        Self {
            rates: Vec::new(),
            tree: vec![0.0],
        }
    }
    
    fn update_tree(&mut self, mut idx: usize, delta: f64) {
        idx += 1; // 0-indexed to 1-indexed
        while idx < self.tree.len() {
            self.tree[idx] += delta;
            idx += idx & (!idx + 1); // Equivalent to idx += idx & -idx
        }
    }

    pub fn set_rate(&mut self, candidate_id: usize, rate: f64) {
        let new_rate = rate.max(0.0);
        let delta = new_rate - self.rates[candidate_id];
        self.rates[candidate_id] = new_rate;
        self.update_tree(candidate_id, delta);
    }

    pub fn get_rate(&self, candidate_id: usize) -> f64 {
        self.rates[candidate_id]
    }

    pub fn add_rate(&mut self, rate: f64) -> usize {
        let idx = self.rates.len();
        self.rates.push(rate);
        self.tree.push(0.0);
        self.update_tree(idx, rate);
        idx
    }

    pub fn swap_remove_rate(&mut self, slot: usize) -> Option<usize> {
        let n = self.rates.len();
        assert!(slot < n, "slot out of range");
        let last_idx = n - 1;
        self.update_tree(slot, -self.rates[slot]);

        let moved_last = slot < last_idx;
        if moved_last {
            // Remove last element's contribution from its old position
            let last_rate = self.rates[last_idx];
            self.update_tree(last_idx, -last_rate);

            // Move last element to the hole
            self.rates[slot] = last_rate;
            self.update_tree(slot, last_rate);
        }

        self.rates.pop();
        self.tree.pop();

        if moved_last { Some(last_idx) } else { None }
    }

    // ── Global queries ──────────────────────────────────────────────────────

    pub fn total(&self) -> f64 {
        if self.tree.len() <= 1 { return 0.0; }
        // The total sum is stored in specific nodes or can be calculated via prefix_sum
        // but for a Fenwick tree, the total is the prefix sum of the last index.
        self.prefix_sum(self.rates.len())
    }

    fn prefix_sum(&self, mut idx: usize) -> f64 {
        let mut sum = 0.0;
        while idx > 0 {
            sum += self.tree[idx];
            idx -= idx & (!idx + 1);
        }
        return sum;
    }

    /// Binary lifting on Fenwick tree to find first index where prefix sum >= target.
    pub fn select(&self, u: f64) -> usize {
        let target = u * self.total();
        let mut pos = 0;
        let mut current_sum = 0.0;
        
        // Start from the largest power of 2 less than or equal to the number of elements
        let mut bit_mask = 1 << (usize::BITS - self.rates.len().leading_zeros().saturating_sub(1) - 1);
        if self.rates.is_empty() { return 0; }

        while bit_mask > 0 {
            let next_pos = pos + bit_mask;
            if next_pos < self.tree.len() && current_sum + self.tree[next_pos] < target {
                pos = next_pos;
                current_sum += self.tree[pos];
            }
            bit_mask >>= 1;
        }

        // pos is the last index where prefix_sum < target, so pos is the 0-based index we need
        pos.min(self.rates.len().saturating_sub(1))
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    // Linear ground-truth implementation (pure, no particle/type awareness)
    struct LinearRateCatalog {
        rates: Vec<f64>,
    }

    impl LinearRateCatalog {
        fn new() -> Self {
            Self { rates: Vec::new() }
        }
        fn add_rate(&mut self, rate: f64) {
            self.rates.push(rate.max(0.0));
        }
        fn set_rate(&mut self, idx: usize, rate: f64) {
            self.rates[idx] = rate.max(0.0);
        }
        fn remove_rate(&mut self, slot: usize) {
            self.rates.swap_remove(slot);
        }
        fn len(&self) -> usize {
            self.rates.len()
        }
        fn total(&self) -> f64 {
            self.rates.iter().sum()
        }
        fn select(&self, u: f64) -> usize {
            let target = u * self.total();
            let mut cumsum = 0.0;
            for (i, &r) in self.rates.iter().enumerate() {
                cumsum += r;
                if cumsum >= target { return i; }
            }
            self.rates.len().saturating_sub(1)
        }
    }

    #[test]
    fn test_fenwick_vs_linear() {
        let mut rng = rand::thread_rng();
        let mut fenwick = RateCatalog::new();
        let mut linear = LinearRateCatalog::new();

        for _ in 0..10_000 {
            let action = rng.gen_range(0..4);
            match action {
                0 if linear.len() > 0 => { // Update a rate
                    let idx = rng.gen_range(0..linear.len());
                    let rate = rng.gen::<f64>() * 10.0;
                    fenwick.set_rate(idx, rate);
                    linear.set_rate(idx, rate);
                }
                1 => { // Add a rate
                    let rate = rng.gen::<f64>() * 10.0;
                    fenwick.add_rate(rate);
                    linear.add_rate(rate);
                }
                2 if linear.len() > 0 => { // Remove a rate
                    let slot = rng.gen_range(0..linear.len());
                    fenwick.swap_remove_rate(slot);
                    linear.remove_rate(slot);
                }
                _ => {}
            }

            // Compare totals (with small epsilon for floating point jitter)
            let diff = (fenwick.total() - linear.total()).abs();
            assert!(diff < 1e-9, "Total mismatch: F: {}, L: {}", fenwick.total(), linear.total());

            // Compare selection
            if fenwick.total() > 0.0 {
                let u = rng.gen::<f64>();
                let f_idx = fenwick.select(u);
                let l_idx = linear.select(u);
                assert_eq!(f_idx, l_idx, "Selection mismatch at u={}: F: {}, L: {}", u, f_idx, l_idx);
            }
        }
    }
}



// pub struct RateCatalog {
//     rates: Vec<f64>,
//     n_types: usize,
// }

// impl RateCatalog {
//     pub fn new(n_types: usize) -> Self {
//         Self { rates: vec![0.0; n_types], n_types }
//     }

//     pub fn set_attach_rate(&mut self, type_id: usize, rate: f64) {
//         self.rates[type_id] = rate.max(0.0);
//     }

//     pub fn attach_rate(&self, type_id: usize) -> f64 {
//         self.rates[type_id]
//     }

//     // ── Detachment ──────────────────────────────────────────────────────────

//     /// Append a new particle slot and return its slot index (0-based among particles).
//     pub fn add_particle(&mut self, detach_rate: f64) -> usize {
//         self.rates.push(detach_rate.max(0.0));
//         self.n_particles() - 1
//     }

//     pub fn set_detach_rate(&mut self, particle_slot: usize, rate: f64) {
//         self.rates[self.n_types + particle_slot] = rate.max(0.0);
//     }

//     pub fn detach_rate(&self, particle_slot: usize) -> f64 {
//         self.rates[self.n_types + particle_slot]
//     }

//     /// Swap-remove particle slot `i`. Returns the old slot of the particle
//     /// moved into `i` (previously the last slot), or `None` if `i` was last.
//     pub fn remove_particle(&mut self, particle_slot: usize) -> Option<usize> {
//         let n = self.n_particles();
//         assert!(particle_slot < n, "particle slot out of range");
//         let idx = self.n_types + particle_slot;
//         let last_idx = self.n_types + n - 1;
//         if particle_slot < n - 1 {
//             self.rates.swap(idx, last_idx);
//         }
//         self.rates.pop();
//         if particle_slot < n - 1 { Some(n - 1) } else { None }
//     }

//     // ── Global queries ──────────────────────────────────────────────────────

//     pub fn total(&self) -> f64 {
//         self.rates.iter().sum()
//     }

//     /// BKL selection: given `u` uniform in [0,1), return the event index such
//     /// that it is the first index where the prefix sum reaches u * R_total.
//     /// Index < n_types → attachment event for that type.
//     /// Index >= n_types → detachment for particle (index - n_types).
//     pub fn select(&self, u: f64) -> usize {
//         let target = u * self.total();
//         let mut cumsum = 0.0;
//         for (i, &r) in self.rates.iter().enumerate() {
//             cumsum += r;
//             if cumsum >= target {
//                 return i;
//             }
//         }
//         self.rates.len().saturating_sub(1)
//     }

//     pub fn n_types(&self) -> usize {
//         self.n_types
//     }

//     pub fn n_particles(&self) -> usize {
//         self.rates.len().saturating_sub(self.n_types)
//     }
// }
