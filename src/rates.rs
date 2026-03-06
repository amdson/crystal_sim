/// Flat rate catalog for BKL event selection.
/// Uses a Fenwick tree (binary indexed tree) to maintain prefix sums for
/// O(log n) updates and selection.
///
/// Internal tree is 0-indexed using the `i | (i+1)` / `i & (i-1)` convention.

pub struct RateCatalog {
    rates: Vec<f64>,
    tree: Vec<f64>,
}

impl RateCatalog {
    pub fn new() -> Self {
        Self {
            rates: Vec::new(),
            tree: Vec::new(),
        }
    }

    /// Add `delta` to the tree node at `index` and propagate upward.
    fn add_at(&mut self, mut index: usize, delta: f64) {
        while index < self.tree.len() {
            self.tree[index] += delta;
            index |= index + 1; // move to parent
        }
    }

    /// Prefix sum of elements [0, index) (exclusive upper bound).
    fn prefix_sum(&self, mut index: usize) -> f64 {
        let mut sum = 0.0;
        while index > 0 {
            sum += self.tree[index - 1];
            index &= index - 1; // strip lowest set bit
        }
        sum
    }

    pub fn set_rate(&mut self, slot: usize, rate: f64) {
        let new_rate = rate.max(0.0);
        let delta = new_rate - self.rates[slot];
        self.rates[slot] = new_rate;
        self.add_at(slot, delta);
    }

    pub fn get_rate(&self, slot: usize) -> f64 {
        self.rates[slot]
    }

    /// Append a new rate and return its index.
    pub fn add_rate(&mut self, rate: f64) -> usize {
        let index = self.rates.len();
        self.rates.push(rate);
        self.tree.push(0.0); // placeholder

        // Absorb children: for each set bit below `index`, add the child node
        let lower_one_bits = (!index).trailing_zeros();
        for i in 0..lower_one_bits {
            let child = index & !(1 << i);
            self.tree[index] += self.tree[child];
        }
        // Finally add this element's own value
        self.tree[index] += rate;

        index
    }

    /// Remove the last element from the tree (inverse of `push`).
    fn pop(&mut self) {
        if self.tree.is_empty() {
            return;
        }
        let last_idx = self.tree.len() - 1;
        // Reconstruct the raw value at last_idx from prefix sums
        let val = self.prefix_sum(last_idx + 1) - self.prefix_sum(last_idx);
        // Subtract it from tree nodes so the pop is clean
        self.add_at(last_idx, -val);
        self.tree.pop();
        self.rates.pop();
    }

    /// Swap-remove: delete the rate at `slot` by moving the last element into
    /// its place (mirrors `Vec::swap_remove`).  Returns `Some(old_last_index)`
    /// if a move happened, or `None` if `slot` was already the last element.
    pub fn swap_remove_rate(&mut self, slot: usize) -> Option<usize> {
        let n = self.rates.len();
        assert!(slot < n, "slot out of range");
        let last_idx = n - 1;

        if slot == last_idx {
            self.pop();
            return None;
        }

        // Grab the last element's rate, then pop it off
        let last_rate = self.rates[last_idx];
        self.pop();

        // Overwrite the hole at `slot` with the former last element's rate
        let old_rate = self.rates[slot];
        self.rates[slot] = last_rate;
        self.add_at(slot, last_rate - old_rate);

        Some(last_idx)
    }

    // ── Global queries ──────────────────────────────────────────────────────

    pub fn total(&self) -> f64 {
        self.prefix_sum(self.rates.len())
    }

    /// BKL selection: given `u` uniform in [0,1), return the index such that
    /// prefix_sum(index) <= u * total < prefix_sum(index + 1).
    ///
    /// Uses the `index_of` binary-descent algorithm from the reference.
    pub fn select(&self, u: f64) -> usize {
        if self.rates.is_empty() {
            return 0;
        }
        let mut target = u * self.total();
        let mut index = 0;
        let mut probe = most_significant_bit(self.tree.len()) * 2;

        while probe > 0 {
            let lsb = probe & probe.wrapping_neg();
            let half_lsb = lsb / 2;
            let other_half_lsb = lsb - half_lsb;

            if let Some(&value) = self.tree.get(probe - 1) {
                if value < target {
                    index = probe;
                    target -= value;
                    probe += half_lsb;
                    if half_lsb > 0 {
                        continue;
                    }
                }
            }

            if lsb % 2 > 0 {
                break;
            }
            probe -= other_half_lsb;
        }

        index.min(self.rates.len().saturating_sub(1))
    }
}

const fn most_significant_bit(n: usize) -> usize {
    if n == 0 {
        0
    } else {
        1 << (usize::BITS - 1 - n.leading_zeros())
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

    #[test]
    fn test_add_only_total() {
        let mut f = RateCatalog::new();
        f.add_rate(1.0);
        assert!((f.total() - 1.0).abs() < 1e-12, "total after 1 add: {}", f.total());
        f.add_rate(2.0);
        assert!((f.total() - 3.0).abs() < 1e-12, "total after 2 adds: {}", f.total());
        f.add_rate(3.0);
        assert!((f.total() - 6.0).abs() < 1e-12, "total after 3 adds: {}", f.total());
        f.add_rate(4.0);
        assert!((f.total() - 10.0).abs() < 1e-12, "total after 4 adds: {}", f.total());
        f.add_rate(5.0);
        assert!((f.total() - 15.0).abs() < 1e-12, "total after 5 adds: {}", f.total());
    }

    #[test]
    fn test_add_only_select() {
        let mut f = RateCatalog::new();
        // rates: [1, 2, 3]  total=6
        f.add_rate(1.0);
        f.add_rate(2.0);
        f.add_rate(3.0);

        // u=0.0 → target=0.0, should select index 0
        assert_eq!(f.select(0.0), 0, "select(0.0)");
        // u just past 1/6 → target just past 1.0, should select index 1
        assert_eq!(f.select(0.2), 1, "select(0.2)");
        // u=0.5 → target=3.0, prefix[0..=1]=3.0 → index 1
        assert_eq!(f.select(0.5), 1, "select(0.5)");
        // u=0.99 → target=5.94 → index 2
        assert_eq!(f.select(0.99), 2, "select(0.99)");
    }

    #[test]
    fn test_set_rate_total() {
        let mut f = RateCatalog::new();
        f.add_rate(1.0);
        f.add_rate(2.0);
        f.add_rate(3.0);
        // Change middle rate: [1, 5, 3] total=9
        f.set_rate(1, 5.0);
        assert!((f.total() - 9.0).abs() < 1e-12, "total after set_rate: {}", f.total());
        // Change first rate: [0, 5, 3] total=8
        f.set_rate(0, 0.0);
        assert!((f.total() - 8.0).abs() < 1e-12, "total after zeroing: {}", f.total());
    }

    #[test]
    fn test_swap_remove_last() {
        let mut f = RateCatalog::new();
        f.add_rate(1.0);
        f.add_rate(2.0);
        f.add_rate(3.0);
        // Remove last element: [1, 2] total=3
        let moved = f.swap_remove_rate(2);
        assert!(moved.is_none(), "removing last should return None");
        assert!((f.total() - 3.0).abs() < 1e-12, "total after removing last: {}", f.total());
        assert!((f.get_rate(0) - 1.0).abs() < 1e-12);
        assert!((f.get_rate(1) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_swap_remove_first() {
        let mut f = RateCatalog::new();
        f.add_rate(1.0);
        f.add_rate(2.0);
        f.add_rate(3.0);
        // Remove slot 0: last (3.0) moves into slot 0 → [3, 2] total=5
        let moved = f.swap_remove_rate(0);
        assert_eq!(moved, Some(2), "should report old index of moved element");
        assert!((f.total() - 5.0).abs() < 1e-12, "total after swap_remove(0): {}", f.total());
        assert!((f.get_rate(0) - 3.0).abs() < 1e-12, "slot 0 should now be 3.0, got {}", f.get_rate(0));
        assert!((f.get_rate(1) - 2.0).abs() < 1e-12, "slot 1 should still be 2.0, got {}", f.get_rate(1));
    }

    #[test]
    fn test_swap_remove_middle() {
        let mut f = RateCatalog::new();
        f.add_rate(10.0);
        f.add_rate(20.0);
        f.add_rate(30.0);
        f.add_rate(40.0);
        // Remove slot 1: last (40.0) moves into slot 1 → [10, 40, 30] total=80
        f.swap_remove_rate(1);
        assert!((f.total() - 80.0).abs() < 1e-12, "total after swap_remove(1): {}", f.total());
        assert!((f.get_rate(0) - 10.0).abs() < 1e-12);
        assert!((f.get_rate(1) - 40.0).abs() < 1e-12, "slot 1 should be 40.0, got {}", f.get_rate(1));
        assert!((f.get_rate(2) - 30.0).abs() < 1e-12);
    }

    #[test]
    fn test_swap_remove_then_select() {
        let mut f = RateCatalog::new();
        f.add_rate(1.0);
        f.add_rate(2.0);
        f.add_rate(3.0);
        // Remove slot 0 → [3, 2] total=5
        f.swap_remove_rate(0);
        // u=0.0 → index 0 (rate 3.0)
        assert_eq!(f.select(0.0), 0);
        // u=0.7 → target=3.5 → index 1 (prefix at 0 = 3.0 < 3.5)
        assert_eq!(f.select(0.7), 1);
    }

    #[test]
    fn test_add_after_remove() {
        let mut f = RateCatalog::new();
        f.add_rate(1.0);
        f.add_rate(2.0);
        f.swap_remove_rate(0); // → [2.0] total=2
        assert!((f.total() - 2.0).abs() < 1e-12, "total after remove: {}", f.total());
        f.add_rate(5.0); // → [2.0, 5.0] total=7
        assert!((f.total() - 7.0).abs() < 1e-12, "total after re-add: {}", f.total());
        assert!((f.get_rate(0) - 2.0).abs() < 1e-12);
        assert!((f.get_rate(1) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_single_element_remove() {
        let mut f = RateCatalog::new();
        f.add_rate(42.0);
        f.swap_remove_rate(0);
        assert!((f.total()).abs() < 1e-12, "total should be 0 after removing only element");
    }
}
