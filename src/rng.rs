/// splitmix64 PRNG — no OS entropy, no getrandom, safe for WASM.
pub struct Rng {
    state: u64,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self { state: if seed == 0 { 0x9e3779b97f4a7c15 } else { seed } }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }

    /// Uniform [0, 1)
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform integer in [0, n)
    pub fn next_usize(&mut self, n: usize) -> usize {
        if n <= 1 { return 0; }
        (self.next_f64() * n as f64) as usize
    }
}
