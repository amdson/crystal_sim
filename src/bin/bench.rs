//! Headless benchmark: run the simulation for a fixed wall-clock duration and
//! report throughput and per-stage timing.
//!
//! Usage:
//!   cargo run --bin bench --release -- [config.json] [--steps N] [--seconds S]
//!
//! Defaults: config/checkerboard.json, run for 5 seconds.

use std::time::Instant;
use crystal_sim::config::SimConfig;
use crystal_sim::kmc::Simulation;

fn main() {
    let mut config_path = "config/checkerboard.json".to_string();
    let mut target_steps: Option<u64> = None;
    let mut target_secs: f64 = 5.0;

    // Minimal arg parsing.
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--steps" => {
                i += 1;
                target_steps = Some(args[i].parse().expect("--steps expects an integer"));
            }
            "--seconds" => {
                i += 1;
                target_secs = args[i].parse().expect("--seconds expects a number");
            }
            other => config_path = other.to_string(),
        }
        i += 1;
    }

    let json = std::fs::read_to_string(&config_path)
        .unwrap_or_else(|e| panic!("Cannot read {config_path}: {e}"));
    let mut config: SimConfig = serde_json::from_str(&json)
        .unwrap_or_else(|e| panic!("Invalid config: {e}"));
    config.init_cache();

    println!("Config:      {config_path}");
    println!("Particle types: {}", config.n_types());

    // Warmup: 1 second, not counted.
    let chunk = config.steps_per_frame.max(1);
    let mut sim = Simulation::new(config.clone());

    print!("Warming up...");
    let warmup = Instant::now();
    while warmup.elapsed().as_secs_f64() < 1.0 {
        sim.step(chunk);
    }
    println!(" done ({} particles after warmup)", sim.particle_count());

    // Timed run.
    sim = Simulation::new(config.clone());
    let start = Instant::now();
    let mut total_steps: u64 = 0;

    loop {
        sim.step(chunk);
        total_steps += chunk as u64;

        let elapsed = start.elapsed().as_secs_f64();
        if let Some(s) = target_steps {
            if total_steps >= s { break; }
        } else if elapsed >= target_secs {
            break;
        }
    }

    let elapsed = start.elapsed().as_secs_f64();
    let steps_per_sec = total_steps as f64 / elapsed;

    println!();
    println!("── Results ──────────────────────────────────────────");
    println!("Elapsed:     {elapsed:.2}s");
    println!("Steps:       {total_steps}");
    println!("Throughput:  {:.0} steps/s", steps_per_sec);
    println!("Particles:   {}", sim.particle_count());
    println!("Sim time:    {:.3}", sim.simulation_time());
    println!();
    println!("── Per-stage timing (entire run) ────────────────────");
    sim.print_timers();
}
