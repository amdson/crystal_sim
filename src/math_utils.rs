//! Mathematical utilities for the crystal simulation.
//!
//! This module provides common mathematical functions and constants
//! used throughout the codebase.

/// Mathematical constants
pub mod constants {
    pub const PI: f64 = std::f64::consts::PI;
    pub const TAU: f64 = std::f64::consts::TAU;
    pub const E: f64 = std::f64::consts::E;
}

/// Vector and matrix operations
pub mod vectors {
    // Vector functions will be implemented here
}

/// Geometry calculations
pub mod geometry {
    // Geometry functions will be implemented here
}

/// Utility functions for common mathematical operations
pub mod utils {
    /// Clamps a value between a minimum and maximum
    pub fn clamp(value: f64, min: f64, max: f64) -> f64 {
        value.max(min).min(max)
    }

    /// Checks if two floating point numbers are approximately equal
    pub fn approx_equal(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    /// Computes the exponential of a value using a third-order approximation
    pub fn exp3(x: f64) -> f64 {
        (120.0 + x * (120.0 + x * (60.0 + x * (20.0 + x * (5.0 + x))))) * 0.0083333333
    }
}
