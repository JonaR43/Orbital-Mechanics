use nalgebra::{Vector3, Point3};
use std::f64::consts::PI;

const G: f64 = 6.67430e-11; // Gravitational constant (m^3·kg^-1·s^-2)
const M: f64 = 5.972e24; // Mass of Earth (kg)
const R_EARTH: f64 = 6.3781e6; // Radius of Earth (m)
const J2: f64 = 1.08263e-3; // Second zonal harmonic (J2 term for Earth)

#[derive(Debug)]
struct Satellite {
    position: Point3<f64>,
    velocity: Vector3<f64>,
    mass: f64,
    semi_major_axis: f64,
    eccentricity: f64,
    inclination: f64,
}

impl Satellite {
    fn new(semi_major_axis: f64, eccentricity: f64, inclination: f64, true_anomaly: f64) -> Self {
        let mean_motion = self.calculate_mean_motion(semi_major_axis);
        let velocity = self.calculate_velocity(semi_major_axis, true_anomaly);
        let position = self.calculate_position(semi_major_axis, eccentricity, true_anomaly);
        
        Satellite {
            position,
            velocity,
            mass: 500.0, // Example mass in kg
            semi_major_axis,
            eccentricity,
            inclination,
        }
    }

    fn calculate_mean_motion(&self, semi_major_axis: f64) -> f64 {
        let period = 2.0 * PI * (semi_major_axis.powf(3.0) / G / M).sqrt();
        (2.0 * PI) / period
    }

    fn calculate_velocity(&self, semi_major_axis: f64, true_anomaly: f64) -> Vector3<f64> {
        let distance = semi_major_axis * (1.0 - 0.01 * true_anomaly.cos()); // Simplified approximation
        let velocity_magnitude = (G * M / distance).sqrt();
        let velocity_direction = Vector3::new(true_anomaly.cos(), true_anomaly.sin(), 0.0);
        velocity_direction * velocity_magnitude
    }

    fn calculate_position(&self, semi_major_axis: f64, eccentricity: f64, true_anomaly: f64) -> Point3<f64> {
        let distance = semi_major_axis * (1.0 - eccentricity * true_anomaly.cos());
        Point3::new(distance * true_anomaly.cos(), distance * true_anomaly.sin(), 0.0)
    }

    // J2 Perturbation effect (simplified model for an equatorial orbit)
    fn apply_j2_effect(&mut self, dt: f64) {
        let r = self.position.norm(); // Radial distance
        let theta = 0.0; // For simplicity, assuming equatorial orbit

        // J2 perturbation acceleration formula
        let j2_acceleration_magnitude = (3.0 / 2.0) * J2 * (R_EARTH.powi(2)) / r.powi(4)
            * (5.0 * theta.sin().powi(2) - 1.0);

        // Perturbing force direction
        let perturbation = Vector3::new(0.0, 0.0, j2_acceleration_magnitude); // Z-axis perturbation for equatorial orbit

        // Update velocity due to J2 effect
        self.velocity += perturbation * dt;
    }

    fn update_position(&mut self, dt: f64) {
        // Update position using velocity
        self.position += Point3::from(self.velocity * dt);
    }

    fn print_status(&self) {
        println!("Position: {:?}", self.position);
        println!("Velocity: {:?}", self.velocity);
    }
}

fn main() {
    let semi_major_axis = 7.0e6; // 7000 km orbit (in meters)
    let eccentricity = 0.01; // Low eccentricity for a near-circular orbit
    let inclination = 0.0; // Equatorial orbit (0° inclination)
    let true_anomaly = 0.0; // Assume the satellite is at the periapsis (0°)

    let mut satellite = Satellite::new(semi_major_axis, eccentricity, inclination, true_anomaly);

    let dt = 1000.0; // Time step in seconds (1 second per update)

    // Simulate for 10 orbits
    let total_orbit_time = 2.0 * PI * (semi_major_axis.powf(3.0) / G / M).sqrt();
    let num_steps = (total_orbit_time / dt) as usize;

    for _ in 0..num_steps {
        satellite.apply_j2_effect(dt); // Apply J2 effect
        satellite.update_position(dt); // Update position
        satellite.print_status(); // Output current state
    }
}
