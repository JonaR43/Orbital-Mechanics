extern crate nalgebra as na;

use std::f64::consts::PI;

// Constants
const G: f64 = 6.67430e-11; // Gravitational constant (m^3 kg^-1 s^-2)
const M: f64 = 5.972e24; // Mass of the Earth (kg)
const R_E: f64 = 6.3781e6; // Earth's radius (m)
const J2: f64 = 1.08262668e-3; // J2 coefficient (dimensionless)
const J3: f64 = -2.532e-6; // J3 coefficient (dimensionless)
const SOLAR_CONSTANT: f64 = 1361.0; // Solar constant (W/m^2)
const ALBEDO: f64 = 0.3; // Albedo of the satellite (reflectivity)
const CROSS_SECTIONAL_AREA: f64 = 10.0; // Cross-sectional area of the satellite (m^2)
const MASS: f64 = 500.0; // Mass of the satellite (kg)
const SOLAR_PRESSURE: f64 = 4.5e-6; // Solar radiation pressure constant (N/m^2)

// Orbital Elements Structure
struct OrbitalElements {
    semi_major_axis: f64, // in meters
    eccentricity: f64,     // dimensionless
    inclination: f64,      // in radians
    omega: f64,            // argument of perigee in radians
    ascending_node: f64,   // in radians
    true_anomaly: f64,     // in radians
}

impl OrbitalElements {
    // Method to update the orbital elements based on perturbations
    fn update_elements(&mut self, dt: f64) {
        let a = self.semi_major_axis;
        let e = self.eccentricity;
        let i = self.inclination;
        let omega = self.omega;
        let omega_dot = self.nodal_precession_rate(a, i);
        let perigee_dot = self.argument_of_perigee_precession_rate(a, i);
        let eccentricity_dot = self.eccentricity_precession_rate(a, e, i);

        // Update orbital elements
        self.ascending_node += omega_dot * dt;
        self.omega += perigee_dot * dt;
        self.eccentricity += eccentricity_dot * dt;

        // Normalize angles to [0, 2*pi]
        self.ascending_node = self.normalize_angle(self.ascending_node);
        self.omega = self.normalize_angle(self.omega);
        self.true_anomaly += self.true_anomaly_dot(a, e) * dt;

        // Normalize true anomaly
        self.true_anomaly = self.normalize_angle(self.true_anomaly);
    }

    // Normalize an angle to the range [0, 2*pi]
    fn normalize_angle(&self, angle: f64) -> f64 {
        let mut angle = angle % (2.0 * PI);
        if angle < 0.0 {
            angle += 2.0 * PI;
        }
        angle
    }

    // Calculate the rate of nodal precession (dΩ/dt) due to J2 and J3
    fn nodal_precession_rate(&self, a: f64, i: f64) -> f64 {
        let j2_effect = (-3.0 / 2.0) * J2 * (R_E.powi(2)) / a.powi(2) * (1.0 - e * e) * i.cos();
        let j3_effect = (-5.0 / 4.0) * J3 * (R_E.powi(3)) / a.powi(3) * (1.0 - e * e) * i.cos();
        j2_effect + j3_effect
    }

    // Calculate the rate of argument of perigee precession (dω/dt) due to J2
    fn argument_of_perigee_precession_rate(&self, a: f64, i: f64) -> f64 {
        (3.0 / 4.0) * J2 * (R_E.powi(2)) / a.powi(2) * (4.0 - 5.0 * i.sin().powi(2))
    }

    // Calculate the rate of eccentricity precession (de/dt) due to solar radiation pressure
    fn eccentricity_precession_rate(&self, a: f64, e: f64, i: f64) -> f64 {
        // Simplified model of solar radiation pressure effect
        let solar_pressure_effect = SOLAR_PRESSURE * CROSS_SECTIONAL_AREA * ALBEDO / (2.0 * a.powi(2));
        solar_pressure_effect * (1.0 - e)
    }

    // Calculate the rate of true anomaly (dθ/dt) based on the orbital parameters
    fn true_anomaly_dot(&self, a: f64, e: f64) -> f64 {
        let n = (G * M / a.powi(3)).sqrt(); // Mean motion
        n * (1.0 - e.powi(2)) / (1.0 + e * (self.true_anomaly).cos()) // True anomaly rate
    }

    // Convert orbital elements to position in 3D space (simplified version)
    fn calculate_position(&self) -> na::Vector3<f64> {
        let a = self.semi_major_axis;
        let e = self.eccentricity;
        let i = self.inclination;
        let omega = self.omega;
        let omega_node = self.ascending_node;
        let theta = self.true_anomaly;

        // Calculate position in orbital plane (r, theta) in polar coordinates
        let r = a * (1.0 - e.powi(2)) / (1.0 + e * theta.cos());

        // Convert to 3D coordinates in the inertial frame
        let x = r * (omega + omega_node).cos();
        let y = r * (omega + omega_node).sin();
        let z = r * i.sin();

        na::Vector3::new(x, y, z)
    }
}

fn main() {
    // Define initial orbital elements
    let mut orbital_elements = OrbitalElements {
        semi_major_axis: 7000.0e3, // Example: 7000 km altitude circular orbit
        eccentricity: 0.001,       // Slightly elliptical orbit
        inclination: 0.0,          // Example: Equatorial orbit
        omega: 0.0,                // Initial argument of perigee
        ascending_node: 0.0,       // Initial ascending node
        true_anomaly: 0.0,         // Initial true anomaly
    };

    // Time step (in seconds) for the simulation
    let dt = 3600.0; // 1 hour

    // Run the simulation for 24 hours
    for t in 0..24 {
        orbital_elements.update_elements(dt);
        let position = orbital_elements.calculate_position();
        println!(
            "Time: {} hours => Position: ({:.3}, {:.3}, {:.3}) km",
            t,
            position.x / 1000.0, // Convert meters to kilometers
            position.y / 1000.0,
            position.z / 1000.0
        );
    }
}

