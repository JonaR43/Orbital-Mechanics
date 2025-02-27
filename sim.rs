extern crate nalgebra as na;

use std::f64::consts::PI;

// Constants
const G: f64 = 6.67430e-11; // Gravitational constant (m^3 kg^-1 s^-2)
const M: f64 = 5.972e24; // Mass of the Earth (kg)
const R_E: f64 = 6.3781e6; // Earth's radius (m)
const J2: f64 = 1.08262668e-3; // J2 coefficient (dimensionless)

// Orbital Elements Structure
struct OrbitalElements {
    semi_major_axis: f64, // in meters
    inclination: f64,     // in radians
    omega: f64,           // argument of perigee in radians
    ascending_node: f64,  // in radians
}

impl OrbitalElements {
    // Method to update the orbital elements based on J2 perturbations
    fn update_elements(&mut self, dt: f64) {
        let a = self.semi_major_axis;
        let i = self.inclination;
        let omega = self.omega;
        let omega_dot = self.nodal_precession_rate(a, i);
        let perigee_dot = self.argument_of_perigee_precession_rate(a, i);

        // Update the orbital elements over time (assuming small changes)
        self.ascending_node += omega_dot * dt;
        self.omega += perigee_dot * dt;

        // Normalize angles to [0, 2*pi]
        self.ascending_node = self.normalize_angle(self.ascending_node);
        self.omega = self.normalize_angle(self.omega);
    }

    // Normalize an angle to the range [0, 2*pi]
    fn normalize_angle(&self, angle: f64) -> f64 {
        let mut angle = angle % (2.0 * PI);
        if angle < 0.0 {
            angle += 2.0 * PI;
        }
        angle
    }

    // Calculate the rate of nodal precession (dΩ/dt)
    fn nodal_precession_rate(&self, a: f64, i: f64) -> f64 {
        (-3.0 / 2.0) * J2 * (R_E.powi(2)) / a.powi(2) * (1.0 - 0.0) * i.cos()
    }

    // Calculate the rate of argument of perigee precession (dω/dt)
    fn argument_of_perigee_precession_rate(&self, a: f64, i: f64) -> f64 {
        (3.0 / 4.0) * J2 * (R_E.powi(2)) / a.powi(2) * (4.0 - 5.0 * i.sin().powi(2))
    }
}

fn main() {
    // Define initial orbital elements
    let mut orbital_elements = OrbitalElements {
        semi_major_axis: 7000.0e3, // Example: 7000 km altitude circular orbit
        inclination: 0.0,          // Example: Equatorial orbit
        omega: 0.0,                // Initial argument of perigee
        ascending_node: 0.0,       // Initial ascending node
    };

    // Time step (in seconds) for the simulation
    let dt = 3600.0; // 1 hour

    // Run the simulation for 24 hours
    for t in 0..24 {
        orbital_elements.update_elements(dt);
        println!(
            "Time: {} hours => Ascending Node: {:.6} rad, Argument of Perigee: {:.6} rad",
            t, orbital_elements.ascending_node, orbital_elements.omega
        );
    }
}
