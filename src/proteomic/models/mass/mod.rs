pub mod neutral_loss;

const MASS_CONVERT_FACTOR: f64 = 1000000.0;
const HYDROGEN_MONO_MASS: f64 = 1.007276;

pub fn convert_mass_to_int(mass: f64) -> i64 {
    return (mass * MASS_CONVERT_FACTOR) as i64;
}

pub fn convert_mass_to_float(mass: i64) -> f64 {
    return mass as f64 / MASS_CONVERT_FACTOR;
}

pub fn thomson_to_dalton(thomson: f64, charge: u8) -> f64 {
    return thomson * charge as f64 - HYDROGEN_MONO_MASS * charge as f64;
}