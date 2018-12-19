pub mod neutral_loss;

const MASS_CONVERT_FACTOR: f64 = 1000000.0;

pub fn convert_mass_to_int(mass: f64) -> i64 {
    return (mass * MASS_CONVERT_FACTOR) as i64;
}

pub fn convert_mass_to_float(mass: i64) -> f64 {
    return mass as f64 / MASS_CONVERT_FACTOR;
}