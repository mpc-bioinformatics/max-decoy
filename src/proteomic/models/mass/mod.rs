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

pub fn calculate_upper_and_lower_tolerance(precursor_mass: i64, upper_mass_limit_ppm: i64, lower_mass_limit_ppm: i64) -> (i64, i64) {
    return (
        (precursor_mass as f64 / 1000000.0 * lower_mass_limit_ppm as f64) as i64,
        (precursor_mass as f64 / 1000000.0 * upper_mass_limit_ppm as f64) as i64
    )
}

/// Returns a tuple like (lower_precursor_tolerance_limit, upper_precursor_tolerance_limit)
pub fn calculate_precursor_tolerance(precursor_mass: i64, upper_mass_limit_ppm: i64, lower_mass_limit_ppm: i64) -> (i64, i64) {
    let tolerances = calculate_upper_and_lower_tolerance(precursor_mass, upper_mass_limit_ppm, lower_mass_limit_ppm);
    return (
        precursor_mass - tolerances.0,
        precursor_mass + tolerances.1
    )
}