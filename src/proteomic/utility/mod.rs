pub mod input_file_digester;
pub mod database_connection;
pub mod logger;
pub mod combinations;
pub mod decoy_generator;
pub mod mz_ml;

pub fn parts_per_million_of(value: f64, ppm: i64) -> f64 {
    return  value / 1000000.0 * ppm as f64;
}