use proteomic::models::mass::convert_mass_to_int;

const WATER_LOSS: (&'static str, f64, f64) = ("H2O", 18.010565, 18.015);
const NONE_LOSS: (&'static str, f64, f64) = ("NONE", 0.0, 0.0);

pub struct NeutralLoss {
    name: &'static str,
    mono_mass: i64,
    average_mass: i64
}

impl NeutralLoss {
    // (name, mono_mass. average_mass)
    pub fn new(mass_tupel: (&'static str, f64, f64)) -> NeutralLoss {
        return NeutralLoss {
            name: mass_tupel.0,
            mono_mass: convert_mass_to_int(mass_tupel.1),
            average_mass: convert_mass_to_int(mass_tupel.2)
        }
    }

    pub fn get_name (&self) -> &'static str {
        return self.name;
    }

    pub fn get_mono_mass (&self) -> i64 {
        return self.mono_mass;
    }

    pub fn get_average_mass (&self) -> i64 {
        return self.average_mass;
    }

    pub fn get(name: &str) -> NeutralLoss {
        match name {
            "H2O" => NeutralLoss::new(WATER_LOSS),
            _ => return NeutralLoss::new(NONE_LOSS)
        }
    }
}




