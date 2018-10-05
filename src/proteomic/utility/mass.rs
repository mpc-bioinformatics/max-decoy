pub struct NeutralLoss {
    name: &'static str,
    mono_mass: f32,
    average_mass: f32
}

impl NeutralLoss {
    // (name, mono_mass. average_mass)
    pub fn new(mass_tupel: (&'static str, f32, f32)) -> NeutralLoss {
        return NeutralLoss {
            name: mass_tupel.0,
            mono_mass: mass_tupel.1,
            average_mass: mass_tupel.2
        }
    }

    pub fn get_name (&self) -> &'static str {
        return self.name;
    }
    pub fn get_mono_mass (&self) -> f32 {
        return self.mono_mass;
    }
    pub fn get_average_mass (&self) -> f32 {
        return self.average_mass;
    }
}

const WATER_LOSS: (&'static str, f32, f32) = ("H2O", 18.010565, 18.015);
const NONE_LOSS: (&'static str, f32, f32) = ("NONE", 0.0, 0.0);


pub fn get_neutral_loss(name: &str) -> NeutralLoss {
    match name {
        "H2O" => NeutralLoss::new(WATER_LOSS),
        _ => return NeutralLoss::new(NONE_LOSS)
    }
}