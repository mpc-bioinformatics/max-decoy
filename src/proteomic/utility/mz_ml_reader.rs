extern crate quick_xml;

use std::path::Path;
use std::collections::HashMap;

use self::quick_xml::Reader;
use self::quick_xml::events::Event;

use proteomic::models::mass;

const NAME_OF_MASS_TO_CHARGE_CV_PARAM: &str = "selected ion m/z";
const NAME_OF_CHARGE_CV_PARAM: &str = "charge state";

pub struct MzMlReader {
    file_path: String
}

impl MzMlReader {
    pub fn new(file_path: &str) -> Self {
        return Self {
            file_path: file_path.to_owned()
        };
    }

    pub fn get_precursor_masses(&self) -> Box<Vec<f64>> {
        let mut precursor_masses: Vec<f64> = Vec::new();
        let mut reader = match Reader::from_file(Path::new(self.file_path.as_str())) {
            Ok(reader) => reader,
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_precursor_masses(): error when reading MzMl-File, original error: {:?}", err)
        };
        let mut buf = Vec::new();
        let mut inside_selected_ion: bool = false;
        let mut charge: Option<u8> = None;
        let mut mass_to_charge_ratio: Option<f64> = None;
        loop {
            match reader.read_event(&mut buf) {
                // entering tags: <tag (attr1="" attr2="" ...)>
                Ok(Event::Start(ref tag)) => {
                    match tag.name() {
                        b"selectedIon" => inside_selected_ion = true,
                        _ => (),
                    }
                },
                // empty tags: <tag (attr1="" attr2="" ...)/>
                Ok(Event::Empty(ref tag)) => {
                    match tag.name() {
                        // process cvParam only if it is inside SelectedIon
                        b"cvParam" if inside_selected_ion => {
                            let attributes = Self::parse_event_attributes(&mut tag.attributes());
                            if let Some(cv_param_type) = attributes.get("name") {
                                // check which name the cvParam has
                                match cv_param_type.as_str() {
                                    NAME_OF_MASS_TO_CHARGE_CV_PARAM => mass_to_charge_ratio = Some(Self::get_mass_to_charge_ratio_from_attributes(&attributes)),
                                    NAME_OF_CHARGE_CV_PARAM => charge = Some(Self::get_charge_from_attributes(&attributes)),
                                    _ => () // other cvParams than "selected ion m/z" and "charge state" are uninteresting
                                }
                            }
                        }
                        _ => ()
                    }
                }
                // leaving of tags: </tag>
                Ok(Event::End(ref tag)) => {
                    match tag.name() {
                        b"selectedIon" => {
                            inside_selected_ion = false;
                            if charge.is_some() & mass_to_charge_ratio.is_some() {
                                precursor_masses.push(
                                    mass::thomson_to_dalton(mass_to_charge_ratio.unwrap(), charge.unwrap())
                                );
                                charge = None;
                                mass_to_charge_ratio = None;
                            } else {
                                panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_precursor_masses(): Leaving selectedIon but precursor mass was not calculated. Either charge or mass_to_charge ratios was missing.");
                            }
                        },
                        _ => (),
                    }
                },
                Ok(Event::Eof) => break, // exits the loop when reaching end of file
                Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
                _ => (), // There are several other `Event`s we do not consider here
            }

            // if we don't keep a borrow elsewhere, we can clear the buffer to keep memory usage low
            buf.clear();
        }
        return Box::new(precursor_masses);
    }

    fn parse_event_attributes(attributes: &mut quick_xml::events::attributes::Attributes) -> Box<HashMap<String, String>> {
        let mut attributes_map: HashMap<String, String> = HashMap::new();
        for attribute_result in attributes.into_iter() {
            if let Ok(attribute) = attribute_result {
                attributes_map.insert(
                    match std::str::from_utf8(attribute.key) {
                        Ok(key) => key.to_owned(),
                        Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.parse_event_attributes(): Error at parsing the attribute's key: {}", err)
                    },
                    match attribute.unescaped_value() {
                        Ok(value) => match std::str::from_utf8(value.as_ref()) {
                            Ok(key) => key.to_owned(),
                            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.parse_event_attributes(): Error at parsing the attribute's key: {}", err)
                        },
                        Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.parse_event_attributes(): Error at parsing the attribute's value: {}", err)
                    }
                );
            }
        }
        return Box::new(attributes_map);
    }

    fn get_charge_from_attributes(attributes: &HashMap<String, String>) -> u8 {
        match attributes.get("value") {
            Some(value) => match value.parse::<u8>() {
                Ok(charge) => return charge,
                Err(_) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_charge_from_attributes(): cvParam name is '{}' but could not parse charge unsigned integet (8 bit)", NAME_OF_CHARGE_CV_PARAM)
            },
            None => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_charge_from_attributes(): cvParam name is '{}' there is no attribute 'value'", NAME_OF_CHARGE_CV_PARAM)
        }
    }

    fn get_mass_to_charge_ratio_from_attributes(attributes: &HashMap<String, String>) -> f64 {
        match attributes.get("value") {
            Some(value) => match value.parse::<f64>() {
                Ok(charge) => return charge,
                Err(_) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_charge_from_attributes(): cvParam name is '{}' but could not parse mass to charge ratio to float (64 bit)", NAME_OF_MASS_TO_CHARGE_CV_PARAM)
            },
            None => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_charge_from_attributes(): cvParam name is '{}' there is no attribute 'value'", NAME_OF_MASS_TO_CHARGE_CV_PARAM)
        }
    }
}