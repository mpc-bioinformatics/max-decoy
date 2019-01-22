pub mod mz_ml_reader;
pub mod spectrum;

use std::collections::HashMap;

const INDENTION: &str = "    ";

fn parse_tag_attributes(tag: &quick_xml::events::BytesStart) -> Box<HashMap<String, String>> {
    let mut attributes_map: HashMap<String, String> = HashMap::new();
    for attribute_result in tag.attributes().into_iter() {
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

fn indent(level: usize) -> String {
    return (0..level).map(|_| INDENTION).collect::<String>();
}