use quick_xml::Reader;
use quick_xml::events::Event;

use proteomic::utility::mz_ml;

pub struct Chromatrogram {
    id_ref: String,
    xml: String,
    indent_level: usize
}

impl Chromatrogram {
    pub fn new(chromatogram_xml: &str, indent_level: usize) -> Self {
        let mut reader = Reader::from_str(chromatogram_xml);
        let mut buf = Vec::new();
        let mut id_ref: Option<String> = None;
        loop {
            match reader.read_event(&mut buf) {
                // entering tags: <tag (attr1="" attr2="" ...)>
                Ok(Event::Start(ref tag)) => {
                    match tag.name() {
                        b"chromatogram" => {
                            let attributes = mz_ml::parse_tag_attributes(tag);
                            if let Some(id_value) = attributes.get("id") {
                                id_ref = Some(id_value.clone());
                                break;
                            }
                        }
                        _ => (),
                    }
                },
                Ok(Event::Eof) => {
                    if id_ref.is_some() {
                        break;
                    } else {
                        panic!("proteomic::utility::mz_ml::chromatogram::Chromatogram::new(): Reached end of chromatogram xml but id-ref is missing:\n{}", chromatogram_xml);
                    }
                }, // exits the loop when reaching end of file
                Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
                _ => (), // There are several other `Event`s we do not consider here
            }
            // if we don't keep a borrow elsewhere, we can clear the buffer to keep memory usage low
            buf.clear();
        }
        return Self {
            id_ref: id_ref.unwrap(),
            xml: chromatogram_xml.to_owned(),
            indent_level: indent_level
        };
    }

    pub fn len(&self) -> usize {
        return self.xml.len();
    }

    pub fn get_id_ref(&self) -> &str {
        return self.id_ref.as_str();
    }

    pub fn get_xml(&self) -> &str {
        return self.xml.as_str();
    }

    pub fn get_indent_level(&self) -> usize {
        return self.indent_level;
    }
}