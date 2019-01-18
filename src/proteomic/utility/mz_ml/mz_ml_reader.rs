use std::path::Path;

use quick_xml::Reader;
use quick_xml::events::Event;

use proteomic::utility::mz_ml;
use proteomic::utility::mz_ml::spectrum::Spectrum;
use proteomic::utility::mz_ml::chromatogram::Chromatrogram;

const NAME_OF_MS_LEVEL_CV_PARAM: &str = "ms level";

pub struct MzMlReader {
    file_path: String
}

impl MzMlReader {
    pub fn new(file_path: &str) -> Self {
        return Self {
            file_path: file_path.to_owned()
        };
    }

    pub fn get_ms_two_spectra(&self) -> Box<Vec<Spectrum>> {
        let mut spectra: Vec<Spectrum> = Vec::new();
        let mut reader = match Reader::from_file(Path::new(self.file_path.as_str())) {
            Ok(reader) => reader,
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_precursor_masses(): error when reading MzMl-File, original error: {:?}", err)
        };
        let mut buf = Vec::new();
        let mut inside_spectrum: bool = false;
        let mut indent_level = 0;
        let mut spectrum_lines: Vec<String> = Vec::new();
        loop {
            match reader.read_event(&mut buf) {
                // entering tags: <tag (attr1="" attr2="" ...)>
                Ok(Event::Start(ref tag)) => {
                    match tag.name() {
                        b"spectrum" => {
                            spectrum_lines.push(format!("{}<{}>", mz_ml::indent(indent_level), Self::bytes_start_to_string(tag)));
                            inside_spectrum = true;
                        },
                        _ if inside_spectrum => spectrum_lines.push(format!("{}<{}>", mz_ml::indent(indent_level), Self::bytes_start_to_string(tag))),
                        _ => (),
                    }
                    indent_level += 1;
                },
                // empty tags: <tag (attr1="" attr2="" ...)/>
                Ok(Event::Empty(ref tag)) => {
                    match tag.name() {
                        _ if inside_spectrum => spectrum_lines.push(format!("{}<{}/>", mz_ml::indent(indent_level), Self::bytes_start_to_string(tag))),
                        _ => ()
                    }
                }
                // leaving of tags: </tag>
                Ok(Event::End(ref tag)) => {
                    indent_level -= 1;
                    match tag.name() {
                        b"spectrum" => {
                            spectrum_lines.push(format!("{}</{}>", mz_ml::indent(indent_level), Self::bytes_end_to_string(&tag)));
                            inside_spectrum = false;
                            let spectrum: String = spectrum_lines.join("\n");
                            if Self::is_ms_two_spectrum(spectrum.as_str()) {
                                spectra.push(Spectrum::new(spectrum.as_str(), indent_level));
                            }
                            spectrum_lines.clear();
                        },
                        _ if inside_spectrum => spectrum_lines.push(format!("{}</{}>", mz_ml::indent(indent_level), Self::bytes_end_to_string(&tag))),
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
        return Box::new(spectra);
    }

    pub fn get_content_before_spectrum_list(&self) -> String {
         let mut reader = match Reader::from_file(Path::new(self.file_path.as_str())) {
            Ok(reader) => reader,
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_precursor_masses(): error when reading MzMl-File, original error: {:?}", err)
        };
        let mut buf = Vec::new();
        let mut indent_level = 0;
        let mut before_spectrum_list_lines: Vec<String> = Vec::new();
        loop {
            match reader.read_event(&mut buf) {
                // entering tags: <tag (attr1="" attr2="" ...)>
                Ok(Event::Start(ref tag)) => {
                    match tag.name() {
                       b"spectrumList" => break,
                        _ => before_spectrum_list_lines.push(format!("{}<{}>", mz_ml::indent(indent_level), Self::bytes_start_to_string(tag)))
                    }
                    indent_level += 1;
                },
                // empty tags: <tag (attr1="" attr2="" ...)/>
                Ok(Event::Empty(ref tag)) => {
                    match tag.name() {
                        _ => before_spectrum_list_lines.push(format!("{}<{}/>", mz_ml::indent(indent_level), Self::bytes_start_to_string(tag)))
                    }
                }
                // leaving of tags: </tag>
                Ok(Event::End(ref tag)) => {
                    indent_level -= 1;
                    match tag.name() {
                        _ => before_spectrum_list_lines.push(format!("{}</{}>", mz_ml::indent(indent_level), Self::bytes_end_to_string(&tag)))
                    }
                },
                Ok(Event::Eof) => break, // exits the loop when reaching end of file
                Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
                _ => (), // There are several other `Event`s we do not consider here
            }

            // if we don't keep a borrow elsewhere, we can clear the buffer to keep memory usage low
            buf.clear();
        }
        return before_spectrum_list_lines.join("\n");
    }

    pub fn get_chromatograms(&self) -> Box<Vec<Chromatrogram>> {
        let mut chromatograms: Vec<Chromatrogram> = Vec::new();
         let mut reader = match Reader::from_file(Path::new(self.file_path.as_str())) {
            Ok(reader) => reader,
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_precursor_masses(): error when reading MzMl-File, original error: {:?}", err)
        };
        let mut buf = Vec::new();
        let mut inside_chromatogram: bool = false;
        let mut indent_level = 0;
        let mut chromatogram_list_lines: Vec<String> = Vec::new();
        loop {
            match reader.read_event(&mut buf) {
                // entering tags: <tag (attr1="" attr2="" ...)>
                Ok(Event::Start(ref tag)) => {
                    match tag.name() {
                       b"chromatogram" => {
                           inside_chromatogram = true;
                           chromatogram_list_lines.push(format!("{}<{}>", mz_ml::indent(indent_level), Self::bytes_start_to_string(tag)));
                       },  
                        _ if inside_chromatogram => chromatogram_list_lines.push(format!("{}<{}>", mz_ml::indent(indent_level), Self::bytes_start_to_string(tag))),
                        _ => ()
                    }
                    indent_level += 1;
                },
                // empty tags: <tag (attr1="" attr2="" ...)/>
                Ok(Event::Empty(ref tag)) => {
                    match tag.name() {
                        _ if inside_chromatogram => chromatogram_list_lines.push(format!("{}<{}/>", mz_ml::indent(indent_level), Self::bytes_start_to_string(tag))),
                        _ => ()
                    }
                }
                // leaving of tags: </tag>
                Ok(Event::End(ref tag)) => {
                    indent_level -= 1;
                    match tag.name() {
                        b"chromatogram" => {
                            chromatogram_list_lines.push(format!("{}</{}>", mz_ml::indent(indent_level), Self::bytes_end_to_string(&tag)));
                            chromatograms.push(Chromatrogram::new(chromatogram_list_lines.join("\n").as_str(), indent_level));
                            chromatogram_list_lines.clear();
                        },
                        b"chromatogramList" => break,
                        _ if inside_chromatogram => chromatogram_list_lines.push(format!("{}</{}>", mz_ml::indent(indent_level), Self::bytes_end_to_string(&tag))),
                        _ => ()
                    }
                },
                Ok(Event::Eof) => break, // exits the loop when reaching end of file
                Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
                _ => (), // There are several other `Event`s we do not consider here
            }

            // if we don't keep a borrow elsewhere, we can clear the buffer to keep memory usage low
            buf.clear();
        }
        return Box::new(chromatograms);
    }

    fn bytes_start_to_string(bytes: &quick_xml::events::BytesStart) -> String {
        match bytes.unescaped() {
            Ok(ref bytes_array) => match std::str::from_utf8(bytes_array) {
                Ok(string) => string.to_owned(),
                Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.bytes_start_to_string(): Error at std::str::from_utf8(): {}", err)
            }
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.bytes_start_to_string(): Error at bytes.unescaped(): {}", err)
        }
    }

    fn bytes_end_to_string(bytes: &quick_xml::events::BytesEnd) -> String {
        match std::str::from_utf8(bytes.name()) {
            Ok(string) => string.to_owned(),
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.bytes_end_to_string(): Error at std::str::from_utf8(): {}", err)
        }
    }

    /// Checks if inside the given XML is a tag 'cvParam' with attributes name="ms level" and value="2".
    /// Returns bool.
    ///
    /// # Arguments
    ///
    /// * `spectrum` - A spectrum tag from mzML-file, make sure it is realy a spectrum tag because no validation is done.
    fn is_ms_two_spectrum(spectrum: &str) -> bool {
        let mut reader = Reader::from_str(spectrum);
        let mut buf = Vec::new();
        loop {
            match reader.read_event(&mut buf) {
                // empty tags: <tag (attr1="" attr2="" ...)/>
                Ok(Event::Empty(ref tag)) => {
                    match tag.name() {
                        b"cvParam" => {
                            let attributes = mz_ml::parse_tag_attributes(tag);
                            if let Some(cv_param_type) = attributes.get("name") {
                                if cv_param_type == NAME_OF_MS_LEVEL_CV_PARAM {
                                    if let Some(ms_level) = attributes.get("value") {
                                        return ms_level.eq("2"); // without further conversion, attributes are Strings
                                    }
                                }
                            }
                        },
                        _ => ()
                    }
                }
                Ok(Event::Eof) => break, // exits the loop when reaching end of file
                Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
                _ => (), // There are several other `Event`s we do not consider here
            }

            // if we don't keep a borrow elsewhere, we can clear the buffer to keep memory usage low
            buf.clear();
        }
        return false;
    }
}