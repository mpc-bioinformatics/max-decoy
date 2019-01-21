use std::path::Path;
use std::io::Cursor;

use quick_xml::Reader;
use quick_xml::events::Event;
use quick_xml::Writer;

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
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_ms_two_spectra(): error when reading MzMl-File, original error: {:?}", err)
        };
        let mut buf = Vec::new();
        let mut inside_spectrum: bool = false;
        let mut indent_level = 0;
        let mut spectrum_writer = Writer::new(Cursor::new(Vec::new())); // chromatogram list
        loop {
            match reader.read_event(&mut buf) {
                // entering tags: <tag (attr1="" attr2="" ...)>
                Ok(Event::Start(ref tag)) => {
                    match tag.name() {
                        b"spectrum" => {
                            inside_spectrum = true;
                            Self::write_event(&mut spectrum_writer, &Event::Start(tag.to_owned()), indent_level);
                        },
                        _ if inside_spectrum => Self::write_event(&mut spectrum_writer, &Event::Start(tag.to_owned()), indent_level),
                        _ => (),
                    }
                    indent_level += 1;
                },
                // empty tags: <tag (attr1="" attr2="" ...)/>
                Ok(Event::Empty(ref tag)) => {
                    match tag.name() {
                        _ if inside_spectrum => Self::write_event(&mut spectrum_writer, &Event::Empty(tag.to_owned()), indent_level),
                        _ => ()
                    }
                }
                // leaving of tags: </tag>
                Ok(Event::End(ref tag)) => {
                    indent_level -= 1;
                    match tag.name() {
                        b"spectrum" if inside_spectrum => {
                            Self::write_event(&mut spectrum_writer, &Event::End(tag.to_owned()), indent_level);
                            let spectrum_xml = match std::str::from_utf8(spectrum_writer.into_inner().into_inner().as_slice()) {
                                Ok(string) => string.to_owned(),
                                Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_ms_two_spectra(): Error at std::str::from_utf8(): {}", err)
                            };
                            if Self::is_ms_two_spectrum(spectrum_xml.as_str()){
                                spectra.push(Spectrum::new(spectrum_xml.as_str(), indent_level));
                            }
                            spectrum_writer = Writer::new(Cursor::new(Vec::new()));
                            inside_spectrum = false;
                        },
                        _ if inside_spectrum => Self::write_event(&mut spectrum_writer, &Event::End(tag.to_owned()), indent_level),
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
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_content_before_spectrum_list(): error when reading MzMl-File, original error: {:?}", err)
        };
        let mut buf = Vec::new();
        let mut indent_level = 0;
        let mut writer = Writer::new(Cursor::new(Vec::new()));
        loop {
            match reader.read_event(&mut buf) {
                // entering tags: <tag (attr1="" attr2="" ...)>
                Ok(Event::Start(ref tag)) => {
                    match tag.name() {
                        b"spectrumList" => break,
                        _ => Self::write_event(&mut writer, &Event::Start(tag.to_owned()), indent_level)
                    };
                    indent_level += 1;
                },
                // empty tags: <tag (attr1="" attr2="" ...)/>
                Ok(Event::Empty(ref tag)) => Self::write_event(&mut writer, &Event::Empty(tag.to_owned()), indent_level),
                // declaration tags: <?xml ... ?>
                Ok(Event::Decl(ref tag)) => Self::write_event(&mut writer, &Event::Decl(tag.to_owned()), indent_level),
                // leaving of tags: </tag>
                Ok(Event::End(ref tag)) => {
                    indent_level -= 1;
                    Self::write_event(&mut writer, &Event::End(tag.to_owned()), indent_level);
                },
                Ok(Event::Eof) => break, // exits the loop when reaching end of file
                Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
                _ => (), // There are several other `Event`s we do not consider here
            }

            // if we don't keep a borrow elsewhere, we can clear the buffer to keep memory usage low
            buf.clear();
        }

        return match std::str::from_utf8(writer.into_inner().into_inner().as_slice()) {
            Ok(string) => string.to_owned(),
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_content_before_spectrum_list(): Error at std::str::from_utf8(): {}", err)
        }
    }

    pub fn get_chromatograms(&self) -> Box<Vec<Chromatrogram>> {
        let mut chromatograms: Vec<Chromatrogram> = Vec::new();
        let mut reader = match Reader::from_file(Path::new(self.file_path.as_str())) {
            Ok(reader) => reader,
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_chromatograms(): error when reading MzMl-File, original error: {:?}", err)
        };
        let mut buf = Vec::new();
        let mut inside_chromatogram: bool = false;
        let mut indent_level = 0;
        let mut chromatogram_writer = Writer::new(Cursor::new(Vec::new()));
        loop {
            match reader.read_event(&mut buf) {
                // entering tags: <tag (attr1="" attr2="" ...)>
                Ok(Event::Start(ref tag)) => {
                    match tag.name() {
                        b"chromatogram" => {
                            inside_chromatogram = true;
                            Self::write_event(&mut chromatogram_writer, &Event::Start(tag.to_owned()), indent_level);
                        },
                        _ if inside_chromatogram => Self::write_event(&mut chromatogram_writer, &Event::Start(tag.to_owned()), indent_level),
                        _ => ()
                    }
                    indent_level += 1;
                },
                // empty tags: <tag (attr1="" attr2="" ...)/>
                Ok(Event::Empty(ref tag)) => {
                    match tag.name() {
                        _ if inside_chromatogram => Self::write_event(&mut chromatogram_writer, &Event::Empty(tag.to_owned()), indent_level),
                        _ => ()
                    }
                }
                // leaving of tags: </tag>
                Ok(Event::End(ref tag)) => {
                    indent_level -= 1;
                    match tag.name() {
                        b"chromatogram" if inside_chromatogram => {
                            Self::write_event(&mut chromatogram_writer, &Event::End(tag.to_owned()), indent_level);
                            let chromatogram_xml = match std::str::from_utf8(chromatogram_writer.into_inner().into_inner().as_slice()) {
                                Ok(string) => string.to_owned(),
                                Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.get_chromatograms(): Error at std::str::from_utf8(): {}", err)
                            };
                            chromatograms.push(Chromatrogram::new(chromatogram_xml.as_str(), indent_level));
                            chromatogram_writer = Writer::new(Cursor::new(Vec::new()));
                            inside_chromatogram = false;
                        },
                        b"chromatogramList" => break,
                        _ if inside_chromatogram => Self::write_event(&mut chromatogram_writer, &Event::End(tag.to_owned()), indent_level),
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

    fn write_new_line(writer: &mut quick_xml::Writer<Cursor<Vec<u8>>>) {
        match writer.write(b"\n") {
            Ok(_) => (),
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.write_new_line(): {}", err)
        }
    }

    fn write_indention(writer: &mut quick_xml::Writer<Cursor<Vec<u8>>>, indention_level: usize) {
        match writer.write(mz_ml::indent(indention_level).as_bytes()) {
            Ok(_) => (),
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader.write_indention(): {}", err)
        }
    }

    fn write_event(writer: &mut quick_xml::Writer<Cursor<Vec<u8>>>, event: &Event, indention_level: usize) {
        Self::write_indention(writer, indention_level);
        match writer.write_event(event) {
            Ok(_) => (),
            Err(err) => panic!("proteomic::utility::mz_ml_reader::MzMlReader::write_event(): {}", err)
        };
        Self::write_new_line(writer);
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