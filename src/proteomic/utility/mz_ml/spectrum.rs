use std::collections::HashMap;
use std::fmt;
use std::io::prelude::*;
use std::fs::OpenOptions;
use std::path::{Path, PathBuf};

use quick_xml::Reader;
use quick_xml::events::Event;
use url::percent_encoding::{utf8_percent_encode, percent_decode, DEFAULT_ENCODE_SET};

use proteomic::utility::mz_ml;
use proteomic::models::mass;

const NAME_OF_MASS_TO_CHARGE_CV_PARAM: &str = "selected ion m/z";
const NAME_OF_CHARGE_CV_PARAM: &str = "charge state";

pub struct Spectrum {
    spectrum_id: String,
    scan_id: String,
    mass_to_charge_ratio: f64,
    charge: u8,
    xml: String,
    indent_level: usize
}

impl Spectrum {
    /// Creates new spectrum from spectrum xml.
    /// Panics if cvParam-tags with attribute name="charge state" or ="selected ion m/z" or spectrum-tag with arrribute name="scan" not found.
    ///
    /// # Arguments
    ///
    /// * `spectrum_xml` - A spectrum tag from mzML-file, make sure it is realy a spectrum tag because no validation is done.
    pub fn new(spectrum_xml: &str, indent_level: usize) -> Spectrum {
        let mut reader = Reader::from_str(spectrum_xml);
        let mut buf = Vec::new();
        let mut inside_selected_ion: bool = false;
        let mut charge: Option<u8> = None;
        let mut mass_to_charge_ratio: Option<f64> = None;
        let mut spectrum_id: Option<String> = None;
        let mut scan_id: Option<String> = None;
        loop {
            match reader.read_event(&mut buf) {
                // entering tags: <tag (attr1="" attr2="" ...)>
                Ok(Event::Start(ref tag)) => {
                    match tag.name() {
                        b"spectrum" => {
                            let attributes = mz_ml::parse_tag_attributes(tag);
                            if let Some(id_value) = attributes.get("id") {
                                scan_id = Some(Self::exctract_scan_id_from_spectrum_id(id_value.as_str()));
                                spectrum_id = Some(id_value.clone());
                            }
                        }
                        b"selectedIon" => inside_selected_ion = true,
                        _ => (),
                    }
                },
                // empty tags: <tag (attr1="" attr2="" ...)/>
                Ok(Event::Empty(ref tag)) => {
                    match tag.name() {
                        // process cvParam only if it is inside SelectedIon
                        b"cvParam" if inside_selected_ion => {
                            let attributes = mz_ml::parse_tag_attributes(tag);
                            if let Some(cv_param_type) = attributes.get("name") {
                                // check which name the cvParam has
                                match cv_param_type.as_str() {
                                    NAME_OF_MASS_TO_CHARGE_CV_PARAM => mass_to_charge_ratio = Some(Self::get_mass_to_charge_ratio_from_attributes(&attributes)),
                                    NAME_OF_CHARGE_CV_PARAM => charge = Some(Self::get_charge_from_attributes(&attributes)),
                                    _ => () // other cvParams than "selected ion m/z" and "charge state" are uninteresting
                                }
                            }
                        },
                        _ => ()
                    }
                }
                // leaving of tags: </tag>
                Ok(Event::End(ref tag)) => {
                    match tag.name() {
                        b"selectedIon" => inside_selected_ion = false,
                        _ => (),
                    }
                },
                Ok(Event::Eof) => {
                    if charge.is_some() & mass_to_charge_ratio.is_some() & scan_id.is_some() & spectrum_id.is_some() {
                        break;
                    } else {
                        panic!("proteomic::utility::mz_ml::spectrum::Spectrum::new(): Reached end of spectrum xml but either charge, mass-to-charge-ratio, title or id-ref was not found in:\n{}", spectrum_xml);
                    }
                }, // exits the loop when reaching end of file
                Err(e) => panic!("Error at position {}: {:?}", reader.buffer_position(), e),
                _ => (), // There are several other `Event`s we do not consider here
            }
            // if we don't keep a borrow elsewhere, we can clear the buffer to keep memory usage low
            buf.clear();
        }
        return Self {
            spectrum_id: spectrum_id.unwrap(),
            scan_id: scan_id.unwrap(),
            mass_to_charge_ratio: mass_to_charge_ratio.unwrap(),    // unwrap should be save at this point, because this function will panics if no mass to charge ratio is found
            charge: charge.unwrap(),                                // unwrap should be save at this point, because this function will panics if no charge is found
            xml: spectrum_xml.to_owned(),
            indent_level: indent_level
        }
    }

    fn get_charge_from_attributes(attributes: &HashMap<String, String>) -> u8 {
        match attributes.get("value") {
            Some(value) => match value.parse::<u8>() {
                Ok(charge) => return charge,
                Err(_) => panic!("proteomic::utility::mz_ml::spectrum::get_charge_from_attributes(): cvParam name is '{}' but could not parse charge unsigned integet (8 bit)", NAME_OF_CHARGE_CV_PARAM)
            },
            None => panic!("proteomic::utility::mz_ml::spectrum::get_charge_from_attributes(): cvParam name is '{}' there is no attribute 'value'", NAME_OF_CHARGE_CV_PARAM)
        }
    }

    fn get_mass_to_charge_ratio_from_attributes(attributes: &HashMap<String, String>) -> f64 {
        match attributes.get("value") {
            Some(value) => match value.parse::<f64>() {
                Ok(charge) => return charge,
                Err(_) => panic!("proteomic::utility::mz_ml::spectrum::get_charge_from_attributes(): cvParam name is '{}' but could not parse mass to charge ratio to float (64 bit)", NAME_OF_MASS_TO_CHARGE_CV_PARAM)
            },
            None => panic!("proteomic::utility::mz_ml::spectrum::get_charge_from_attributes(): cvParam name is '{}' there is no attribute 'value'", NAME_OF_MASS_TO_CHARGE_CV_PARAM)
        }
    }
    pub fn get_scan_id(&self) -> &str {
        return self.scan_id.as_str();
    }

    pub fn get_spectrum_id(&self) -> &str {
        return self.spectrum_id.as_str();
    }

    pub fn get_mass_to_charge_ratio(&self) -> &f64 {
        return &self.mass_to_charge_ratio;
    }

    pub fn get_charge(&self) -> &u8 {
        return &self.charge;
    }

    pub fn get_xml(&self) -> &str {
        return self.xml.as_str();
    }

    pub fn get_indent_level(&self) -> usize {
        return self.indent_level;
    }

    pub fn to_string(&self) -> String {
        return format!(
            "proteomic::utility::mz_ml::spectrum::Spectrum\n\ttitle => {}\n\tm/z => {}\n\tcharge => {}\n\txml => {}",
            if !self.scan_id.is_empty() { self.scan_id.as_str() } else { self.spectrum_id.as_str() },
            self.mass_to_charge_ratio,
            self.charge,
            self.xml
        );
    }

    /// Returns the title but will replace all character which are no diget, letter or a hyphen with an underscore.
    pub fn get_filename(&self) -> String {
        let filename: &str = if !self.scan_id.is_empty() { self.scan_id.as_str() } else { self.spectrum_id.as_str() };
        let iter = utf8_percent_encode(filename, DEFAULT_ENCODE_SET);
        return iter.collect::<String>();
    }

    /// Creates mzML-file with this spectrum only.
    ///
    /// # Arguments
    ///
    /// * `content_before_spectrum_list` - Content from original mzML-file, before spectrumListTag.
    pub fn to_mz_ml(&self, content_before_spectrum_list: &str, destination_folder: &Path, file_suffix: &str) {
        // create index list for mzML with open indexList-tag
        let mut index_list = format!("{}<indexList count=\"1\">\n",mz_ml::indent(1));
        // open index-tag for spectrums
        index_list.push_str(format!("{}<index name=\"spectrum\">\n", mz_ml::indent(2)).as_str());
        // begin with xml
        let mut mz_ml_content: String = content_before_spectrum_list.to_owned();
        // open spectrumList-tag
        mz_ml_content.push_str(format!("{}<spectrumList count=\"1\" defaultDataProcessingRef=\"pwiz_Reader_conversion\">\n",mz_ml::indent(self.indent_level - 1)).as_str());
        // get offset to the end of content
        let mut offset = mz_ml_content.len();
        // and add number of indention characters of this spectrum to get spectrum's offset
        offset += Self::count_chars_before_tag_starts(self.xml.as_str());
        // add offset-tag for this spectrum
        index_list.push_str(format!("{}<offset idRef=\"{}\">{}</offset>\n", mz_ml::indent(3), self.spectrum_id, offset).as_str());
        // add spectrum's xml to mzML
        mz_ml_content.push_str(format!("{}\n", self.xml.as_str()).as_str());
        // close spectrumList-tag
        mz_ml_content.push_str(format!("{}</spectrumList>\n", mz_ml::indent(self.indent_level - 1)).as_str());
        // close index-tag for spectrums
        index_list.push_str(format!("{}</index>\n", mz_ml::indent(2)).as_str());
        // close run-tag (opened in content_before_spectrum_list)
        mz_ml_content.push_str(format!("{}</run>\n", mz_ml::indent(2)).as_str());
        // close mzML-tag (opened in content_before_spectrum_list)
        mz_ml_content.push_str(format!("{}</mzML>\n", mz_ml::indent(1)).as_str());
        // close indexList-tag
        index_list.push_str(format!("{}</indexList>\n", mz_ml::indent(1)).as_str());
        // get offset to the end of content
        offset = mz_ml_content.len();
        // and add number of indention characters of index_list (used a bit later)
        offset += Self::count_chars_before_tag_starts(index_list.as_str());
        // add index_list
        mz_ml_content.push_str(index_list.as_str());
        // add indexList-offset
        mz_ml_content.push_str(format!("{}<indexListOffset>{}</indexListOffset>\n", mz_ml::indent(1), offset).as_str());
        // open fileChecksum-tag (as in mzML-specification defined, the opening tag must included in sha1-hash)
        mz_ml_content.push_str(format!("{}<fileChecksum>", mz_ml::indent(1)).as_str());
        let mut sha1_hash = sha1::Sha1::new();
        sha1_hash.update(mz_ml_content.as_bytes());
        // add sha1_hash and close fileChecksum-tag
        mz_ml_content.push_str(format!("{}</fileChecksum>\n", sha1_hash.digest().to_string()).as_str());
        // close indexedmzML-tag
        mz_ml_content.push_str("</indexedmzML>");
        // create mzML-file-path
        let mut mz_ml_file_path = destination_folder.to_path_buf();
        let mut filename = self.get_filename();
        if file_suffix.len() > 0 {
            filename.push_str("_");
            filename.push_str(file_suffix);
        }
        mz_ml_file_path.push(filename.as_str());
        mz_ml_file_path.set_extension("mzML");
        // open and write file
        let mut file = match OpenOptions::new().read(true).write(true).create(true).open(mz_ml_file_path) {
            Ok(file) => file,
            Err(err) => panic!("proteomic::utility::mz_ml::spectrum::Spectrum.to_mz_ml(): Coul not opening mzML-file: {}", err)
        };
        match file.write_all(mz_ml_content.as_bytes()) {
            Ok(_) => (),
            Err(err) => panic!("proteomic::utility::mz_ml::spectrum::Spectrum.to_mz_ml(): Could not write mzML-file: {}", err)
        }
    }

    /// Count chars before tag begins. Designed to count the indention characters to add them to offset.
    fn count_chars_before_tag_starts(tag: &str) -> usize {
        for (idx, letter) in tag.chars().enumerate() {
            if letter == '<' {
                return idx;
            }
        }
        panic!("proteomic::utility::mz_ml::spectrum::Spectrum::count_chars_before_tag_starts(): tag does not start at all")
    }

    fn exctract_scan_id_from_spectrum_id(spectrum_id: &str) -> String {
        let pos = onig::Regex::new(r"scan=.*\b").unwrap().find(spectrum_id);
        match pos {
            Some((beg, end)) => {
                return match String::from(&spectrum_id[beg..end]).split('=').into_iter().last() {
                    Some(id) => id.to_owned(),
                    None => "".to_owned()
                };
            },
            None =>
                return String::new()
        }
    }
}

impl fmt::Display for Spectrum {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.to_string())
    }
}