use std::io::prelude::*;
use std::fs::File;
use std::io::LineWriter;
use std::fs::OpenOptions;


pub struct SimpleLogger {
    log_file_path: String,
    log_file: LineWriter<File>
}

impl SimpleLogger {
    pub fn new(log_file_path: &str) -> SimpleLogger {
        let log_file = match OpenOptions::new().read(true).write(true).create(true).open(&log_file_path) {
            Ok(file) => file,
            Err(err) => panic!("SimpleLogger [{}] ERROR: {:?}", &log_file_path, err)
        };
        // create and return logger
        return SimpleLogger {
            log_file_path: log_file_path.to_owned(),
            log_file: LineWriter::new(log_file)
        }
    }

    pub fn write(&mut self, message: &String){
        match self.log_file.write_all(message.as_bytes()) {
            Ok(_byte_count) => (),
            Err(err) => println!("SimpleLogger [{}] ERROR: {:?}", self.log_file_path, err)
        }
    }

    pub fn flush(&mut self) {
        match self.log_file.flush() {
            Ok(_ok) => (),
            Err(err) => println!("SimpleLogger [{}] ERROR: {:?}", self.log_file_path, err)
        }
    }
}

impl Drop for SimpleLogger {
    fn drop(&mut self) {
        self.flush();
    }
}