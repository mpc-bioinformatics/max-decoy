use std::sync::Mutex;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::time::Duration;
use std::io::prelude::*;
use std::io::LineWriter;
use std::fs::OpenOptions;


const LOG_INTERVALL: i64 = 900;
const WAIT_DURATION: Duration = Duration::from_millis(1000);

pub struct AsyncPerformanceLogger {
    log_file_path: String,
    created_protein_counter: Arc<AtomicUsize>,
    created_peptide_counter: Arc<AtomicUsize>,
    created_peptide_protein_association_counter: Arc<AtomicUsize>,
    processed_protein_counter: Arc<AtomicUsize>,
    processed_peptide_counter: Arc<AtomicUsize>,
    processed_peptide_protein_association_counter: Arc<AtomicUsize>,
    stop_flag: Arc<Mutex<bool>>,
    thread_handle: Option<thread::JoinHandle<()>>   // 'option dance' see: https://users.rust-lang.org/t/spawn-threads-and-join-in-destructor/1613
}

impl AsyncPerformanceLogger {
    pub fn new(log_file_path: &str) -> AsyncPerformanceLogger {
        return AsyncPerformanceLogger {
            log_file_path: log_file_path.to_owned(),
            created_protein_counter: Arc::new(AtomicUsize::new(0)),
            created_peptide_counter: Arc::new(AtomicUsize::new(0)),
            created_peptide_protein_association_counter: Arc::new(AtomicUsize::new(0)),
            processed_protein_counter: Arc::new(AtomicUsize::new(0)),
            processed_peptide_counter: Arc::new(AtomicUsize::new(0)),
            processed_peptide_protein_association_counter: Arc::new(AtomicUsize::new(0)),
            stop_flag: Arc::new(Mutex::new(false)),
            thread_handle: None
        }
    }

    pub fn start_logging(&mut self) {
        // clone Arcs for moving into thread
        let created_protein_counter_ptr = self.created_protein_counter.clone();
        let created_peptide_counter_ptr = self.created_peptide_counter.clone();
        let created_peptide_protein_association_counter_ptr = self.created_peptide_protein_association_counter.clone();
        let processed_protein_counter_ptr = self.processed_protein_counter.clone();
        let processed_peptide_counter_ptr = self.processed_peptide_counter.clone();
        let processed_peptide_protein_association_counter_ptr = self.processed_peptide_protein_association_counter.clone();
        let stop_flag_ptr = self.stop_flag.clone();
        let log_file_path = self.log_file_path.clone();
        self.thread_handle = Some(
            thread::spawn(move || {
                //let mut loop_counter: usize = 0;
                let performance_file = match OpenOptions::new().read(true).write(true).create(true).open(log_file_path.as_str()) {
                    Ok(file) => file,
                    Err(err) => panic!("proteomic::utility::logger::performance_logger::PerformanceLogger.thread: error at opening file {}: {}", log_file_path, err)
                };
                let mut performance_file = LineWriter::new(performance_file);
                match performance_file.write(b"\"seconds\",\"created proteins\",\"processed proteins\",\"created peptides\",\"processed peptides\",\"created protein/peptide-association\",\"processed protein/peptide-association\"\n") {
                    Ok(_) => (),
                    Err(err) => println!("proteomic::utility::logger::performance_logger::PerformanceLogger.thread: Could not write to file: {}", err)
                }
                let start_at_sec = time::now().to_timespec().sec;
                let mut next_write_at = start_at_sec + LOG_INTERVALL;
                let mut stop_loop = false;
                // outer loop runs until stop_flag is set
                while !stop_loop {
                    let now_in_seconds = time::now().to_timespec().sec;
                    let mut write_intervall_expired = now_in_seconds >= next_write_at;
                    match stop_flag_ptr.lock() {
                        Ok(ref stop) if **stop => {
                            stop_loop = true;
                            write_intervall_expired = true;  // overwrite write_intervall_expired, so the thread will write a last time, regardless if the write intervall is expired or not
                        },
                        Ok(_) => (),
                        Err(_) => panic!("proteomic::utility::logger::performance_logger::PerformanceLogger.thread: tried to lock a poisoned mutex for 'stop_flag'")
                    }
                    if write_intervall_expired  {
                        next_write_at = now_in_seconds + LOG_INTERVALL;
                        match performance_file.write(
                            format!(
                                "{},{},{},{},{},{},{}\n",
                                now_in_seconds - start_at_sec,
                                created_protein_counter_ptr.load(Ordering::Relaxed),
                                processed_protein_counter_ptr.load(Ordering::Relaxed),
                                created_peptide_counter_ptr.load(Ordering::Relaxed),
                                processed_peptide_counter_ptr.load(Ordering::Relaxed),
                                created_peptide_protein_association_counter_ptr.load(Ordering::Relaxed),
                                processed_peptide_protein_association_counter_ptr.load(Ordering::Relaxed)
                            ).as_bytes()
                        ) {
                            Ok(_) => (),
                            Err(err) => println!("proteomic::utility::logger::performance_logger::PerformanceLogger.thread: Could not write to file: {}", err)
                        }
                    }
                    thread::sleep(WAIT_DURATION);
                }
            })
        )
    }

    pub fn increase_counter_by(&self, created_protein_counter_by: usize, created_peptide_counter_by: usize, created_peptide_protein_association_counter_by: usize, processed_protein_counter_by: usize, processed_peptide_counter_by: usize, processed_peptide_protein_association_counter_by: usize) {
        self.created_protein_counter.fetch_add(created_protein_counter_by, Ordering::Relaxed);
        self.created_peptide_counter.fetch_add(created_peptide_counter_by, Ordering::Relaxed);
        self.created_peptide_protein_association_counter.fetch_add(created_peptide_protein_association_counter_by, Ordering::Relaxed);
        self.processed_protein_counter.fetch_add(processed_protein_counter_by, Ordering::Relaxed);
        self.processed_peptide_counter.fetch_add(processed_peptide_counter_by, Ordering::Relaxed);
        self.processed_peptide_protein_association_counter.fetch_add(processed_peptide_protein_association_counter_by, Ordering::Relaxed);
    }
}

// 'Destructor' which waits for thread to stop
impl Drop for AsyncPerformanceLogger {
    fn drop(&mut self) {
        match self.stop_flag.lock() {
            Ok(mut stop) => *stop = true,
            Err(_) => println!("proteomic::utility::logger::performance_logger::PerformanceLogger.drop(): tried to lock a poisoned mutex for 'stop_flag'")
        }
        match self.thread_handle.take() {
            Some(handle) => match handle.join() {
                Ok(_) => (),
                Err(err) => panic!("proteomic::utility::logger::performance_logger::PerformanceLogger.drop(): error at handle.join(): {:?}", err)
            },
            None => ()
        }
    }
}