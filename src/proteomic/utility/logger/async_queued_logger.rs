use std::sync::Mutex;
use std::sync::Arc;
use std::thread;
use std::time::Duration;
use std::collections::VecDeque;
use std::io::prelude::*;
use std::io::LineWriter;
use std::fs::OpenOptions;


const WAIT_DURATION: Duration = Duration::from_millis(100);

pub struct AsyncQueuedLogger {
    log_file_path: String,
    queue: Arc<Mutex<VecDeque<String>>>,
    stop_flag: Arc<Mutex<bool>>,
    thread_handle: Option<thread::JoinHandle<()>>   // 'option dance' see: https://users.rust-lang.org/t/spawn-threads-and-join-in-destructor/1613
}

impl AsyncQueuedLogger {
    pub fn new(log_file_path: &str) -> AsyncQueuedLogger {
        // create atomic pointer for logger with clones for thread (prefixed with 't')
        let stop_flag = Arc::new(Mutex::new(false));
        let t_stop_flag = stop_flag.clone();
        let queue = Arc::new(Mutex::new(VecDeque::new()));
        let t_queue = queue.clone();
        let t_log_file_path = log_file_path.to_owned();
        // create and return logger
        return AsyncQueuedLogger {
            log_file_path: log_file_path.to_owned(),
            queue: queue,
            stop_flag: stop_flag,
            thread_handle: Some(
                thread::spawn(move || {
                    //let mut loop_counter: usize = 0;
                    let log_file = match OpenOptions::new().read(true).write(true).create(true).open(&t_log_file_path) {
                        Ok(file) => file,
                        Err(err) => panic!("AsyncQueuedLogger [{}] ERROR: {:?}", &t_log_file_path, err)
                    };
                    let mut log_file = LineWriter::new(log_file);
                    // outer loop runs until stop_flag is set
                    'outer: loop {
                        // inner loop runs until queue is empty respectively until pop_front returns None
                        'inner: loop {
                            if let Ok(mut queue) = t_queue.lock() {
                                match queue.pop_front() {
                                    Some(message) =>{
                                        println!("{}", &message);
                                        match log_file.write_all(message.as_bytes()) {
                                            Ok(_byte_count) => (),
                                            Err(err) => println!("AsyncQueuedLogger [{}] ERROR: {:?}", &t_log_file_path, err)
                                        }
                                    }
                                    None => break 'inner
                                }
                            }
                            match log_file.flush() {
                                Ok(_ok) => (),
                                Err(err) => println!("AsyncQueuedLogger [{}] ERROR: {:?}", &t_log_file_path, err)
                            }
                        }
                        // break outer loop if stop_flag is true else threads sleeps for a short duration and continues with outer loop
                        if let Ok(stop) = t_stop_flag.lock() {
                            if *stop {
                                break 'outer;
                            }
                        } else {
                            panic!("AsyncQueuedLogger::thread tried to lock a poisoned mutex for 'stop_flag'");
                        }
                        thread::sleep(WAIT_DURATION);
                    }
                })
            )
        }
    }

    pub fn push_back(&self, message: String) -> usize {
        if let Ok(mut queue) = self.queue.lock() {
            queue.push_back(message);
            queue.len()
        } else {
            panic!("AsyncQueuedLogger::add_work() tried to lock a poisoned mutex for 'queue'");
        }
    }
}

// 'Destructor' which waits for thread to stop
impl Drop for AsyncQueuedLogger {
    fn drop(&mut self) {
        match self.stop_flag.lock() {
            Ok(mut stop) => *stop = true,
            Err(err) => println!("{:?}", err)
        }
        match self.thread_handle.take() {
            Some(handle) => match handle.join() {
                Ok(_ok) => (),
                Err(err) => panic!("{:?}", err)
            },
            None => ()
        }
        let mut temp = String::new();
        for message in self.queue.lock().unwrap().iter() {
            temp.push_str(message);
        }
        if temp.len() > 0 {
            println!("AsyncQueuedLogger [{}]: Remaining messages\n{}", self.log_file_path, temp);
        }
    }
}