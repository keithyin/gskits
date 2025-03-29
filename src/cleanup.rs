use std::{fs, path};

#[derive(Debug)]
pub struct AutoCleanFile {
    filename: String,
}

impl AutoCleanFile {
    pub fn new(fname: String) -> Self {
        Self { filename: fname }
    }
}

impl From<String> for AutoCleanFile {
    fn from(value: String) -> Self {
        Self { filename: value }
    }
}

impl Drop for AutoCleanFile {
    fn drop(&mut self) {
        if path::Path::new(&self.filename).exists() {
            let _ = fs::remove_file(&self.filename);
        }
    }
}
