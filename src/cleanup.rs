use std::{fs, ops::Deref, path};


/// if the object is dropped, the file will be deleted
#[derive(Debug)]
pub struct AutoCleanFile {
    filename: String,
}

impl AutoCleanFile {
    pub fn new(fname: String) -> Self {
        Self { filename: fname }
    }
}

impl Deref for AutoCleanFile {
    type Target = String;

    fn deref(&self) -> &Self::Target {
        &self.filename
    }
}

impl From<String> for AutoCleanFile {
    fn from(value: String) -> Self {
        Self { filename: value }
    }
}

impl From<&String> for AutoCleanFile {
    fn from(value: &String) -> Self {
        Self {
            filename: value.to_string(),
        }
    }
}

impl From<&str> for AutoCleanFile {
    fn from(value: &str) -> Self {
        Self {
            filename: value.to_string(),
        }
    }
}

impl Drop for AutoCleanFile {
    fn drop(&mut self) {
        if path::Path::new(&self.filename).exists() {
            let _ = fs::remove_file(&self.filename);
        }
    }
}
