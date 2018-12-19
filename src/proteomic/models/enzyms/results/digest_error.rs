pub struct DigestError {
    message: String
}

impl DigestError {
    pub fn new(message: &str) -> DigestError {
        return DigestError {
            message: message.to_owned()
        }
    }

    pub fn get_message(&self) -> &String {
        return &self.message;
    }

}