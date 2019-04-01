# MaxDecoy

## Installation

1. Install rustup: `curl https://sh.rustup.rs -sSf | sh`
2. Install the Rust toolchain version 1.32.0 using rustup: `rustup install 1.32.0`
3. Go to the MaxDecoy-code-directory and build it with: `cargo build --release`

The binary file ends up in `target/release/max_decoy`.

## Database preparation
The folder `db` contains a SQL-schema for PostgreSQL which defines the Databasestructure. For production coment the testing partitions and uncomment the 100 partition for production.

## Usage
Copy the file `.env.example` to the folder where you start MaxDecoy, rename it `.env` and adjust it to your needs.   
Use `max_decoy --help` to show the subcommands of MaxDecoy and `max_decoy subcommand --help` to show the parameters of a single subcommand.
