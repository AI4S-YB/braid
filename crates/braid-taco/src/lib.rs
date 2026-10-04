//! TACO meta-assembly.
//!
//! Independent Rust implementation of tacorna/taco 0.7.3
//! (commit `eeaeb879b8622365123edbc61ebc100d84194b80`).
//! Upstream is MIT, copyright Matthew Iyer, Yashar Niknafs, and Balaji Pandian.
//! This crate is not a copy of those sources.

mod assemble;
mod changepoint;
mod graph;
mod path;
mod splice;
mod types;
mod util;

pub use assemble::{assemble, write_taco_texts, TacoFailure, TacoSettings, TacoTexts};
