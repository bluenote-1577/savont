pub mod seeding;
pub mod types;
pub mod kmer_comp;
//pub mod mapping;
pub mod cli;
pub mod constants;
//pub mod polishing_mod;
pub mod seq_parse;
pub mod map_processing;
pub mod utils;
pub mod asv_cluster;
pub mod alignment;

//pub mod cbloom;
//
//#[cfg(target_arch = "x86_64")]
//pub mod avx2_seeding;
//#[cfg(target_arch = "x86_64")]
//pub mod avx2_chaining;

// Use of a mod or pub mod is not actually necessary.
pub mod built_info {
   // The file has been placed there by the build script.
   // include!(concat!(env!("OUT_DIR"), "/built.rs"));
}
