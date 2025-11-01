use clap::{Parser, ValueEnum};
use crate::constants::CLI_HEADINGS;

#[derive(Parser, Debug)]
#[command(
    name = "savont",
    about =
    "savont - high-resolution 16S rRNA clustering and denoising for long-read metagenomic data\n\nEXAMPLE: savont reads.fq.gz -o output_directory -t 50",
    version,
    author
)]

#[derive(Default, Clone)]
pub struct Cli {
    /// Input read file(s) in FASTQ or FASTA format (.gz supported). Multiple files are concatenated.
    #[arg(num_args = 1.., required = true, value_name = "FASTQ/FASTA")]
    pub input_files: Vec<String>,

    /// Output directory for results (created if it does not exist)
    #[arg(short, long, default_value = "savont-out")]
    pub output_dir: String,

    /// Number of threads to use for parallel processing
    #[arg(short, long, default_value = "20")]
    pub threads: usize,

    /// Delete intermediate files after completion to save disk space
    #[arg(long, hide = true)]
    pub clean_dir: bool,

    /// K-mer sampling rate: select 1 out of every C k-mers (higher = faster, less memory, slightly less sensitive)
    #[arg(short, long, default_value = "11", help_heading = CLI_HEADINGS[1], hide = true)]
    pub c: usize,

    /// Use only forward strand k-mers (for strand-specific protocols)
    #[arg(short, long, help_heading = CLI_HEADINGS[1])]
    pub single_strand: bool,

    /// Minimum number of reads required to keep a cluster (ASV)
    #[arg(long, default_value_t=12, help_heading = CLI_HEADINGS[1])]
    pub min_cluster_size: usize,

    /// Minimum estimated read accuracy (%) to include in clustering
    #[arg(long, default_value_t=99., help_heading = CLI_HEADINGS[1])]
    pub quality_value_cutoff: f64,

    /// Bloom filter size in GB for k-mer filtering (0 = auto, increase for very large datasets)
    #[arg(short, long, default_value_t=0., help_heading = CLI_HEADINGS[1])]
    pub bloom_filter_size: f64,

    /// Use more aggressive k-mer filtering (faster but may be non-deterministic)
    #[arg(long, help_heading = CLI_HEADINGS[1], hide = true)]
    pub aggressive_bloom: bool,

    /// Logging verbosity level
    #[arg(short, long, value_enum, default_value = "debug")]
    pub log_level: LogLevel,

    /// Skip chimera detection step (not recommended)
    #[arg(long)]
    pub skip_chimera_detection: bool,

    /// Disable SNPmer clustering (not recommended, uses only k-mers)
    #[arg(long, default_value_t=false, help_heading = CLI_HEADINGS[2], hide = true)]
    pub no_snpmers: bool,

    /// K-mer size for clustering (must be odd and < 24)
    #[arg(short, long, default_value = "17", help_heading = CLI_HEADINGS[1], hide = true)]
    pub kmer_size: usize,

    /// Print help in markdown format
    #[arg(long, hide = true)]
    pub markdown_help: bool,

    // Legacy fields kept for compatibility with unused assembly code
    #[arg(skip)]
    pub hifi: bool,

    /// Do not assume 16S full-length reads
    #[arg(long, help_heading = CLI_HEADINGS[2])]
    pub not_full_16s: bool,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

impl Default for LogLevel {
    fn default() -> Self {
        LogLevel::Debug
    }
}


impl Cli {
    pub fn log_level_filter(&self) -> log::LevelFilter {
        match self.log_level {
            LogLevel::Error => log::LevelFilter::Error,
            LogLevel::Warn => log::LevelFilter::Warn,
            LogLevel::Info => log::LevelFilter::Info,
            LogLevel::Debug => log::LevelFilter::Debug,
            LogLevel::Trace => log::LevelFilter::Trace,
        }
    }

    pub fn to_string(&self) -> String {
        format!("{:?}", self)
    }
}
