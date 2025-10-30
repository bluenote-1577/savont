use clap::{Parser, ValueEnum};
use crate::constants::CLI_HEADINGS;

use crate::constants::{IDENTITY_THRESHOLDS, ID_THRESHOLD_ITERS};

#[derive(Parser, Debug)]
#[command(
    name = "savont",
    about = 
    "savont - high-resolution metagenomic assembly with noisy long reads. See online documentation for full options. \n\nEXAMPLE (Nanopore R10): savont nanopore_reads.fq.gz -o output_directory -t 50\nEXAMPLE (PacBio HiFi): savont pacbio_reads.fq.gz -o output_directory -t 50 --hifi",
    version,
    author
)]

#[derive(Default, Clone)]
pub struct Cli {
    /// Input read file(s) -- multiple files are concatenated
    #[arg(num_args = 1.., required = true, value_name = "FASTQ/FASTA (.gz)")]
    pub input_files: Vec<String>,

    /// (DEFAULT) R10 nanopore mode for sup/hac data (> ~97% median accuracy). Specifying this flag does not do anything for now.
    #[arg(long, help_heading = CLI_HEADINGS[0])]
    pub nano_r10: bool,

    /// PacBio HiFi mode -- assumes less chimericism and higher accuracy
    #[arg(long, help_heading = CLI_HEADINGS[0])]
    pub hifi: bool,

    /// 16s rRNA mode -- optimizes parameters for full-length 16s assembly
    #[arg(long, help_heading = CLI_HEADINGS[0], default_value_t=true)]
    pub full_length_16s: bool,
        
    /// Output directory for results; created if it does not exist
    #[arg(short, long, default_value = "savont-out")]
    pub output_dir: String,

    /// Number of threads to use for processing
    #[arg(short, long, default_value = "20")]
    pub threads: usize,

     /// Do not dump large intermediate data to disk (intermediate data is useful for rerunning)
    #[arg(long)]
    pub clean_dir: bool,
   
    /// Compression ratio (1/c k-mers selected). Must be <= 15  
    #[arg(short, long, default_value = "11", help_heading = CLI_HEADINGS[1])]
    pub c: usize,

    /// Compression ratio (1/c k-mers selected). Must be <= 15  
    #[arg(short, long,  help_heading = CLI_HEADINGS[1])]
    pub single_strand: bool,

    /// Minimum cluster size for k-mer clustering
    #[arg(long, default_value_t=12, help_heading = CLI_HEADINGS[1])]
    pub min_cluster_size: usize,
        
    /// Disallow reads with < % identity for graph building (estimated from base qualities) 
    #[arg(long, default_value_t=98., help_heading = CLI_HEADINGS[1])]
    pub quality_value_cutoff: f64,

    /// Minimum overlap length for graph construction
    #[arg(long, default_value_t=500, help_heading = CLI_HEADINGS[1])]
    pub min_ol: usize,
        
    /// Bloom filter size in GB. Increase for massive datasets
    #[arg(short, long, default_value_t=0., help_heading = CLI_HEADINGS[1])]
    pub bloom_filter_size: f64,

    /// More aggressive filtering of low-abundance k-mers. May be non-deterministic
    #[arg(long, help_heading = CLI_HEADINGS[1])]
    pub aggressive_bloom: bool,
    
    /// Verbosity level. Warning: trace is very verbose
    #[arg(short, long, value_enum, default_value = "debug")]
    pub log_level: LogLevel,

    /// Output contigs with >= this number of reads
    #[arg(long, default_value_t = 1, help_heading = "Output thresholds")]
    pub min_reads_contig: usize,

    /// Remove singleton contigs with <= this estimated coverage depth (DP1 coverage; 99% identity coverage)
    #[arg(long, default_value_t=3., help_heading = "Output thresholds")]
    pub singleton_coverage_threshold: f64,

    /// Remove contigs with <= this estimated coverage depth and <= 2 reads (DP1 coverage; 99% identity coverage)
    #[arg(long, default_value_t=1., help_heading = "Output thresholds")]
    pub secondary_coverage_threshold: f64,

    /// Remove all contigs with <= this estimated coverage depth (DP1 coverage; 99% identity coverage)
    #[arg(long, default_value=None, help_heading = "Output thresholds")]
    pub absolute_coverage_threshold: Option<f64>,

    /// No polishing (not recommended)
    #[arg(long, default_value_t=false, help_heading = CLI_HEADINGS[2], hide = true)]
    pub no_polish: bool,

    /// Disable usage of SNPmers (not recommended)
    #[arg(long, default_value_t=false, help_heading = CLI_HEADINGS[2], hide = true)]
    pub no_snpmers: bool,

    /// K-mer size (must be odd and < 24)
    #[arg(short, long, default_value = "17", help_heading = CLI_HEADINGS[1], hide = true)]
    pub kmer_size: usize,

    /// Soft clips with < this # of bases are allowed for alignment
    #[arg(long, default_value_t=300, help_heading = CLI_HEADINGS[3], hide = true)]
    pub maximal_end_fuzz: usize, 

    /// Maximum bubble length to pop; keep alternates
    #[arg(long, default_value_t=500000, help_heading = CLI_HEADINGS[4], hide = true)]
    pub max_bubble_threshold: usize,

    /// Print this markdown document
    #[arg(long, hide = true)]
    pub markdown_help: bool,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum Preset{
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
