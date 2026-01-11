use clap::{Parser, Subcommand, ValueEnum, Args};
use crate::constants::CLI_HEADINGS;

#[derive(Parser, Debug)]
#[command(
    name = "savont",
    about = "savont - high-resolution ASV (Amplicon Sequence Variant) generation and taxonomic profiling for ONT R10.4/HiFi long-read amplicon sequencing",
    version,
    author,
    disable_help_subcommand = true,

)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,

    /// Logging verbosity level
    #[arg(short, long, value_enum, default_value = "debug", global = true)]
    pub log_level: LogLevel,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Cluster long reads of >~ 98% accuracy into ASVs (Amplicon Sequence Variants)
    #[command(name = "asv")]
    Cluster(ClusterArgs),

    /// Classify ASVs against a reference database and generate taxonomy abundance table at species/genus level
    #[command(name = "classify")]
    Classify(ClassifyArgs),

    /// Download reference databases for savont (EMU or SILVA)
    #[command(name = "download")]
    Download(DownloadArgs),
}

#[derive(Parser, Debug, Clone)]
pub struct ClusterArgs {
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
    #[arg(short, long, default_value = "11", help_heading = CLI_HEADINGS[0], hide = true)]
    pub c: usize,

    /// Minimum read length for reads 
    #[arg(long, default_value = "1100", help_heading = CLI_HEADINGS[0])]
    pub min_read_length: usize,

    /// Maximum read length for reads
    #[arg(long, default_value = "2000", help_heading = CLI_HEADINGS[0])]
    pub max_read_length: usize,

    /// Minimum estimated read accuracy (%) to include in clustering
    #[arg(long, default_value_t=98., help_heading = CLI_HEADINGS[0])]
    pub quality_value_cutoff: f64,

    /// Minimum base quality to be considered high-quality for SNPmer detection. Set lower for older reads (~18 for R9). 
    #[arg(long, default_value_t=25, help_heading = CLI_HEADINGS[0])]
    pub minimum_base_quality: u8,

    /// Use only forward strand k-mers (for strand-specific protocols)
    #[arg(short, long, help_heading = CLI_HEADINGS[1])]
    pub single_strand: bool,

    /// Minimum number of reads required to keep a cluster (ASV)
    #[arg(long, default_value_t=12, help_heading = CLI_HEADINGS[1])]
    pub min_cluster_size: usize,

    
    /// Bloom filter size in GB for k-mer filtering (0 = auto, increase for very large datasets)
    #[arg(short, long, default_value_t=0., help_heading = CLI_HEADINGS[1], hide=true)]
    pub bloom_filter_size: f64,

    /// Minimum depth required for sequences with ambiguous bases to be included in output
    #[arg(short, long, default_value_t=250, help_heading = CLI_HEADINGS[2])]
    pub n_depth_cutoff: usize,

    /// Use homopolymer-compressed sequences for clustering and consensus generation
    #[arg(short, long, default_value_t=false, help_heading = CLI_HEADINGS[2], hide=true)]
    pub use_hpc: bool,


    /// Mask low-quality bases in consensus sequences (set to 'N' if below posterior probability threshold)
    #[arg(long, help_heading = CLI_HEADINGS[2])]
    pub mask_low_quality: bool,

    /// Negative alternate posterior probability threshold (natural log scale) for base consensus. Higher = more stringent for low-quality consensuses. Do not set higher than min_depth * ln(error_rate). 
    #[arg(short, long, default_value_t=30., help_heading = CLI_HEADINGS[2])]
    pub posterior_threshold_ln: f64,

    
    /// Maximum number of reclustering iterations
    #[arg(long, default_value_t=10, help_heading = CLI_HEADINGS[1])]
    pub max_iterations_recluster: usize,

    /// Use more aggressive k-mer filtering (faster but may be non-deterministic)
    #[arg(long, help_heading = CLI_HEADINGS[1], hide = true)]
    pub aggressive_bloom: bool,

    /// Skip chimera detection step (not recommended)
    #[arg(long, hide=true)]
    pub skip_chimera_detection: bool,

    /// Disable SNPmer clustering (not recommended, uses only k-mers)
    #[arg(long, default_value_t=false, help_heading = CLI_HEADINGS[2], hide = true)]
    pub no_snpmers: bool,

    /// K-mer size for clustering (must be odd and < 24)
    #[arg(short, long, default_value = "17", help_heading = CLI_HEADINGS[1], hide = true)]
    pub kmer_size: usize,

    /// Blockmer suffix length for polymorphic marker detection (experimental)
    #[arg(long, default_value = "3", help_heading = CLI_HEADINGS[1], hide = true)]
    pub blockmer_length: usize,

    /// Use blockmers instead of SNPmers for polymorphic marker clustering (experimental)
    #[arg(long, default_value_t = false, help_heading = CLI_HEADINGS[1], hide = true)]
    pub use_blockmers: bool,

    /// Allowable errors for bi-chimeric detection (higher = more sensitive, slower)
    #[arg(long, default_value_t=1, help_heading = CLI_HEADINGS[3])]
    pub chimera_allowable_errors: usize,

    /// Length of near-perfect asv segment matches to consider for chimera detection (higher = less sensitive). Default is 1/10 of the minimum read length.
    #[arg(long, help_heading = CLI_HEADINGS[3])]
    pub chimera_detect_length: Option<usize>,


    /// Print help in markdown format
    #[arg(long, hide = true)]
    pub markdown_help: bool,

    // Legacy fields kept for compatibility with unused assembly code
    #[arg(skip)]
    pub hifi: bool,

    /// Try phasing heterogeneous clusters
    #[arg(long, help_heading = CLI_HEADINGS[2], hide=true)]
    pub phase_heterogeneous: bool,
}

#[derive(Parser, Debug, Clone)]
pub struct ClassifyArgs {
    /// Directory containing clustering results
    #[arg(short, long, required = true)]
    pub input_dir: String,

    /// Output directory for classification results. Default: same as the input directory
    #[arg(short, long)]
    pub output_dir: Option<String>,

    #[command(flatten)]
    pub db_type: DatabaseType,

    /// Number of threads to use for parallel processing
    #[arg(short, long, default_value = "20")]
    pub threads: usize,

    /// Minimum identity threshold for species-level classification (default: 99%)
    #[arg(long, default_value_t = 99.0)]
    pub species_threshold: f64,

    /// Minimum identity threshold for genus-level classification (default: 94.5%)
    #[arg(long, default_value_t = 94.5)]
    pub genus_threshold: f64,
}

#[derive(Args, Debug, Clone)]
#[group(required = true, multiple = false)]
pub struct DatabaseType {
    /// Emu database path
    #[arg(long, help_heading = "Database")]
    pub emu_db: Option<String>,

    /// SILVA database path
    #[arg(long, help_heading = "Database")]
    pub silva_db: Option<String>,
}

#[derive(Parser, Debug, Clone)]
pub struct DownloadArgs {
    /// Download location directory
    #[arg(short, long, required = true)]
    pub location: String,

    #[command(flatten)]
    pub db_type: DownloadDatabaseType,
}

#[derive(Args, Debug, Clone)]
#[group(required = true, multiple = true)]
pub struct DownloadDatabaseType {
    /// Download EMU database
    #[arg(long)]
    pub emu_db: bool,

    /// Download SILVA database
    #[arg(long)]
    pub silva_db: bool,
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
