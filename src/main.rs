use bincode;
use flexi_logger::style;
use clap::Parser;
use flexi_logger::{DeferredNow, Duplicate, FileSpec, Record};
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use savont::asv_cluster;
use savont::cli;
use savont::constants::*;
use savont::kmer_comp;
use savont::map_processing;
use savont::seq_parse;
use savont::types;
use savont::types::HeavyCutOptions;
use savont::types::OverlapAdjMap;
use savont::utils::*;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;
use sysinfo::System;
fn main() {
    let total_start_time = Instant::now();
    let mut args = cli::Cli::parse();

    let output_dir = initialize_setup(&mut args);

    log::info!("Starting assembly...");

    // Step 1: Process k-mers, count k-mers, and get SNPmers
    let mut kmer_info = get_kmers_and_snpmers(&args, &output_dir);
    log_memory_usage(true, "STAGE 1: Obtained SNPmers");

    // Step 1.5: Get twin reads from SNPmers
    let twin_reads = get_twin_reads_from_kmer_info(
        &mut kmer_info,
        &args,
        &output_dir,
        &output_dir.join("binary_temp"),
    );
    log_memory_usage(true, "STAGE 1.5: Obtained twin reads from SNPmers");

    let cluster = asv_cluster::cluster_reads_by_kmers(&twin_reads, &args, &output_dir);
    asv_cluster::cluster_reads_by_snpmers(&twin_reads, &cluster, &args, &output_dir);
}

fn my_own_format_colored(
    w: &mut dyn std::io::Write,
    now: &mut DeferredNow,
    record: &Record,
) -> Result<(), std::io::Error> {
    let mut paintlevel = record.level();
    if paintlevel == log::Level::Info {
        paintlevel = log::Level::Debug;
    }
    write!(
        w,
        "({}) {} [{}] {}",
        now.format(TS_DASHES_BLANK_COLONS_DOT_BLANK),
        style(paintlevel).paint(record.level().to_string()),
        record.module_path().unwrap_or(""),
        &record.args()
    )
}

fn my_own_format(
    w: &mut dyn std::io::Write,
    now: &mut DeferredNow,
    record: &Record,
) -> Result<(), std::io::Error> {
    write!(
        w,
        "({}) {} [{}] {}",
        now.format(TS_DASHES_BLANK_COLONS_DOT_BLANK),
        record.level(),
        record.module_path().unwrap_or(""),
        &record.args()
    )
}

fn initialize_setup(args: &mut cli::Cli) -> PathBuf {

    if args.markdown_help {
        let markdown_options = clap_markdown::MarkdownOptions::default();
        markdown_options.show_table_of_contents(true);
        clap_markdown::print_help_markdown::<cli::Cli>();
        std::process::exit(0);
    }


    for file in &args.input_files {
        if !Path::new(file).exists() && file != MAGIC_EXIST_STRING{
            eprintln!(
                "ERROR [savont] Input file {} does not exist. Exiting.",
                file
            );
            std::process::exit(1);
        }
    }

    let output_dir = Path::new(args.output_dir.as_str());

    if !output_dir.exists() {
        std::fs::create_dir_all(output_dir).expect("Could not create output directory. Exiting.");
    } else {
        if !output_dir.is_dir() {
            eprintln!(
                "ERROR [savont] Output directory specified by `-o` exists and is not a directory."
            );
            std::process::exit(1);
        }
    }

    let binary_temp_dir = output_dir.join("binary_temp");
    if !binary_temp_dir.exists(){
        std::fs::create_dir_all(&binary_temp_dir).expect("Could not create temp directory for binary files");
    } else {
        if !binary_temp_dir.is_dir() {
            panic!("Could not create temp directory for binary files. Exiting.");
        }
    }

    // Initialize logger with CLI-specified level
    let log_spec = format!("{},skani=info", args.log_level_filter().to_string());
    let filespec = FileSpec::default()
        .directory(output_dir)
        .basename("savont");
    flexi_logger::Logger::try_with_str(log_spec)
        .expect("Something went wrong with logging")
        .log_to_file(filespec) // write logs to file
        .duplicate_to_stderr(Duplicate::Info) // print warnings and errors also to the console
        .format(my_own_format_colored) // use a simple colored format
        .format_for_files(my_own_format)
        .start()
        .expect("Something went wrong with creating log file");

    let cli_args: Vec<String> = std::env::args().collect();
    log::info!("COMMAND: {}", cli_args.join(" "));
    log::info!("VERSION: {}", env!("CARGO_PKG_VERSION"));
    log::info!("SYSTEM NAME: {}", System::name().unwrap_or(format!("Unknown")));
    log::info!("SYSTEM HOST NAME: {}", System::host_name().unwrap_or(format!("Unknown")));
    //log::debug!("BINARY BUILD DATE: {}",  built_info::BUILT_TIME_UTC);
        // The built info is available in the `built` module

    // Validate k-mer size
    if args.kmer_size % 2 == 0 {
        log::error!("K-mer size must be odd");
        std::process::exit(1);
    }
    // Initialize thread pool, bigger stack size because sorting k-mers fails otherwise...
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .stack_size(16 * 1024 * 1024)
        .build_global()
        .unwrap();

    return output_dir.to_path_buf();
}

fn get_kmers_and_snpmers(args: &cli::Cli, output_dir: &PathBuf) -> types::KmerGlobalInfo {
    let saved_input = args.input_files == [MAGIC_EXIST_STRING];

    let binary_temp_dir = output_dir.join("binary_temp");
    let snpmer_info_path = binary_temp_dir.join("snpmer_info.bin");

    let kmer_info;
    if saved_input {
        if !snpmer_info_path.exists() {
            log::error!("No input files provided. See --help for usage.");
            std::process::exit(1);
        }
    }

    if saved_input && snpmer_info_path.exists() {
        kmer_info = bincode::deserialize_from(BufReader::new(
            File::open(snpmer_info_path).unwrap(),
        ))
        .unwrap();
        log::info!("Loaded snpmer info from file.");
    } else {
        let start = Instant::now();
        let big_kmer_map = seq_parse::read_to_split_kmers(args.kmer_size, args.threads, &args);
        log::info!(
            "Time elapsed in for counting k-mers is: {:?}",
            start.elapsed()
        );

        let start = Instant::now();
        //kmer_info = kmer_comp::get_snpmers(big_kmer_map, args.kmer_size, &args);
        kmer_info = kmer_comp::get_snpmers_inplace_sort(big_kmer_map, args.kmer_size, &args);
        log::info!(
            "Time elapsed in for parsing snpmers is: {:?}",
            start.elapsed()
        );

        if !args.clean_dir{
            bincode::serialize_into(
                BufWriter::new(File::create(snpmer_info_path).unwrap()),
                &kmer_info,
            )
            .unwrap();
        }
    }
    return kmer_info;
}

fn get_twin_reads_from_kmer_info(
    kmer_info: &mut types::KmerGlobalInfo,
    args: &cli::Cli,
    output_dir: &PathBuf,
    cleaning_temp_dir: &PathBuf,
) -> Vec<types::TwinRead>{
    log::info!("Getting twin reads from snpmers...");
    let mut twin_reads_raw = kmer_comp::twin_reads_from_snpmers(kmer_info, &args);
    twin_reads_raw.sort_by(|a,b| b.est_id.unwrap_or(100.0).partial_cmp(&a.est_id.unwrap_or(100.0)).unwrap());
    log_memory_usage(true, "STAGE 1.5: Initially obtained dirty twin reads");
    return twin_reads_raw;
}