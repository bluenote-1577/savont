use bincode;
use flexi_logger::style;
use clap::Parser;
use flexi_logger::{DeferredNow, Duplicate, FileSpec, Record};
use savont::asv_cluster;
use savont::cli;
use savont::constants::*;
use savont::kmer_comp;
use savont::classify;
use savont::seq_parse;
use savont::types;
use savont::alignment;
use savont::utils::*;
use savont::chimera;
use savont::taxonomy;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;
use sysinfo::System;
fn main() {
    let args = cli::Cli::parse();

    match &args.command {
        cli::Commands::Cluster(cluster_args) => {
            run_cluster(cluster_args, &args);
        }
        cli::Commands::Classify(classify_args) => {
            run_classify(classify_args, &args);
        }
    }
}

fn run_cluster(args: &cli::ClusterArgs, cli_args: &cli::Cli) {
    let output_dir = initialize_setup_cluster(args, cli_args);

    log::info!("Starting clustering...");

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

    let clusters = asv_cluster::cluster_reads_by_kmers(&twin_reads, &args, &output_dir);
    let clusters = asv_cluster::cluster_reads_by_snpmers(&twin_reads, &clusters, &args, &output_dir);
    let mut consensuses = alignment::align_and_consensus(&twin_reads, clusters, &args, &output_dir);

    // Generate pileups for quality estimation
    let mut pileups = alignment::generate_consensus_pileups(&twin_reads, &consensuses, &args);

    // Update consensus HP lengths from modal values calculated from pileups
    log::info!("Updating consensus HP lengths from pileup alignments");
    for (consensus, pileup) in consensuses.iter_mut().zip(pileups.iter()) {
        // Extract modal HP lengths from pileup
        let modal_hp_lengths: Vec<u8> = pileup.iter().map(|p| p.ref_hp_length).collect();
        consensus.hp_lengths = modal_hp_lengths;
    }

    // Estimate quality error rates from top 10% of clusters
    let quality_error_map = alignment::estimate_quality_error_rates(&pileups, &consensuses, 0.1);

    // Polish consensus sequences using Bayesian inference
    let mut low_qual_consensus = alignment::polish_consensuses(&mut pileups, &mut consensuses, &quality_error_map, &twin_reads, &args);
    log_memory_usage(true, "STAGE 1.9: Polished consensus sequences");

    // Decompress HPC sequences before merging and chimera detection
    log::info!("Decompressing {} HPC consensus sequences for merging and chimera detection", consensuses.len());
    for consensus in &mut consensuses {
        consensus.decompress();
    }

    // Decompress low quality consensus sequences as well
    log::info!("Decompressing {} low quality HPC consensus sequences", low_qual_consensus.len());
    for consensus in &mut low_qual_consensus {
        consensus.decompress();
    }

    alignment::write_consensus_fasta(&low_qual_consensus, &output_dir.join("low_quality_consensus_sequences.fasta"), "lowqual")
        .expect("Failed to write low_quality_consensus_sequences.fasta");

    // Merge similar consensus sequences based on alignment and depth (using decompressed sequences)
    // This also merges low quality consensuses into high quality ones
    let consensuses = alignment::merge_similar_consensuses(&twin_reads, consensuses, low_qual_consensus, &args);
    log_memory_usage(true, "STAGE 2: Merged similar consensus sequences");

    // Detect and filter chimeric consensus sequences
    if args.skip_chimera_detection {
        log::info!("Skipping chimera detection as per user request.");
        log_memory_usage(true, "STAGE 2.5: Skipped chimera detection");
        return;
    }
    let chimeras = chimera::detect_chimeras(&consensuses, &args);
    let consensuses = chimera::filter_chimeras(consensuses, &chimeras);
    log_memory_usage(true, "STAGE 2.5: Filtered chimeric consensus sequences");

    log::info!("Final consensus count after all filtering: {}", consensuses.len());

    // Write final consensus sequences after chimera filtering
    let output_dir = std::path::PathBuf::from(&args.output_dir);
    let final_fasta = output_dir.join(ASV_FILE);
    alignment::write_consensus_fasta(&consensuses, &final_fasta, "final")
        .expect(format!("Failed to write {}", ASV_FILE).as_str());
    log::info!("Wrote {} final consensus sequences to {}", consensuses.len(), ASV_FILE);

    // Write final cluster information
    let final_clusters = output_dir.join("final_clusters.tsv");
    alignment::write_clusters_tsv(&consensuses, &twin_reads, &final_clusters, "final")
        .expect("Failed to write final_clusters.tsv");
    log::info!("Wrote final cluster information to final_clusters.tsv");
}

fn run_classify(args: &cli::ClassifyArgs, cli_args: &cli::Cli) {
    // Initialize output directory and logger
    let _output_dir = initialize_setup_classify(args, cli_args);

    log::info!("Starting classification...");
    let db;
    if let Some(db_path) = &args.db_type.emu_db {
        log::info!("Using EMU database at {}", db_path);
        db = taxonomy::Database::load_emu(Path::new(db_path))
            .expect("Failed to load EMU database");
    }
    else if let Some(silva_path) = &args.db_type.silva_db {
        log::info!("Using Silva database at {}", silva_path);
        db = taxonomy::Database::load_silva(Path::new(silva_path))
            .expect("Failed to load Silva database");
    }
    else {
        log::error!("No database specified for classification. Use --emu-db or --silva-db.");
        std::process::exit(1);

    }

    classify::classify(&args, &db);

}

fn initialize_setup_classify(args: &cli::ClassifyArgs, cli_args: &cli::Cli) -> PathBuf {
    let output_dir = Path::new(&args.output_dir);

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

    // Initialize logger
    let log_spec = format!("{},skani=info", cli_args.log_level_filter().to_string());
    let filespec = FileSpec::default()
        .directory(output_dir)
        .basename("savont_classify");
    let _logger_handle = flexi_logger::Logger::try_with_str(log_spec)
        .expect("Something went wrong with logging")
        .log_to_file(filespec)
        .duplicate_to_stderr(Duplicate::Info)
        .format(my_own_format_colored)
        .format_for_files(my_own_format)
        .create_symlink("savont_classify_latest.log")
        .start()
        .expect("Something went wrong with creating log file");

    let command_args: Vec<String> = std::env::args().collect();
    log::info!("COMMAND: {}", command_args.join(" "));
    log::info!("VERSION: {}", env!("CARGO_PKG_VERSION"));

    // Initialize thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .stack_size(16 * 1024 * 1024)
        .build_global()
        .unwrap();

    output_dir.to_path_buf()
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

fn initialize_setup_cluster(args: &cli::ClusterArgs, cli_args: &cli::Cli) -> PathBuf {

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

    // Initialize logger with CLI-specified level
    let log_spec = format!("{},skani=info", cli_args.log_level_filter().to_string());
    let filespec = FileSpec::default()
        .directory(output_dir)
        .basename("savont");
    let _logger_handle = flexi_logger::Logger::try_with_str(log_spec)
        .expect("Something went wrong with logging")
        .log_to_file(filespec) // write logs to file
        .duplicate_to_stderr(Duplicate::Info) // print warnings and errors also to the console
        .format(my_own_format_colored) // use a simple colored format
        .format_for_files(my_own_format)
        .create_symlink("savont_latest.log")
        .start()
        .expect("Something went wrong with creating log file");

    let command_args: Vec<String> = std::env::args().collect();
    log::info!("COMMAND: {}", command_args.join(" "));
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

fn get_kmers_and_snpmers(args: &cli::ClusterArgs, output_dir: &PathBuf) -> types::KmerGlobalInfo {
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
    }
    return kmer_info;
}

fn get_twin_reads_from_kmer_info(
    kmer_info: &mut types::KmerGlobalInfo,
    args: &cli::ClusterArgs,
    _output_dir: &PathBuf,
    _cleaning_temp_dir: &PathBuf,
) -> Vec<types::TwinRead>{
    log::info!("Getting twin reads from snpmers...");
    let mut twin_reads_raw = kmer_comp::twin_reads_from_snpmers(kmer_info, &args);
    twin_reads_raw.sort_by(|a,b| b.est_id.unwrap_or(100.0).partial_cmp(&a.est_id.unwrap_or(100.0)).unwrap());
    return twin_reads_raw;
}