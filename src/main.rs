use flexi_logger::style;
use clap::Parser;
use flexi_logger::{DeferredNow, Duplicate, FileSpec, Record};
use savont::asv_cluster;
use savont::cli;
use savont::constants::*;
use savont::kmer_comp;
use savont::classify;
use savont::seeding;
use savont::seq_parse;
use savont::types;
use savont::alignment;
use savont::types::ConsensusSequence;
use savont::types::decode_kmer48;
use savont::utils::*;
use savont::chimera;
use savont::taxonomy;
use savont::download;
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;
use sysinfo::System;
use fxhash::FxHashSet;
fn main() {
    let args = cli::Cli::parse();

    match &args.command {
        cli::Commands::Cluster(cluster_args) => {
            run_cluster(cluster_args, &args);
        }
        cli::Commands::Classify(classify_args) => {
            run_classify(classify_args, &args);
        }
        cli::Commands::Download(download_args) => {
            run_download(download_args, &args);
        }
    }
}

fn run_cluster(args: &cli::ClusterArgs, cli_args: &cli::Cli) {
    let output_dir = initialize_setup_cluster(args, cli_args);

    // Create temp directory for intermediate files
    let temp_dir = output_dir.join("temp");
    std::fs::create_dir_all(&temp_dir).expect("Failed to create temp directory");
    log::info!("Created temp directory for intermediate files: {}", temp_dir.display());

    log::info!("Starting clustering...");

    // Step 1: Process k-mers, count k-mers, and get SNPmers
    let (mut kmer_info, mut blockmer_info) = get_kmers_and_snpmers(&args, &temp_dir);
    log_memory_usage(true, "STAGE 1: Obtained SNPmers");
    log::info!("Using blockmers: {}", args.use_blockmers);

    // Step 1.5: Get twin reads from SNPmers
    let twin_reads = get_twin_reads_from_kmer_info(
        &mut kmer_info,
        &mut blockmer_info,
        &args,
        &temp_dir,
        &temp_dir.join("binary_temp"),
    );

    let clusters = asv_cluster::cluster_reads_by_kmers(&twin_reads, &args, &temp_dir);
    let clusters = asv_cluster::cluster_reads_by_snpmers(&twin_reads, &clusters, &args, &temp_dir);
    let mut consensuses = alignment::align_and_consensus(&twin_reads, clusters, &args, &temp_dir);

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

    // Analyze pileup consensuses 
    let mut low_qual_consensus = alignment::analyze_pileup_consensuses(&mut pileups, &mut consensuses, &quality_error_map, &twin_reads, &args, &temp_dir);
    log_memory_usage(true, "STAGE 2: Analyzed pileup consensus sequences");

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

    alignment::write_consensus_fasta(&low_qual_consensus, &temp_dir.join("low_quality_consensus_sequences.fasta"), "lowqual")
        .expect("Failed to write low_quality_consensus_sequences.fasta");

    // Merge similar consensus sequences based on alignment and depth (using decompressed sequences)
    // This also merges low quality consensuses into high quality ones
    let mut consensuses = alignment::merge_similar_consensuses(&twin_reads, consensuses, low_qual_consensus, &args, &temp_dir);
    log_memory_usage(true, "STAGE 3: Merged similar consensus sequences");

    // Detect and filter chimeric consensus sequences
    if args.skip_chimera_detection {
        log::info!("Skipping chimera detection as per user request.");
        return;
    }
    let chimeras = chimera::detect_chimeras(&mut consensuses, &args);
    let mut consensuses = chimera::filter_chimeras(consensuses, &chimeras);
    log_memory_usage(true, "STAGE 4: Filtered chimeric consensus sequences");

    log::info!("Final consensus count after chimera filtering: {}", consensuses.len());

    // Check for within-ASV heterogeneity
    // if args.phase_heterogeneous {
    //     log::info!("Checking for heterogeneous ASVs to phase...");
    //     alignment::check_asv_heterogeneity(&twin_reads, &mut consensuses, &args);
    // }

    // Refine ASV depths using EM algorithm on read-level mappings
    alignment::refine_asv_depths_with_em(&twin_reads, &mut consensuses, &kmer_info, &args, &temp_dir);
    consensuses.sort_by(|a, b| b.depth.partial_cmp(&a.depth).unwrap());
    log_memory_usage(true, "STAGE 5: Refined ASV depths with EM algorithm");

    log::info!("Final consensus count after EM refinement: {}", consensuses.len());

    // Write final consensus sequences after EM refinement
    let output_dir = std::path::PathBuf::from(&args.output_dir);
    let final_fasta = output_dir.join(ASV_FILE);
    alignment::write_consensus_fasta(&consensuses, &final_fasta, "final")
        .expect(format!("Failed to write {}", ASV_FILE).as_str());
    log::info!("Wrote {} final consensus sequences to {}", consensuses.len(), ASV_FILE);

    debug_consensus_twin_read(&kmer_info, &consensuses, args);

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

fn run_download(args: &cli::DownloadArgs, cli_args: &cli::Cli) {
    // Initialize simple console logger for download command
    let log_spec = format!("{}", cli_args.log_level_filter().to_string());
    let _logger_handle = flexi_logger::Logger::try_with_str(log_spec)
        .expect("Something went wrong with logging")
        .duplicate_to_stderr(flexi_logger::Duplicate::All)
        .format(my_own_format_colored)
        .start()
        .expect("Something went wrong with creating logger");

    log::info!("Starting database download...");

    download::download(&args);

    log::info!("Download complete!");
}

fn initialize_setup_classify(args: &cli::ClassifyArgs, cli_args: &cli::Cli) -> PathBuf {
    let output_dir = if let Some(dir) = args.output_dir.as_ref() {
        Path::new(dir)
    } else {
        Path::new(&args.input_dir)
    };

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

fn get_kmers_and_snpmers(args: &cli::ClusterArgs, output_dir: &PathBuf) -> (types::KmerGlobalInfo, types::BlockmerGlobalInfo) {
    let saved_input = args.input_files == [MAGIC_EXIST_STRING];

    let binary_temp_dir = output_dir.join("binary_temp");
    let snpmer_info_path = binary_temp_dir.join("snpmer_info.bin");

    let kmer_info;
    let blockmer_info;
    if saved_input {
        if !snpmer_info_path.exists() {
            log::error!("No input files provided. See --help for usage.");
            std::process::exit(1);
        }
    }

    let start = Instant::now();
    let (big_snpmer_map, big_blockmer_map) = seq_parse::read_to_split_kmers(args.kmer_size, args.blockmer_length, args.threads, &args);
    log::info!(
        "Time elapsed in for counting k-mers is: {:?}",
        start.elapsed()
    );

    let start = Instant::now();
    if args.use_blockmers {
        blockmer_info = kmer_comp::get_blockmers_inplace_sort(big_blockmer_map, &big_snpmer_map, args.kmer_size, args.blockmer_length, &args);
        log::info!(
            "Time elapsed in for parsing blockmers is: {:?}",
            start.elapsed()
        );
    }
    else{
        blockmer_info = types::BlockmerGlobalInfo::default();
    }

    //kmer_info = kmer_comp::get_snpmers(big_snpmer_map, args.kmer_size, &args);
    kmer_info = kmer_comp::get_snpmers_inplace_sort(big_snpmer_map, args.kmer_size, &args);
    log::info!(
        "Time elapsed in for parsing snpmers is: {:?}",
        start.elapsed()
    );

    return (kmer_info, blockmer_info);
}

fn get_twin_reads_from_kmer_info(
    kmer_info: &mut types::KmerGlobalInfo,
    blockmer_info: &mut types::BlockmerGlobalInfo,
    args: &cli::ClusterArgs,
    _output_dir: &PathBuf,
    _cleaning_temp_dir: &PathBuf,
) -> Vec<types::TwinRead>{
    log::info!("Getting reads...");
    let mut twin_reads_raw = kmer_comp::twin_reads_from_snpmers(kmer_info, blockmer_info, &args);
    twin_reads_raw.sort_by(|a,b| b.est_id.unwrap_or(100.0).partial_cmp(&a.est_id.unwrap_or(100.0)).unwrap());
    return twin_reads_raw;
}

fn debug_consensus_twin_read(kmer_info: &types::KmerGlobalInfo, consensuses: &[ConsensusSequence], args: &cli::ClusterArgs) {

    use std::collections::HashSet;
    let mut snpmer_set = HashSet::default();
    for snpmer_i in kmer_info.snpmer_info.iter(){
        let k = snpmer_i.k as usize;
        let snpmer1 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[0] as u64) << (k-1) );
        let snpmer2 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[1] as u64) << (k-1) );
        snpmer_set.insert(snpmer1);
        snpmer_set.insert(snpmer2);
    }

    for (i,consensus) in consensuses.iter().enumerate() {
        log::trace!("Consensus ID: {}, Index {}, Depth: {}, Length: {}", consensus.id, i, consensus.depth, consensus.decompressed_sequence.as_ref().unwrap().len());
        let tr_rep = seeding::get_twin_read_syncmer(consensus.decompressed_sequence.as_ref().unwrap().clone(), None, args.kmer_size, args.c, args.blockmer_length, &snpmer_set, &FxHashSet::default(), String::new(), args.minimum_base_quality).unwrap();
        let snpmers = tr_rep.snpmers_vec().into_iter().map(|(pos, kmer48)| (pos, decode_kmer48(kmer48, args.kmer_size as u8))).collect::<Vec<_>>();
        log::trace!("SNPmer bases are: {:?}", snpmers);
    }
}