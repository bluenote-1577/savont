use crate::{types::*, utils};
use crate::constants::*;
use std::io::Write;
use std::io::BufWriter;
use std::sync::Mutex;
use crate::cli::Cli;
use std::path::PathBuf;
use minimap2::Aligner;
use rayon::prelude::*;
use bio_seq::prelude::*;
use std::collections::HashMap;

/// Represents a base in the pileup at a reference position
#[derive(Debug, Clone)]
pub enum PileupBase {
    /// A matched/mismatched base with quality score
    Base(u8, u8),
    /// A deletion at this reference position
    Deletion,
    /// An insertion before this reference position (stores inserted bases)
    Insertion(Vec<u8>, Vec<u8>),
}

/// Pileup information at a single reference position
#[derive(Debug, Clone)]
pub struct Pileup {
    pub ref_pos: usize,
    pub ref_base: u8,
    pub bases: Vec<PileupBase>,
    pub alt_posterior: Option<f64>,
}

impl Pileup {
    pub fn new(ref_pos: usize, ref_base: u8) -> Self {
        Self {
            ref_pos,
            ref_base,
            bases: Vec::new(),
            alt_posterior: None,
        }
    }

    pub fn add_base(&mut self, base: u8, quality: u8) {
        self.bases.push(PileupBase::Base(base, quality));
    }

    pub fn add_deletion(&mut self) {
        self.bases.push(PileupBase::Deletion);
    }

    pub fn add_insertion(&mut self, bases: Vec<u8>, qualities: Vec<u8>) {
        self.bases.push(PileupBase::Insertion(bases, qualities));
    }

    pub fn depth(&self) -> usize {
        self.bases.len()
    }

    pub fn depth_nodeletion(&self) -> (usize, usize, usize) {
        let with_bases = self.bases.iter().filter(|b| !matches!(b, PileupBase::Deletion)).count();
        let with_insertions = self.bases.iter().filter(|b| matches!(b, PileupBase::Insertion(_, _))).count();
        let with_deletions = self.bases.iter().filter(|b| matches!(b, PileupBase::Deletion)).count();
        (with_bases, with_insertions, with_deletions)
    }
}

/// Check if there's a homopolymer run of length > 2 around a given position
/// Looks at both directions from the position
fn has_homopolymer_context(seq: &[u8], pos: usize, window: usize) -> bool {
    if seq.is_empty() {
        return false;
    }

    let start = pos.saturating_sub(window);
    let end = (pos + window + 1).min(seq.len());

    if end <= start + 2 {
        return false;
    }

    // Check for runs of length > 2 in the window
    for i in start..=end.saturating_sub(3) {
        if i + 2 < seq.len() && seq[i] == seq[i + 1] && seq[i + 1] == seq[i + 2] {
            return true;
        }
    }

    false
}

/// Calculate adjusted error count from CIGAR alignment
/// Counts mismatches and indels, but only counts indels if they're NOT
/// surrounded by homopolymer runs of length > 2
fn calculate_adjusted_errors(
    cigar: &[(u32, u8)],
    query_seq: &[u8],
    target_seq: &[u8],
    query_start: usize,
    target_start: usize,
) -> usize {
    let mut error_count = 0;
    let buffer = 35;
    let mut query_pos = query_start;
    let mut target_pos = target_start;

    for &(length, op) in cigar {
        let len = length as usize;

        match op {
            0 => {
                // For matches, count actual mismatches
                for _ in 0..len {
                    if query_pos < query_seq.len() && target_pos < target_seq.len() {
                        if query_seq[query_pos] != target_seq[target_pos]  && (query_seq[query_pos] != b'N' && target_seq[target_pos] != b'N') {
                            if query_pos > buffer && query_pos + buffer < query_seq.len() {
                                error_count += 1;
                            }
                        }
                    }
                    query_pos += 1;
                    target_pos += 1;
                }
            }
            100 => {
                // Explicit mismatch
                error_count += len;
                query_pos += len;
                target_pos += len;
            }
            1 => {
                // Insertion in query
                // Only count if NOT in homopolymer context
                let in_homopolymer = has_homopolymer_context(query_seq, query_pos, 2)
                    || has_homopolymer_context(target_seq, target_pos, 2);

                if !in_homopolymer {
                    if query_pos > buffer && query_pos + len + buffer < query_seq.len() {
                        error_count += len;
                    }
                }
                query_pos += len;
            }
            2 => {
                // Deletion in query
                // Only count if NOT in homopolymer context
                let in_homopolymer = has_homopolymer_context(query_seq, query_pos, 2)
                    || has_homopolymer_context(target_seq, target_pos, 2);

                if !in_homopolymer {
                    if target_pos > buffer && target_pos + len + buffer < target_seq.len() {
                        error_count += len;
                    }
                }
                target_pos += len;
            }
            101 => {
                // Soft clip
                query_pos += len;
            }
            102 => {
                // Hard clip - no position change
            }
            _ => {
                // Other operations
                log::warn!("Unexpected CIGAR operation: {}", op as char);
            }
        }
    }

    error_count
}

/// Generate consensus from aligned sequences using POA (Partial Order Alignment)
/// Takes sequences, qualities, and a coverage threshold
/// Returns the consensus sequence as Vec<u8>
fn generate_consensus_poa(
    sequences: &[Vec<u8>],
    qualities: &[Vec<u8>],
    _coverage_threshold: i32,
) -> Vec<u8> {
    if sequences.is_empty() {
        return Vec::new();
    }

    // SPOA implementation
    let mut engine = spoa_rs::AlignmentEngine::new_affine(spoa_rs::AlignmentType::kSW, 3, -6, -8, -4);
    let mut graph = spoa_rs::Graph::new();


    for i in 0..sequences.len() {
        let str_seq = String::from_utf8_lossy(&sequences[i]);
        let u32_qual = qualities[i].iter().map(|&q| q as u32).collect::<Vec<u32>>();
        let (_score, spoa_align) = engine.align(&str_seq, &graph);
        graph.add_alignment_with_weights(spoa_align, &str_seq, &u32_qual);
    }

    let consensus = graph.generate_consensus();
    //return consensus;

    let msa_seqs = graph.generate_msa();
    let msa_string = msa_seqs.join("\n");
    log::debug!("Final Consensus {} has MSA\n{}", sequences.len(), msa_string);

    return consensus.into_bytes();

    // Generate consensus with coverage threshold
    //graph.consensus()
}

pub fn align_and_consensus(twin_reads: &[TwinRead], clusters: Vec<Vec<usize>>, args: &Cli, output_dir: &PathBuf) -> Vec<ConsensusSequence> {
    let max_seqs_consensus = 75;

    // Log which POA implementation is being used
    log::info!("Using SPOA for consensus generation");

    let consensus_seqs = Mutex::new(Vec::new());
    // Implementation of alignment and consensus generation
    clusters.par_iter().enumerate().for_each(|(cluster_idx, cluster)| {
        let mut sequences: Vec<Vec<u8>> = Vec::new();
        let mut qualities: Vec<Vec<u8>> = Vec::new();
        for &read_idx in cluster {
            let twin_read = &twin_reads[read_idx];
            let seq_u8 : Vec<u8> = twin_read.dna_seq.iter().map(|x| x.to_char().to_ascii_uppercase() as u8).collect();
            let qual_u8 = if let Some(qual_seq) = &twin_read.qual_seq {
                qual_seq.iter().map(|x| (x as u8) * 3 + 33).collect()
            } else {
                vec![33; twin_read.dna_seq.len()]
            };

            let bin_size = QUALITY_SEQ_BIN;
            let mut query_quals_u8: Vec<u8> = qual_u8
                .iter()
                .flat_map(|x| vec![*x; bin_size])
                .collect::<Vec<u8>>();

            if query_quals_u8.len() > seq_u8.len() {
                query_quals_u8.truncate(seq_u8.len());
            }
            else if query_quals_u8.len() < seq_u8.len() {
                let last_qual = query_quals_u8[query_quals_u8.len() - 1];
                query_quals_u8.extend(vec![last_qual; seq_u8.len() - query_quals_u8.len()]);
            }
            sequences.push(seq_u8);
            qualities.push(query_quals_u8);
        }
        //let largest_sequence_index = sequences.iter().enumerate().max_by_key(|(i, seq)| seq.len() * (twin_reads[*i].est_id.unwrap() * 100.) as usize).map(|(i, _)| i).unwrap();
        let largest_sequence_index = sequences.iter().enumerate().max_by_key(|(i, seq)| seq.len()).map(|(i, _)| i).unwrap();

        // Create an aligner with appropriate preset
        let aligner = Aligner::builder().map_ont().with_index_threads(args.threads).with_cigar().with_seq(&sequences[largest_sequence_index]).expect("Failed to create aligner");
        let mappings = Mutex::new(Vec::new());

        sequences.iter().enumerate().for_each(|(i, seq)| {
            let alignment = aligner.map(seq, true, false, None, None, None);
            mappings.lock().unwrap().push((i,alignment.unwrap()));
            if i > max_seqs_consensus {
                return;
            }
        });

        // Prepare aligned sequences for POA consensus
        let mut aligned_sequences = Vec::new();
        let mut aligned_qualities = Vec::new();

        let mut mappings = mappings.into_inner().unwrap();
        mappings.sort_by_key(|k| k.0);

        // Add the largest sequence first as seed
        aligned_sequences.push(sequences[largest_sequence_index].clone());
        aligned_qualities.push(qualities[largest_sequence_index].clone());

        for (i, mappings) in mappings.iter(){
            let i = *i;
            if i == largest_sequence_index {
                continue; // Skip the seed sequence
            }

            let best_mapping = &mappings.first().unwrap();
            let cigar_str = &best_mapping.alignment.as_ref().unwrap().cigar_str.as_ref().unwrap();
            log::debug!("Read {}: CIGAR: {}, Query start: {}, Query end: {}, Target start: {}, Target end: {}, Cluster IDX: {}",
            i, cigar_str, best_mapping.query_start, best_mapping.query_end, best_mapping.target_start, best_mapping.target_end, cluster_idx);

            let final_seq;
            let final_qual;
            let qstart;
            let qend;
            if best_mapping.strand == minimap2::Strand::Reverse {
                qstart = sequences[i].len() as i32 - best_mapping.query_end;
                qend = sequences[i].len() as i32 - best_mapping.query_start;
                final_seq = utils::reverse_complement(&sequences[i]);
                final_qual = qualities[i].iter().rev().cloned().collect();
            }
            else{
                qstart = best_mapping.query_start;
                qend = best_mapping.query_end;
                final_seq = sequences[i].clone();
                final_qual = qualities[i].clone();
            }

            let mapped_seq = final_seq[qstart as usize..qend as usize].to_vec();
            let mapped_qual = final_qual[qstart as usize..qend as usize].to_vec();

            aligned_sequences.push(mapped_seq);
            aligned_qualities.push(mapped_qual);

            if aligned_sequences.len() > max_seqs_consensus {
                break;
            }
        }

        let coverage_threshold = (cluster.len().min(max_seqs_consensus) / 10 ) as i32;
        let coverage_threshold = coverage_threshold.max((args.min_cluster_size as i32));

        // Generate consensus using POA (choose implementation based on CLI flag)
        let consensus = generate_consensus_poa(&aligned_sequences, &aligned_qualities, coverage_threshold);
        let buffer = 20;
        if consensus.len() < 2 * buffer {
            log::warn!("Consensus sequence for cluster {} is too short (length {}). Skipping trimming.", cluster_idx, consensus.len());
            return;
        }
        //let consensus = consensus[20..consensus.len()-20].to_vec(); //trim 20bp from each end
        //let msa_string = graph.multiple_sequence_alignment(true)[0..10].iter().map(|x| String::from_utf8_lossy(x).to_string()).collect::<Vec<_>>().join("\n")   ;
        //log::debug!("MSA Cluster {}\n{}", cluster_idx, msa_string);
        let depth = cluster.len();
        consensus_seqs.lock().unwrap().push((cluster_idx, consensus.clone(), depth, cluster.clone()));

        log::info!("Completed alignment for cluster of size {}", cluster.len());
    });

    let mut consensus_seqs = consensus_seqs.into_inner().unwrap();
    consensus_seqs.sort_by_key(|k| k.0);
    let consensus_seqs: Vec<ConsensusSequence> = consensus_seqs.into_iter().map(|(id, seq, depth, cluster)| ConsensusSequence::new(seq, depth, id, cluster)).collect();
    //write consensus sequences to file
    let consensus_path = output_dir.join("consensus_sequences.fasta");
    let mut writer = BufWriter::new(
        std::fs::File::create(consensus_path).expect("Failed to create consensus output file"),
    );
    for (i, consensus) in consensus_seqs.iter().enumerate(){
        let header = format!(">consensus_{}_depth_{}", i, consensus.depth);
        writeln!(writer, "{}", header).expect("Failed to write consensus header");
        writeln!(writer, "{}", String::from_utf8_lossy(&consensus.sequence)).expect("Failed to write consensus sequence");
    }
    consensus_seqs
}

/// Generate pileups for consensus sequences by aligning reads back to them
/// Aligns up to max_seqs_consensus reads per cluster and builds position-wise pileups
pub fn generate_consensus_pileups(
    twin_reads: &[TwinRead],
    consensuses: &[ConsensusSequence],
    args: &Cli,
) -> Vec<Vec<Pileup>> {
    let max_seqs_consensus = 500;

    let pileups = Mutex::new(Vec::new());

    // Process each consensus and its reads in parallel
    consensuses.par_iter().enumerate().for_each(|(cluster_idx, consensus)| {
        let cluster = &consensus.cluster;
        let consensus_seq = &consensus.sequence;

        // Initialize pileup for this consensus
        let mut cluster_pileup: Vec<Pileup> = consensus_seq
            .iter()
            .enumerate()
            .map(|(pos, &base)| Pileup::new(pos, base))
            .collect();

        // Create aligner with this consensus as reference
        let aligner = Aligner::builder()
            .map_ont()
            .with_index_threads(1) // Use 1 thread per aligner since we parallelize over consensuses
            .with_cigar()
            .with_seq(consensus_seq)
            .expect("Failed to create aligner");

        // Align reads from this cluster back to consensus
        let reads_to_align = cluster.len().min(max_seqs_consensus);

        for i in 0..reads_to_align {
            let read_idx = cluster[i];
            let twin_read = &twin_reads[read_idx];

            // Get sequence and quality
            let seq_u8: Vec<u8> = twin_read.dna_seq.iter()
                .map(|x| x.to_char().to_ascii_uppercase() as u8)
                .collect();

            let qual_u8 = if let Some(qual_seq) = &twin_read.qual_seq {
                qual_seq.iter().map(|x| (x as u8) * 3 + 33).collect()
            } else {
                vec![33; twin_read.dna_seq.len()]
            };

            // Bin qualities similar to align_and_consensus
            let bin_size = QUALITY_SEQ_BIN;
            let mut query_quals_u8: Vec<u8> = qual_u8
                .iter()
                .flat_map(|x| vec![*x; bin_size])
                .collect();

            // Adjust quality length to match sequence length
            if query_quals_u8.len() > seq_u8.len() {
                query_quals_u8.truncate(seq_u8.len());
            } else if query_quals_u8.len() < seq_u8.len() {
                let last_qual = query_quals_u8[query_quals_u8.len() - 1];
                query_quals_u8.extend(vec![last_qual; seq_u8.len() - query_quals_u8.len()]);
            }

            // Align read to consensus
            let alignment = aligner.map(&seq_u8, true, false, None, None, None);

            if let Ok(mappings) = alignment {
                if let Some(best_mapping) = mappings.first() {
                    if let Some(ref alignment_info) = best_mapping.alignment {
                        if let Some(ref cigar) = alignment_info.cigar {
                            // Get aligned portion of sequence and quality
                            let mut final_seq = seq_u8.clone();
                            let mut final_qual = query_quals_u8.clone();

                            // Handle reverse complement if needed
                            let reverse; 
                            if best_mapping.strand == minimap2::Strand::Reverse {
                                reverse = true;
                                final_seq = utils::reverse_complement(&seq_u8);
                                final_qual = query_quals_u8.iter().rev().cloned().collect();
                            }
                            else{
                                reverse = false;
                            }

                            // Extract mapped portion
                            let query_start; 
                            let query_end;
                            if reverse {
                                query_start = final_seq.len() - best_mapping.query_end as usize;
                                query_end = final_seq.len() - best_mapping.query_start as usize;
                            } else {
                                query_start = best_mapping.query_start as usize;
                                query_end = best_mapping.query_end as usize;
                            };
                            let mapped_seq = &final_seq[query_start..query_end];
                            let mapped_qual = &final_qual[query_start..query_end];

                            // Process CIGAR to populate pileup
                            let mut ref_pos = best_mapping.target_start as usize;
                            let mut query_pos = 0;

                            for &(length, op) in cigar.iter() {
                                let len = length as usize;

                                match op {
                                    0 => {
                                        // Match or mismatch - add bases to pileup
                                        for j in 0..len {
                                            if ref_pos + j < cluster_pileup.len() && query_pos + j < mapped_seq.len() {
                                                let base = mapped_seq[query_pos + j];
                                                let qual = mapped_qual[query_pos + j];
                                                cluster_pileup[ref_pos + j].add_base(base, qual);
                                            }
                                        }
                                        ref_pos += len;
                                        query_pos += len;
                                    }
                                    1 => {
                                        // Insertion in read - associate with previous ref position
                                        if ref_pos > 0 && ref_pos - 1 < cluster_pileup.len() && query_pos + len <= mapped_seq.len() {
                                            let inserted_bases = mapped_seq[query_pos..query_pos + len].to_vec();
                                            let inserted_qualities = mapped_qual[query_pos..query_pos + len].to_vec();
                                            cluster_pileup[ref_pos - 1].add_insertion(inserted_bases, inserted_qualities);
                                        }
                                        query_pos += len;
                                    }
                                    2 => {
                                        // Deletion in read - add deletion to pileup
                                        for j in 0..len {
                                            if ref_pos + j < cluster_pileup.len() {
                                                cluster_pileup[ref_pos + j].add_deletion();
                                            }
                                        }
                                        ref_pos += len;
                                    }
                                    _ => {
                                        log::warn!("Unexpected CIGAR operation in pileup: {}", op as char);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        log::trace!("Generated pileup for consensus {} with {} positions", cluster_idx, cluster_pileup.len());
        pileups.lock().unwrap().push((cluster_idx, cluster_pileup));
    });

    let mut pileups = pileups.into_inner().unwrap();
    pileups.sort_by_key(|k| k.0);
    let pileups :Vec<Vec<Pileup>> = pileups.into_iter().map(|(_, pileup)| pileup).collect();

    // debug print out

    for (i, pileup) in pileups.iter().enumerate() {
        log::trace!("Pileup for consensus {}:", i);
        for pos in pileup {
            let bases_str: Vec<String> = pos.bases.iter().map(|b| {
                match b {
                    PileupBase::Base(base, qual) => format!("{}({})", *base as char, *qual),
                    PileupBase::Deletion => String::from("D"),
                    PileupBase::Insertion(bases, quals) => format!("I({})({})", String::from_utf8_lossy(bases), String::from_utf8_lossy(quals)),
                }
            }).collect();
            log::trace!("Pos {}: Ref base: {}, Depth: {}, Bases: {}", pos.ref_pos, pos.ref_base as char, pos.depth(), bases_str.join(", "));
        }
    }

    pileups
}

/// Estimate error rate as a function of quality score from pileup data
/// Uses top N clusters and filters positions with <5% error rate
pub fn estimate_quality_error_rates(
    pileups: &[Vec<Pileup>],
    consensuses: &[ConsensusSequence],
    top_frac: f64,
) -> HashMap<u8, f64> {
    // Select top N clusters by depth
    let mut cluster_depths: Vec<(usize, usize)> = consensuses
        .iter()
        .enumerate()
        .map(|(idx, cons)| (idx, cons.depth))
        .collect();
    cluster_depths.sort_by(|a, b| b.1.cmp(&a.1)); // Sort by depth descending

    let top_clusters: Vec<usize> = cluster_depths
        .iter()
        .take((top_frac * cluster_depths.len() as f64).round() as usize)
        .map(|(idx, _)| *idx)
        .collect();

    log::info!("Analyzing quality error rates from top {} clusters", top_clusters.len());

    // Track errors and total bases per quality score
    let mut quality_stats: HashMap<u8, (usize, usize)> = HashMap::new(); // quality -> (errors, total)

    let prior_count = 1;

    for &cluster_idx in &top_clusters {
        if cluster_idx >= pileups.len() {
            continue;
        }

        let cluster_pileup = &pileups[cluster_idx];

        for pileup in cluster_pileup {
            // Calculate error rate at this position
            let mut total_bases = 0;
            let mut error_bases = 0;

            for base_entry in &pileup.bases {
                match base_entry {
                    PileupBase::Base(base, _qual) => {
                        total_bases += 1;
                        if *base != pileup.ref_base {
                            error_bases += 1;
                        }
                    }
                    PileupBase::Deletion => {
                        total_bases += 1;
                        error_bases += 1; // Count deletions as errors
                    }
                    PileupBase::Insertion(_, _) => {
                        // Insertions are associated with previous position, count as error
                        total_bases += 1;
                        error_bases += 1;
                    }
                }
            }

            // Only use positions with < 5% error rate
            if total_bases > 0 {
                let error_fraction = error_bases as f64 / total_bases as f64;
                if error_fraction < 0.05 {
                    // Add quality stats for this position
                    for base_entry in &pileup.bases {
                        if let PileupBase::Base(base, qual) = base_entry {
                            let entry = quality_stats.entry(*qual).or_insert((prior_count, prior_count));
                            entry.1 += 1; // total count
                            if *base != pileup.ref_base {
                                entry.0 += 1; // error count
                            }
                        }
                    }
                }
            }
        }
    }

    // Add prior count


    // Sort qualities for output
    let mut qualities: Vec<u8> = quality_stats.keys().cloned().collect();
    qualities.sort();

    // Calculate overall statistics
    let total_bases_analyzed: usize = quality_stats.values().map(|(_, total)| total).sum();
    let total_errors: usize = quality_stats.values().map(|(errors, _)| errors).sum();
    let overall_error_rate = if total_bases_analyzed > 0 {
        total_errors as f64 / total_bases_analyzed as f64
    } else {
        0.0
    };

    // Output ASCII histogram
    log::debug!("=================================================================");
    log::debug!("Quality Error Rate Histogram (from {} high-confidence positions)", total_bases_analyzed);
    log::debug!("Overall error rate: {:.4}% ({}/{})", overall_error_rate * 100.0, total_errors, total_bases_analyzed);
    log::debug!("=================================================================");

    for qual in qualities {
        if let Some(&(errors, total)) = quality_stats.get(&qual) {
            let error_rate = errors as f64 / total as f64;
            let bar_length = (error_rate * 100.0).round() as usize; // Scale to 100 chars max
            let bar = "#".repeat(bar_length.min(50));
            let spaces = " ".repeat(50_usize.saturating_sub(bar_length));

            log::debug!(
                "Q{:3}: [{}{}] {:6.3}% ({:7}/{:7} errors)",
                qual,
                bar,
                spaces,
                error_rate * 100.0,
                errors,
                total
            );
        }
    }
    log::debug!("=================================================================");

    return quality_stats.iter().map(|(&q, &(e, t))| {
        let rate = if t > 0 { e as f64 / t as f64 } else { 0.0 };
        (q, rate)
    }).collect();
}

/// Log-sum-exp trick for numerically stable log(exp(a) + exp(b))
fn log_sum_exp(log_a: f64, log_b: f64) -> f64 {
    let max = log_a.max(log_b);
    if max.is_infinite() && max.is_sign_negative() {
        return f64::NEG_INFINITY;
    }
    max + ((log_a - max).exp() + (log_b - max).exp()).ln()
}

/// Polish consensus sequences using Bayesian inference with quality-aware error rates
/// Trims low coverage ends and calculates posterior probabilities for each base
pub fn polish_consensuses(
    pileups: &mut Vec<Vec<Pileup>>,
    consensuses: &mut Vec<ConsensusSequence>,
    quality_error_map: &HashMap<u8, f64>,
    args: &Cli,
) -> Vec<ConsensusSequence> {
    let bad_length_threshold = 100;
    let min_coverage_abs = args.min_cluster_size;
    let deletion_insertion_quality = 48u8; // Fixed quality for indels

    // Get error rate for deletions/insertions
    let indel_error_rate = quality_error_map.get(&deletion_insertion_quality)
        .copied()
        .unwrap_or(0.05); // Default 5% if not in map

    log::info!("Polishing {} consensus sequences with min coverage {}", pileups.len(), min_coverage_abs);

    // Select consensuses to debug: most abundant (highest depth), 10th percentile, 90th percentile
    let mut sorted_by_depth: Vec<(usize, usize)> = consensuses
        .iter()
        .enumerate()
        .map(|(idx, cons)| (idx, cons.depth))
        .collect();
    sorted_by_depth.sort_by(|a, b| b.1.cmp(&a.1));
    
    let debug_indices: Vec<usize> = vec![30];

    // Process each consensus
    for (cluster_idx, cluster_pileup) in pileups.iter_mut().enumerate() {
        let min_coverage = (cluster_pileup.iter().map(|p| p.depth()).max().unwrap_or(0) / 2).max(min_coverage_abs);
        if cluster_pileup.is_empty() {
            continue;
        }

        // 1. Trim low coverage ends
        let mut start_idx = 0;
        let mut end_idx = cluster_pileup.len();

        // Find first position with sufficient coverage
        for (i, pileup) in cluster_pileup.iter().enumerate() {
            if pileup.depth() >= min_coverage {
                start_idx = i;
                break;
            }
        }

        // Find last position with sufficient coverage
        for (i, pileup) in cluster_pileup.iter().enumerate().rev() {
            if pileup.depth() >= min_coverage {
                end_idx = i + 1;
                break;
            }
        }

        if start_idx >= end_idx {
            log::warn!("Consensus {} has no positions with sufficient coverage", cluster_idx);
            continue;
        }

        log::debug!("Consensus {}: Trimming from {}-{} to {}-{} with min depth {}",
            cluster_idx, 0, cluster_pileup.len(), start_idx, end_idx, min_coverage);


        // Update pileup to trimmed version
        *cluster_pileup = cluster_pileup[start_idx..end_idx].to_vec();

        // 2. Calculate posterior probabilities for each position
        let mut posterior_probs = Vec::new();

        for pileup in cluster_pileup.iter_mut() {
            let ref_base = pileup.ref_base;

            // Calculate log P(Z | ref = ref_base)
            let mut log_prob_ref = 0.0;

            // Calculate log P(Z | ref != ref_base)
            let mut log_prob_not_ref = 0.0;

            for base_entry in &pileup.bases {
                match base_entry {
                    PileupBase::Base(obs_base, qual) => {
                        let error_rate = quality_error_map.get(qual).copied().unwrap_or(0.05);
                        let accuracy = 1.0 - error_rate;

                        if *obs_base == ref_base {
                            // Observed base matches reference
                            log_prob_ref += accuracy.ln();
                            log_prob_not_ref += error_rate.ln();
                        } else {
                            // Observed base differs from reference
                            log_prob_ref += error_rate.ln();
                            log_prob_not_ref += accuracy.ln();
                        }
                    }
                    PileupBase::Deletion => {
                        // Treat indels as evidence of "not reference"
                        log_prob_ref += indel_error_rate.ln();
                        log_prob_not_ref += (1.0 - indel_error_rate).ln();
                    }
                    PileupBase::Insertion(_, qualities) => {
                        // Add another single evidence since the base before the insertion is not actually correct
                        let first_qual = qualities.first().copied().unwrap_or(deletion_insertion_quality);
                        let error_rate = quality_error_map.get(&first_qual).copied().unwrap_or(0.05);
                        log_prob_not_ref += (1.0 - error_rate).ln();
                        log_prob_ref += error_rate.ln();

                        let error_rates: Vec<f64> = qualities.iter()
                            .map(|q| quality_error_map.get(q).copied().unwrap_or(0.05))
                            .collect();
                        log_prob_ref += error_rates.iter().map(|&er| er.ln()).sum::<f64>();
                        log_prob_not_ref += error_rates.iter().map(|&er| (1.0 - er).ln()).sum::<f64>();
                    }
                }
            }

            // Calculate posterior: P(ref | Z) = P(Z | ref) / (P(Z | ref) + P(Z | not ref))
            // In log space: log P(ref | Z) = log_prob_ref - log_sum_exp(log_prob_ref, log_prob_not_ref)
            let log_normalizer = log_sum_exp(log_prob_ref, log_prob_not_ref);
            let alt_posterior = log_prob_not_ref - log_normalizer;

            if alt_posterior > -30.0 {
                log::debug!("Low posterior probability at consensus {}, position {}: alternate_posterior {:.6}, log_prob_ref {:.4}, log_prob_not_ref {:.4}, depth {}, range {}-{}",
                    cluster_idx, pileup.ref_pos, alt_posterior, log_prob_ref, log_prob_not_ref, pileup.depth(), start_idx, end_idx);
                pileup.alt_posterior = Some(alt_posterior);
            }

            posterior_probs.push(alt_posterior);
        }

        // 3. Debug output for selected consensuses
        let debug_posterior = true;
        if debug_posterior{
            if debug_indices.contains(&cluster_idx) {
                log::debug!("=================================================================");
                log::debug!("Posterior probabilities for consensus {} (depth {})",
                    cluster_idx, consensuses.get(cluster_idx).map(|c| c.depth).unwrap_or(0));
                log::debug!("Position range: {}-{}", start_idx, end_idx);
                log::debug!("=================================================================");

                // Print in chunks of 80 positions for readability
                for chunk_start in (0..posterior_probs.len()).step_by(80) {
                    let chunk_end = (chunk_start + 80).min(posterior_probs.len());

                    // Print position numbers
                    let positions: Vec<String> = (chunk_start..chunk_end)
                        .map(|i| format!("{:4}", i))
                        .collect();
                    log::debug!("Pos:  {}", positions.join(" "));

                    // Print reference bases
                    let ref_bases: Vec<String> = cluster_pileup[chunk_start..chunk_end]
                        .iter()
                        .map(|p| format!("{:>4}", p.ref_base as char))
                        .collect();
                    log::debug!("Ref:  {}", ref_bases.join(" "));

                    // Print depths
                    let depths: Vec<String> = cluster_pileup[chunk_start..chunk_end]
                        .iter()
                        .map(|p| format!("{:?}", p.depth_nodeletion()))
                        .collect();
                    log::debug!("Cov:  {}", depths.join(" "));

                    // Print posterior probabilities
                    let probs: Vec<String> = posterior_probs[chunk_start..chunk_end]
                        .iter()
                        .map(|p| format!("{:4.2}", p))
                        .collect();
                    log::debug!("Post: {}", probs.join(" "));
                    log::debug!("");
                }
                log::debug!("=================================================================");
            }
        }

    }

    let cons_len = consensuses.len();
    for i in 0..cons_len {
        let mut low_confidence_positions = vec![];
        for pileup in &pileups[i]{
            if let Some(_) = pileup.alt_posterior {
                low_confidence_positions.push(pileup.ref_pos);
            }
        }

        if pileups[i].is_empty() {
            log::warn!("Consensus {} has empty pileup after polishing", i);
            continue;
        }
        let left_start = pileups[i].first().map(|p| p.ref_pos).unwrap();
        let right_end = pileups[i].last().map(|p| p.ref_pos).unwrap();

        let start_polish = bad_length_threshold + left_start;
        let end_polish = right_end - bad_length_threshold;

        let low_conf_region_left = low_confidence_positions.iter()
        .filter(|&&pos| pos < start_polish).map(|x| *x).max().unwrap_or(left_start);
        let low_conf_region_right = low_confidence_positions.iter()
        .filter(|&&pos| pos >= end_polish).map(|x| *x).min().unwrap_or(right_end);

        let consensus = &mut consensuses[i];
        if low_conf_region_left > 0 {
            log::debug!("Consensus {}: Masking low-confidence region at start up to position {}", i, low_conf_region_left);
            for pos in 0..=low_conf_region_left {
                consensus.sequence[pos] = b'N';
            }
        }
        if low_conf_region_right < consensus.sequence.len() {
            log::debug!("Consensus {}: Masking low-confidence region at end from position {} to {}", i, low_conf_region_right, consensus.sequence.len());
            for pos in low_conf_region_right..consensus.sequence.len() {
                consensus.sequence[pos] = b'N';
            }
        }
        let pileups = &pileups[i];
        for pileup in pileups.iter(){
            if let Some(_) = pileup.alt_posterior {
                consensus.sequence[pileup.ref_pos] = b'N';
                if pileup.ref_pos > low_conf_region_left && pileup.ref_pos < low_conf_region_right {
                    log::debug!("Consensus {}: Marking position {} as low quality, ends {}-{}", i, pileup.ref_pos, low_conf_region_left, low_conf_region_right);
                    consensus.low_quality_positions.push(pileup.ref_pos);
                }
            }
        }
    }

    let low_quality_consensuses = consensuses.iter().filter(|c| c.low_quality_positions.len() > 1).map(|c| c.clone()).collect::<Vec<_>>();
    log::info!("Low quality consensus sequences: {:?}", &low_quality_consensuses.iter().map(|c| c.id).collect::<Vec<_>>());
    consensuses.retain(|c| c.low_quality_positions.is_empty());

    log::info!("Polishing complete");

    //Write to new fasta

    let output_fasta_path = format!("{}/polished_consensuses.fasta", args.output_dir);
    let mut writer = BufWriter::new(
        std::fs::File::create(&output_fasta_path).expect("Failed to create polished consensus output file"),
    );
    for (i, consensus) in consensuses.iter().enumerate(){
        let header = format!(">polished_consensus_{}_depth_{}", i, consensus.depth);
        writeln!(writer, "{}", header).expect("Failed to write polished consensus header");
        writeln!(writer, "{}", String::from_utf8_lossy(&consensus.sequence)).expect("Failed to write polished consensus sequence");
    }

    return low_quality_consensuses;
}

/// Merge similar consensus sequences based on alignment and depth criteria
/// For each consensus mapping to a higher depth consensus, if the relative depth
/// ratio is less than (1/2)^(NM+1), merge the lower depth consensus into the higher one
pub fn merge_similar_consensuses(
    twin_reads: &[TwinRead],
    consensuses: Vec<ConsensusSequence>,
    args: &Cli,
) -> Vec<ConsensusSequence> {
    if consensuses.is_empty() {
        return consensuses;
    }

    // Build aligner using the first consensus as reference
    let aligner = Aligner::builder()
        .ava_ont()
        .with_index_threads(args.threads)
        .with_cigar()
        .with_index(&format!("{}/polished_consensuses.fasta", args.output_dir), None)
        .expect("Failed to create aligner");

    // Store mappings: (query_idx, target_idx, nm, target_depth)
    let mappings = Mutex::new(Vec::new());
    // Align all consensus sequences to each other in parallel
    consensuses.par_iter().enumerate().for_each(|(query_idx, query_consensus)| {
                    // Align query to target
        let alignment_result = aligner.map(&query_consensus.sequence, true, false, None, None, None);

        if let Ok(alignments) = alignment_result {
            //if let Some(best_mapping) = alignments.iter().max_by_key(|x| x.alignment.as_ref().unwrap().alignment_score.unwrap()){
            for best_mapping in alignments.iter() {
                let target_idx = best_mapping.target_id as usize;
                let target_consensus = &consensuses[target_idx];
                if let Some(alignment) = &best_mapping.alignment {
                    let query_start = best_mapping.query_start as usize;
                    let query_end = best_mapping.query_end as usize;
                    let target_start = best_mapping.target_start as usize;

                    if query_end - query_start < query_consensus.sequence.len() * 3 / 4  || alignment.nm > 30 {
                        continue; // Skip short alignments
                    }

                    // Calculate adjusted error count using CIGAR
                    //dbg!(&alignment.cigar, &alignment.cigar_str);
                    let mut adjusted_errors = if let Some(ref cigar) = alignment.cigar {
                        if best_mapping.strand == minimap2::Strand::Reverse {
                                let rev_query_seq = utils::reverse_complement(&query_consensus.sequence);
                                calculate_adjusted_errors(
                                    cigar,
                                    &rev_query_seq,
                                    &target_consensus.sequence,
                                    query_consensus.sequence.len() - query_end,
                                    target_start,
                                )
                        } else {
                            calculate_adjusted_errors(
                                cigar,
                                &query_consensus.sequence,
                                &target_consensus.sequence,
                                query_start,
                                target_start,
                            )
                        }
                    } else {
                        // Fall back to NM if no CIGAR available
                        alignment.nm as usize
                    };

                    if (alignment.nm as usize) < adjusted_errors {
                        log::debug!("Adjusted errors ({}) greater than NM ({}) for alignment between consensus {} and {}, using NM as adjusted errors",
                            adjusted_errors, alignment.nm, query_idx, target_idx);
                        adjusted_errors = alignment.nm as usize;
                    }

                    mappings.lock().unwrap().push((
                        query_idx,
                        target_idx,
                        adjusted_errors,
                        target_consensus.depth,
                    ));
                }
            }
        }
    });

    let mappings = mappings.into_inner().unwrap();

    // For each query consensus, find the best target (highest depth) to merge with
    let mut merge_map: HashMap<usize, usize> = HashMap::new(); // query_idx -> target_idx

    for query_idx in 0..consensuses.len() {
        // Get all valid mappings for this query
        let mut valid_targets: Vec<(usize, usize, usize)> = mappings
            .iter()
            .filter(|(q_idx, t_idx, nm, t_depth)| {
                if *q_idx != query_idx {
                    return false;
                }
                if *q_idx == *t_idx {
                    return false;
                }

                let query_depth = consensuses[query_idx].depth;
                let target_depth = *t_depth;

                // Calculate relative depth and threshold
                let relative_depth = query_depth as f64 / target_depth as f64;
                let mut threshold = 0.5_f64.powf((*nm as f64) * 0.75 + 1.25);
                if *nm == 0{
                    threshold = 0.999999; // More lenient for perfect matches
                }

                log::debug!("Considering merge: Query {} (depth {}) -> Target {} (depth {}), Adjusted errors {}, Relative depth {:.4}, Threshold {:.4} => {}",
                    query_idx, query_depth, t_idx, target_depth, nm, relative_depth, threshold,
                    relative_depth < threshold);

                relative_depth < threshold
            })
            .map(|(_, t_idx, nm, t_depth)| (*t_idx, *nm, *t_depth))
            .collect();

        // If there are valid targets, choose the one with highest depth
        if !valid_targets.is_empty() {
            valid_targets.sort_by(|a, b| b.2.cmp(&a.2)); // Sort by depth descending
            let best_target = valid_targets[0].0;
            merge_map.insert(query_idx, best_target);
        }
    }

    // Apply merges: create new clusters and consensus sequences
    let mut new_clusters: Vec<Vec<usize>> = consensuses.iter().map(|c| c.cluster.clone()).collect();
    let mut merged_into: HashMap<usize, usize> = HashMap::new(); // Maps original idx to final idx

    // Resolve merge chains (A->B, B->C should result in A->C, B->C)
    for query_idx in 0..consensuses.len() {
        if let Some(&target_idx) = merge_map.get(&query_idx) {
            let mut final_target = target_idx;
            while let Some(&next_target) = merge_map.get(&final_target) {
                final_target = next_target;
            }
            merged_into.insert(query_idx, final_target);
        }
    }

    // Perform the merges
    for (&query_idx, &target_idx) in &merged_into {
        log::info!(
            "Merging consensus {} (depth {}) into consensus {} (depth {})",
            query_idx,
            consensuses[query_idx].depth,
            target_idx,
            consensuses[target_idx].depth
        );

        // Move all reads from query cluster to target cluster
        let reads_to_move = new_clusters[query_idx].clone();
        new_clusters[target_idx].extend(reads_to_move);
        new_clusters[query_idx].clear();
    }

    // Build new consensus sequences with updated depths
    let mut new_consensuses = Vec::new();

    for (idx, consensus) in consensuses.into_iter().enumerate() {
        if !new_clusters[idx].is_empty() {
            let new_depth = new_clusters[idx].len();
            let new_cluster = new_clusters[idx].clone();
            new_consensuses.push(ConsensusSequence::new(consensus.sequence, new_depth, consensus.id, new_cluster));
        }
    }

    log::info!(
        "Consensus merging complete: {} -> {} consensuses",
        new_clusters.len(),
        new_consensuses.len()
    );


    let output_dir = std::path::PathBuf::from(&args.output_dir);
    let final_file = output_dir.join("final_clusters_merged.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&final_file).unwrap());

    for (snpmer_rep_id, consensus) in new_consensuses.iter().enumerate() {
        let cluster = &consensus.cluster;
        let representative = cluster[0];
        writeln!(
            writer,
            "final_cluster_{}\tsize_{}\trepresentative_{}\tmembers\n{}",
            snpmer_rep_id,
            cluster.len(),
            representative,
            cluster.iter().map(|x| format!("{} {}", &twin_reads[*x].id, &twin_reads[*x].est_id.unwrap())).collect::<Vec<_>>().join("\n")
        ).unwrap();
    }

    let final_file = output_dir.join("final_consensus_sequences.fasta");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&final_file).unwrap());

    for (id, consensus) in new_consensuses.iter().enumerate() {
        writeln!(writer, ">id_{}_consensus_{}_depth_{}\n{}", id, consensus.id, consensus.depth, String::from_utf8_lossy(&consensus.sequence)).unwrap();
    }

    new_consensuses
}
