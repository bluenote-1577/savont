use crate::{seeding, types::*, utils};
use crate::constants::*;
use std::io::Write;
use std::sync::Mutex;
use crate::cli::ClusterArgs as Cli;
use std::path::PathBuf;
use minimap2::Aligner;
use rayon::prelude::*;
use bio_seq::prelude::*;
use std::collections::HashMap;
use fxhash::{FxHashMap, FxHashSet};
use crate::asv_cluster::find_compatible_candidates;
use crate::kmer_comp;

/// Represents a base in the pileup at a reference position
#[derive(Debug, Clone)]
pub enum PileupBase {
    /// A matched/mismatched base with quality score and homopolymer length
    /// (base, quality, hp_length from this read)
    Base(u8, u8, u8),
    /// A deletion at this reference position
    Deletion,
    /// An insertion before this reference position
    /// Vec of (base, quality, hp_length) tuples
    Insertion(Vec<(u8, u8, u8)>),
}

/// Pileup information at a single reference position
#[derive(Debug, Clone)]
pub struct Pileup {
    pub ref_pos: usize,        // Position in HPC consensus
    pub ref_base: u8,          // HPC base
    pub ref_hp_length: u8,     // Modal HP length from aligned reads
    pub bases: Vec<PileupBase>,
    pub alt_posterior: Option<f64>,
}

impl Pileup {
    pub fn new(ref_pos: usize, ref_base: u8, ref_hp_length: u8) -> Self {
        Self {
            ref_pos,
            ref_base,
            ref_hp_length,
            bases: Vec::new(),
            alt_posterior: None,
        }
    }

    pub fn add_base(&mut self, base: u8, quality: u8, hp_length: u8) {
        self.bases.push(PileupBase::Base(base, quality, hp_length));
    }

    pub fn add_deletion(&mut self) {
        self.bases.push(PileupBase::Deletion);
    }

    pub fn add_insertion(&mut self, insertion_data: Vec<(u8, u8, u8)>) {
        self.bases.push(PileupBase::Insertion(insertion_data));
    }

    pub fn depth(&self) -> usize {
        self.bases.len()
    }

    pub fn depth_nodeletion(&self) -> (usize, usize, usize) {
        let with_bases = self.bases.iter().filter(|b| !matches!(b, PileupBase::Deletion)).count();
        let with_insertions = self.bases.iter().filter(|b| matches!(b, PileupBase::Insertion(_))).count();
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
/// Use gap-collapsed NM 
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
                        if len < 10{
                            error_count += 1;
                        }
                        else{
                            error_count += len;
                        }
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
                        if len < 10{
                            error_count += 1;
                        }
                        else{
                            error_count += len;
                        }
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
    cluster_idx: usize,
) -> Vec<u8> {
    if sequences.is_empty() {
        return Vec::new();
    }

    // SPOA implementation
    let mut engine = spoa_rs::AlignmentEngine::new_affine(spoa_rs::AlignmentType::kSW, 3, -8, -6, -6);
    let mut graph = spoa_rs::Graph::new();


    for i in 0..sequences.len() {
        let str_seq = String::from_utf8_lossy(&sequences[i]);
        let u32_qual = qualities[i].iter().map(|&q| q as u32).collect::<Vec<u32>>();
        let (_score, spoa_align) = engine.align(&str_seq, &graph);
        graph.add_alignment_with_weights(spoa_align, &str_seq, &u32_qual);
    }

    let consensus = graph.generate_consensus();

    let msa_seqs = &graph.generate_msa()[0..5];
    let msa_string = msa_seqs.join("\n");
    log::trace!("Final Consensus {} has length {} with MSA\n{}", cluster_idx, consensus.len(), msa_string);

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
        //let largest_sequence_index = sequences.iter().enumerate().max_by_key(|(_, seq)| seq.len()).map(|(i, _)| i).unwrap();
        // largest sequence index is the 90th percentile length sequence to avoid chimeras
        let mut lengths_and_i: Vec<(_,_)> = sequences.iter().enumerate().map(|(i, seq)| (seq.len(), i)).collect();
        lengths_and_i.sort_by_key(|k| k.0);
        let largest_sequence_index = lengths_and_i[(lengths_and_i.len() as f64 * 0.9) as usize].1;


        // Create an aligner with appropriate preset
        let aligner = Aligner::builder().map_ont().with_index_threads(args.threads).with_cigar().with_seq(&sequences[largest_sequence_index]).expect("Failed to create aligner");
        let mappings = Mutex::new(Vec::new());

        sequences.iter().enumerate().for_each(|(i, seq)| {
            let alignment = aligner.map(seq, true, false, None, None, None);
            if alignment.is_err() {
                log::debug!("No alignment found for read {} in cluster {}", i, cluster_idx);
                return;
            }
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
            if mappings.is_empty() {
                log::debug!("No alignment found for read {} in cluster {}", i, cluster_idx);
                continue;
            }

            let best_mapping = &mappings.first().unwrap();
            let cigar_str = &best_mapping.alignment.as_ref().unwrap().cigar_str.as_ref().unwrap();
            log::trace!("Read {} len {}: CIGAR: {}, Query start: {}, Query end: {}, Target start: {}, Target end: {}, Cluster IDX: {}",
            i, best_mapping.query_len.unwrap(), cigar_str, best_mapping.query_start, best_mapping.query_end, best_mapping.target_start, best_mapping.target_end, cluster_idx);

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

        // HPC compress all sequences before POA
        let mut hpc_aligned_sequences = Vec::new();
        let mut hpc_aligned_qualities = Vec::new();
        for i in 0..aligned_sequences.len() {
            let (hpc_seq, hpc_qual, _hp_lens) = utils::homopolymer_compress_with_quality(
                &aligned_sequences[i],
                &aligned_qualities[i],
                args.use_hpc,
            );
            hpc_aligned_sequences.push(hpc_seq);
            hpc_aligned_qualities.push(hpc_qual);
        }

        let coverage_threshold = (cluster.len().min(max_seqs_consensus) / 10 ) as i32;
        let coverage_threshold = coverage_threshold.max(args.min_cluster_size as i32);

        // Generate consensus using POA (working with HPC sequences)
        let hpc_consensus = generate_consensus_poa(&hpc_aligned_sequences, &hpc_aligned_qualities, coverage_threshold, cluster_idx);

        // Compress the consensus again to ensure it's fully HPC
        let (hpc_consensus_seq, _) = utils::homopolymer_compress(&hpc_consensus, args.use_hpc);

        let buffer = 20;
        if hpc_consensus_seq.len() < 2 * buffer {
            log::warn!("HPC consensus sequence for cluster {} is too short (length {}). Skipping trimming.", cluster_idx, hpc_consensus_seq.len());
            return;
        }
        //let consensus = consensus[20..consensus.len()-20].to_vec(); //trim 20bp from each end
        //let msa_string = graph.multiple_sequence_alignment(true)[0..10].iter().map(|x| String::from_utf8_lossy(x).to_string()).collect::<Vec<_>>().join("\n")   ;
        //log::debug!("MSA Cluster {}\n{}", cluster_idx, msa_string);
        let depth = cluster.len();
        // HP lengths will be calculated from pileup alignments, use placeholder (1) for now
        let placeholder_hp_lengths = vec![1u8; hpc_consensus_seq.len()];
        consensus_seqs.lock().unwrap().push((cluster_idx, hpc_consensus_seq.clone(), placeholder_hp_lengths, depth, cluster.clone()));

        log::debug!("Completed alignment for cluster of size {}", cluster.len());
    });

    let mut consensus_seqs = consensus_seqs.into_inner().unwrap();
    consensus_seqs.sort_by_key(|k| k.0);
    let consensus_seqs: Vec<ConsensusSequence> = consensus_seqs.into_iter().map(|(id, seq, hp_lens, depth, cluster)| ConsensusSequence::new(seq, hp_lens, depth, id, cluster)).collect();

    // Write HPC consensus sequences to file
    let consensus_path = output_dir.join("consensus_sequences.fasta");
    write_consensus_fasta(&consensus_seqs, &consensus_path, "initial")
        .expect("Failed to write consensus_sequences.fasta");
    log::info!("Wrote {} consensus sequences to consensus_sequences.fasta", consensus_seqs.len());

    consensus_seqs
}

/// Generate pileups for consensus sequences by aligning reads back to them
/// Aligns up to max_seqs_consensus reads per cluster and builds position-wise pileups
pub fn generate_consensus_pileups(
    twin_reads: &[TwinRead],
    consensuses: &[ConsensusSequence],
    args: &Cli,
) -> Vec<Vec<Pileup>> {
    let max_seqs_consensus = MAX_SEQS_CONSENSUS;

    let pileups = Mutex::new(Vec::new());

    // Process each consensus and its reads in parallel
    consensuses.par_iter().enumerate().for_each(|(cluster_idx, consensus)| {
        let cluster = &consensus.cluster;
        let consensus_seq = &consensus.sequence;
        let consensus_hp_lengths = &consensus.hp_lengths;

        // Initialize pileup for this consensus with placeholder ref_hp_length (will be updated later)
        let mut cluster_pileup: Vec<Pileup> = consensus_seq
            .iter()
            .enumerate()
            .map(|(pos, &base)| Pileup::new(pos, base, consensus_hp_lengths[pos]))
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

            // HPC compress the read before aligning to HPC consensus
            let (hpc_seq, hpc_qual, hp_lens) = utils::homopolymer_compress_with_quality(&seq_u8, &query_quals_u8, args.use_hpc);

            // Align HPC read to HPC consensus
            let alignment = aligner.map(&hpc_seq, true, false, None, None, None);

            if let Ok(mappings) = alignment {
                if let Some(best_mapping) = mappings.first() {
                    if let Some(ref alignment_info) = best_mapping.alignment {
                        if let Some(ref cigar) = alignment_info.cigar {
                            // Get aligned portion of HPC sequence, quality, and HP lengths
                            let final_seq;
                            let final_qual;
                            let final_hp_lens;

                            // Handle reverse complement if needed
                            let reverse;
                            if best_mapping.strand == minimap2::Strand::Reverse {
                                reverse = true;
                                final_seq = utils::reverse_complement(&hpc_seq);
                                final_qual = hpc_qual.iter().rev().cloned().collect();
                                final_hp_lens = hp_lens.iter().rev().cloned().collect();
                            }
                            else{
                                final_seq = hpc_seq;
                                final_qual = hpc_qual;
                                final_hp_lens = hp_lens;
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
                            let mapped_hp_lens = &final_hp_lens[query_start..query_end];

                            // Process CIGAR to populate pileup
                            let mut ref_pos = best_mapping.target_start as usize;
                            let mut query_pos = 0;

                            for &(length, op) in cigar.iter() {
                                let len = length as usize;

                                match op {
                                    0 => {
                                        // Match or mismatch - add bases to pileup with HP lengths
                                        for j in 0..len {
                                            if ref_pos + j < cluster_pileup.len() && query_pos + j < mapped_seq.len() {
                                                let base = mapped_seq[query_pos + j];
                                                let qual = mapped_qual[query_pos + j];
                                                let hp_len = mapped_hp_lens[query_pos + j];
                                                cluster_pileup[ref_pos + j].add_base(base, qual, hp_len);
                                            }
                                        }
                                        ref_pos += len;
                                        query_pos += len;
                                    }
                                    1 => {
                                        // Insertion in read - associate with previous ref position
                                        if ref_pos > 0 && ref_pos - 1 < cluster_pileup.len() && query_pos + len <= mapped_seq.len() {
                                            let mut insertion_data = Vec::new();
                                            for j in 0..len {
                                                let base = mapped_seq[query_pos + j];
                                                let qual = mapped_qual[query_pos + j];
                                                let hp_len = mapped_hp_lens[query_pos + j];
                                                insertion_data.push((base, qual, hp_len));
                                            }
                                            cluster_pileup[ref_pos - 1].add_insertion(insertion_data);
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
    let mut pileups: Vec<Vec<Pileup>> = pileups.into_iter().map(|(_, pileup)| pileup).collect();

    // Calculate modal HP lengths from aligned reads for each position
    for pileup_vec in pileups.iter_mut() {
        for pileup in pileup_vec.iter_mut() {
            // Collect all HP lengths from aligned bases at this position
            let mut hp_lengths: Vec<u8> = Vec::new();
            for base_entry in &pileup.bases {
                if let PileupBase::Base(_base, _qual, hp_len) = base_entry {
                    hp_lengths.push(*hp_len);
                }
            }

            // Calculate modal HP length
            if !hp_lengths.is_empty() {
                // Count occurrences of each HP length
                let mut counts: std::collections::HashMap<u8, usize> = std::collections::HashMap::new();
                for &hp_len in &hp_lengths {
                    *counts.entry(hp_len).or_insert(0) += 1;
                }

                // Find the most common HP length
                let modal_hp_length = counts.iter()
                    .max_by_key(|(_, count)| *count)
                    .map(|(hp_len, _)| *hp_len)
                    .unwrap_or(1);

                pileup.ref_hp_length = modal_hp_length;
            } else {
                // No bases aligned, keep placeholder value
                pileup.ref_hp_length = 1;
            }
        }
    }

    // debug print out

    if log::log_enabled!(log::Level::Trace){
        for (i, pileup) in pileups.iter().enumerate() {
            log::trace!("Pileup for consensus {}:", i);
            for pos in pileup {
                let bases_str: Vec<String> = pos.bases.iter().map(|b| {
                    match b {
                        PileupBase::Base(base, qual, hp_len) => format!("{}(q={},hp={})", *base as char, *qual, *hp_len),
                        PileupBase::Deletion => String::from("D"),
                        PileupBase::Insertion(data) => {
                            let bases: Vec<u8> = data.iter().map(|(b, _, _)| *b).collect();
                            format!("I({})", String::from_utf8_lossy(&bases))
                        },
                    }
                }).collect();
                log::trace!("Pos {}: Ref base: {}, Ref HP len: {}, Depth: {}, Bases: {}",
                    pos.ref_pos, pos.ref_base as char, pos.ref_hp_length, pos.depth(), bases_str.join(", "));
            }
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
                    PileupBase::Base(base, _qual, _hp_len) => {
                        total_bases += 1;
                        if *base != pileup.ref_base {
                            error_bases += 1;
                        }
                    }
                    PileupBase::Deletion => {
                        total_bases += 1;
                        error_bases += 1; // Count deletions as errors
                    }
                    PileupBase::Insertion(_) => {
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
                        if let PileupBase::Base(base, qual, _hp_len) = base_entry {
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

/// Write cluster information to TSV file
/// Format: cluster_id, size, representative, members (one per line with ID and est_id)
pub fn write_clusters_tsv(
    consensuses: &[ConsensusSequence],
    twin_reads: &[TwinRead],
    output_path: &std::path::Path,
    prefix: &str,
) -> std::io::Result<()> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(output_path)?);

    for (cluster_id, consensus) in consensuses.iter().enumerate() {
        let cluster = &consensus.cluster;
        if cluster.is_empty() {
            continue;
        }

        let representative = cluster[0];
        writeln!(
            writer,
            "{}_cluster_{}\tsize_{}\trepresentative_{}\tmembers\n{}",
            prefix,
            cluster_id,
            cluster.len(),
            representative,
            cluster.iter().map(|x| format!("{} {}", &twin_reads[*x].id, &twin_reads[*x].est_id.unwrap_or(100.))).collect::<Vec<_>>().join("\n")
        )?;
    }

    Ok(())
}

/// Write consensus sequences to a FASTA file
/// Uses decompressed sequences if available, otherwise uses HPC sequences
pub fn write_consensus_fasta(
    consensuses: &[ConsensusSequence],
    output_path: &std::path::Path,
    prefix: &str,
) -> std::io::Result<()> {
    let mut writer = std::io::BufWriter::new(std::fs::File::create(output_path)?);

    for (i, consensus) in consensuses.iter().enumerate() {
        let mut cons_clone = consensus.clone();
        cons_clone.decompress();
        let consensus_seq = cons_clone.decompressed_sequence.clone().unwrap();
        let start_non = consensus_seq.iter().enumerate().find(|&(_i, &b)| b != b'N').map(|(i, _)| i).unwrap_or(0);
        let end_non = consensus_seq.iter().enumerate().rfind(|&(_i, &b)| b != b'N').map(|(i, _)| i).unwrap_or(consensus_seq.len());
        let header = format!(">{}_consensus_{}_depth_{} debug_id:{} chimera_score:{} unambiguous_read_assignments:{} ambig_read_assignments:{} num_align_leq_10_mismatches:{}", 
            prefix, i, consensus.depth + consensus.appended_depth, consensus.id, consensus.chimera_score.unwrap_or(0), consensus.unambig_best_read_map_count.unwrap_or(0), 
            consensus.ambig_read_map_count.unwrap_or(0), consensus.num_map_leq_10nm.unwrap_or(0));
        writeln!(writer, "{}", header)?;

        // Use decompressed sequence if available, otherwise use HPC sequence
        let sequence = String::from_utf8_lossy(&consensus_seq[start_non..=end_non]);

        writeln!(writer, "{}", sequence)?;
    }

    Ok(())
}

/// Polish consensus sequences using Bayesian inference with quality-aware error rates
/// Trims low coverage ends and calculates posterior probabilities for each base
pub fn analyze_pileup_consensuses(
    pileups: &mut Vec<Vec<Pileup>>,
    consensuses: &mut Vec<ConsensusSequence>,
    quality_error_map: &HashMap<u8, f64>,
    twin_reads: &[TwinRead],
    args: &Cli,
    temp_dir: &PathBuf,
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
        let min_coverage = (cluster_pileup.iter().map(|p| p.depth()).max().unwrap_or(0) / 3).max(min_coverage_abs);
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

        log::trace!("Consensus {}: Trimming from {}-{} to {}-{} with min depth {}",
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
                    PileupBase::Base(obs_base, qual, _hp_len) => {
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
                    PileupBase::Insertion(insertion_data) => {
                        // Add another single evidence since the base before the insertion is not actually correct
                        let first_qual = insertion_data.first().map(|(_, q, _)| *q).unwrap_or(deletion_insertion_quality);
                        let error_rate = quality_error_map.get(&first_qual).copied().unwrap_or(0.05);
                        log_prob_not_ref += (1.0 - error_rate).ln();
                        log_prob_ref += error_rate.ln();

                        let error_rates: Vec<f64> = insertion_data.iter().take(0)
                            .map(|(_, q, _)| quality_error_map.get(q).copied().unwrap_or(0.05))
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

            //if alt_posterior > -30.0 {
            if alt_posterior > -args.posterior_threshold_ln{
                log::debug!("Low posterior probability at consensus {}, covs: {:?},  position {}: alternate_posterior {:.6}, log_prob_ref {:.4}, log_prob_not_ref {:.4}, depth {}, range {}-{}",
                    cluster_idx, pileup.depth_nodeletion(),  pileup.ref_pos, alt_posterior, log_prob_ref, log_prob_not_ref, pileup.depth(), start_idx, end_idx);
                let ref_count = pileup.bases.iter().filter(|b| {
                    if let PileupBase::Base(base, _, _) = b {
                        *base == ref_base
                    } else {
                        false
                    }
                }).count();
                let non_ref_base_count = pileup.depth() - ref_count;
                log::debug!("    Reference base: {}, Ref count: {}, Non-ref count: {}", ref_base as char, ref_count, non_ref_base_count);
                //print top 20 bases
                let mut dbg_string = String::from("Bases: ");
                for base_entry in pileup.bases.iter().take(20) {
                    match base_entry {
                        PileupBase::Base(obs_base, qual, hp_len) => {
                            dbg_string.push_str(&format!(" Base: {} (q={}, hp={}) ", *obs_base as char, *qual, *hp_len));
                        }
                        PileupBase::Deletion => {
                            dbg_string.push_str(" Deletion ");
                        }
                        PileupBase::Insertion(insertion_data) => {
                            let bases: Vec<u8> = insertion_data.iter().map(|(b, _, _)| *b).collect();
                            dbg_string.push_str(&format!(" Insertion: {} ", String::from_utf8_lossy(&bases)));
                        }
                    }
                }
                log::debug!("{}", dbg_string);
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
            log::trace!("Consensus {}: Masking low-confidence region at start up to position {}", i, low_conf_region_left);
            for pos in 0..=low_conf_region_left {
                consensus.sequence[pos] = b'N';
            }
        }
        if low_conf_region_right < consensus.sequence.len() {
            log::trace!("Consensus {}: Masking low-confidence region at end from position {} to {}", i, low_conf_region_right, consensus.sequence.len());
            for pos in low_conf_region_right..consensus.sequence.len() {
                consensus.sequence[pos] = b'N';
            }
        }
        let pileups = &pileups[i];
        for pileup in pileups.iter(){
            if let Some(_) = pileup.alt_posterior {
                if args.mask_low_quality{
                    consensus.sequence[pileup.ref_pos] = b'N';
                }
                if pileup.ref_pos > low_conf_region_left && pileup.ref_pos < low_conf_region_right {
                    log::debug!("Consensus {}: Marking position {} as low quality, ends {}-{}", i, pileup.ref_pos, low_conf_region_left, low_conf_region_right);
                    consensus.low_quality_positions.push(pileup.ref_pos);
                }
            }
        }
    }

    // Write clusters before filtering out low quality consensuses
    let prefilter_file = temp_dir.join("clusters_before_quality_filter.tsv");
    write_clusters_tsv(consensuses, twin_reads, &prefilter_file, "prefilter")
        .expect("Failed to write clusters_before_quality_filter.tsv");
    log::info!("Wrote cluster information before filtering to clusters_before_quality_filter.tsv");

    let low_quality_consensuses = consensuses.iter().filter(|c| lq_criteria(c, args)).map(|c| c.clone()).collect::<Vec<_>>();
    log::info!("Low quality consensus sequences: {:?}", &low_quality_consensuses.iter().map(|c| c.id).collect::<Vec<_>>());
    consensuses.retain(|c| !lq_criteria(c, args));

    log::info!("Polishing complete");

    //Write to new fasta



    return low_quality_consensuses;
}

fn lq_criteria(consensus: &ConsensusSequence, args: &Cli) -> bool {
    (consensus.low_quality_positions.len() > 0) && 
    (consensus.depth / ((consensus.low_quality_positions.len() * consensus.low_quality_positions.len())) < args.n_depth_cutoff)
}

fn remove_similar_seqs_kmers(mut consensuses: Vec<ConsensusSequence>) -> Vec<ConsensusSequence> {
    let mut kmer_index = std::collections::HashMap::new();
    let mut filtered_consensuses: Vec<ConsensusSequence> = Vec::new();
    let mut consensus_id_to_minis = std::collections::HashMap::new();
    let adapter_buffer = 25;

    for (i, consensus) in consensuses.iter().enumerate() {
        if consensus.sequence.len() < 100{
            continue;
        }
        let mut minimizers = vec![];
        let mut positions = vec![];
        seeding::minimizer_seeds_positions(&consensus.sequence[adapter_buffer..consensus.sequence.len()-adapter_buffer], &mut minimizers, &mut positions, 10, 21);
        for &mini in minimizers.iter() {
            kmer_index.entry(mini).or_insert_with(Vec::new).push(i);
        }
        consensus_id_to_minis.insert(i, minimizers);
    }

    for (&enum_id, minimizers) in consensus_id_to_minis.iter() {
        let mut possible_greater_ids = std::collections::HashSet::new();
        let mut first = true;
        for mini in minimizers {
            if first{
                if let Some(ids) = kmer_index.get(mini) {
                    for id in ids {
                        if consensuses[*id].depth / 2 > consensuses[enum_id].depth {
                            possible_greater_ids.insert(*id);
                        }
                    }
                }
            }
            else{
                if let Some(ids) = kmer_index.get(mini) {
                    let id_set = ids.iter().cloned().collect::<std::collections::HashSet<usize>>();
                    possible_greater_ids = possible_greater_ids.intersection(&id_set).cloned().collect();
                }
            }
            first = false;
        }
        if possible_greater_ids.is_empty() {
            filtered_consensuses.push(std::mem::take(&mut consensuses[enum_id]));
        }
    }

    return filtered_consensuses;
}

/// Merge similar consensus sequences based on alignment and depth criteria
/// For each consensus mapping to a higher depth consensus, if the relative depth
/// ratio is less than (1/2)^(NM+1), merge the lower depth consensus into the higher one
pub fn merge_similar_consensuses(
    twin_reads: &[TwinRead],
    consensuses: Vec<ConsensusSequence>,
    low_qual_consensuses: Vec<ConsensusSequence>,
    args: &Cli,
    temp_dir: &PathBuf,
) -> Vec<ConsensusSequence> {
    if consensuses.is_empty() {
        return consensuses;
    }

    
    // Look at [35,length-35] bases to avoid adapters. Take all k-mers. Remove subsetted reads at lower depth.
    log::info!("Removing duplicate consensus sequences based on k-mer similarity");
    let prev_size = consensuses.len();
    let consensuses = remove_similar_seqs_kmers(consensuses);
    log::info!("Reduced consensus sequences from {} to {} after k-mer based deduplication", prev_size, consensuses.len());

    // Write polished consensus sequences to FASTA for indexing
    let output_fasta_path = temp_dir.join("polished_consensuses.fasta");
    write_consensus_fasta(&consensuses, &output_fasta_path, "polished")
        .expect("Failed to write polished_consensuses.fasta");
    log::info!("Wrote {} polished consensus sequences to polished_consensuses.fasta", consensuses.len());


    // Build aligner using the first consensus as reference
    let mut aligner = Aligner::builder()
        .lrhq()
        .with_index_threads(args.threads)
        .with_cigar()
        .with_index(output_fasta_path.to_str().unwrap(), None)
        .expect("Failed to create aligner");

    aligner.mapopt.set_no_diag();
    
    // Store mappings: (query_idx, target_idx, nm, target_depth)
    let mappings = Mutex::new(Vec::new());

    // First, merge low quality consensuses into high quality consensuses
    log::info!("Merging {} low quality consensuses into high quality consensuses", low_qual_consensuses.len());
    let low_qual_mappings = Mutex::new(Vec::new());

    low_qual_consensuses.par_iter().enumerate().for_each(|(low_qual_idx, low_qual_consensus)| {
        let query_seq = low_qual_consensus.decompressed_sequence.as_ref()
            .expect("Low quality consensus must be decompressed");

        // Align to all high quality consensuses to find best match
        let alignment_result = aligner.map(query_seq, true, false, None, None, None);

        if let Ok(alignments) = alignment_result {
            //print mapqs:
            // Find the best primary alignment
            if let Some(best_mapping) = alignments.first() {
                let target_idx = best_mapping.target_id as usize;

                if let Some(alignment) = &best_mapping.alignment {

                    if alignment.nm > 10 {
                        return; // Skip high error alignments
                    }

                    log::debug!("Low quality consensus {} (id={}, depth={}) maps to consensus {} (depth = {}) with NM={}",
                        low_qual_idx, low_qual_consensus.id, low_qual_consensus.depth, target_idx, consensuses[target_idx].depth, alignment.nm);

                    // Store the mapping: (low_qual_consensus, target_idx)
                    low_qual_mappings.lock().unwrap().push((low_qual_idx, target_idx));
                }
            }
        }
    });

    // Merge low quality consensus reads into their target consensuses
    let low_qual_mappings = low_qual_mappings.into_inner().unwrap();
    let mut consensuses = consensuses;

    for (query_idx, target_idx) in low_qual_mappings {
        // Merge the clusters
        let low_qual_consensus = &low_qual_consensuses[query_idx];
        consensuses[target_idx].appended_depth += low_qual_consensus.depth;
    }

    log::info!("Merging similar consensus sequences based on alignment and depth criteria");

    // Align all consensus sequences to each other in parallel (using decompressed sequences)
    consensuses.par_iter().enumerate().for_each(|(query_idx, query_consensus)| {
        // Use decompressed sequence for alignment
        let query_seq = query_consensus.decompressed_sequence.as_ref()
            .expect("Consensus sequence must be decompressed before merging");

        // Align query to target
        let alignment_result = aligner.map(query_seq, true, false, None, None, None);

        if let Ok(alignments) = alignment_result {
            //if let Some(best_mapping) = alignments.iter().max_by_key(|x| x.alignment.as_ref().unwrap().alignment_score.unwrap()){
            for best_mapping in alignments.iter() {
                let target_idx = best_mapping.target_id as usize;
                let target_consensus = &consensuses[target_idx];
                let target_seq = target_consensus.decompressed_sequence.as_ref()
                    .expect("Consensus sequence must be decompressed before merging");

                if let Some(alignment) = &best_mapping.alignment {
                    let query_start = best_mapping.query_start as usize;
                    let query_end = best_mapping.query_end as usize;
                    let target_start = best_mapping.target_start as usize;

                    if query_end - query_start < query_seq.len() * 3 / 4  || alignment.nm > 30 {
                        continue; // Skip short alignments
                    }

                    // Calculate adjusted error count using CIGAR (on decompressed sequences)
                    //dbg!(&alignment.cigar, &alignment.cigar_str);
                    let mut adjusted_errors = if let Some(ref cigar) = alignment.cigar {
                        if best_mapping.strand == minimap2::Strand::Reverse {
                                let rev_query_seq = utils::reverse_complement(query_seq);
                                calculate_adjusted_errors(
                                    cigar,
                                    &rev_query_seq,
                                    target_seq,
                                    query_seq.len() - query_end,
                                    target_start,
                                )
                        } else {
                            calculate_adjusted_errors(
                                cigar,
                                query_seq,
                                target_seq,
                                query_start,
                                target_start,
                            )
                        }
                    } else {
                        // Fall back to NM if no CIGAR available
                        alignment.nm as usize
                    };

                    if (alignment.nm as usize) < adjusted_errors {
                        log::trace!("Adjusted errors ({}) greater than NM ({}) for alignment between consensus {} and {}, using NM as adjusted errors",
                            adjusted_errors, alignment.nm, query_idx, target_idx);
                        adjusted_errors = alignment.nm as usize;
                    }

                    if adjusted_errors < 2{
                        log::trace!("Consensus {} (depth {}) maps to consensus {} (depth {}) with adjusted errors {}",
                            consensuses[query_idx].id, consensuses[query_idx].depth, consensuses[target_idx].id, consensuses[target_idx].depth, adjusted_errors);
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

                    // Handle special case of identical consensuses
                    if query_depth == target_depth{
                        if query_idx > *t_idx{
                            // To avoid circular merges, only allow one direction for identical consensuses
                            return true;
                        }
                        else{
                            return false;
                        }
                    }
                }

                log::trace!("Considering merge: Query {} (depth {}) -> Target {} (depth {}), Adjusted errors {}, Relative depth {:.4}, Threshold {:.4} => {}",
                    consensuses[query_idx].id, query_depth, consensuses[*t_idx].id, target_depth, nm, relative_depth, threshold,
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
        log::trace!(
            "Merging consensus {} (depth {}) into consensus {} (depth {})",
            consensuses[query_idx].id,
            consensuses[query_idx].depth,
            consensuses[target_idx].id,
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
            let mut new_cons = ConsensusSequence::new(consensus.sequence, consensus.hp_lengths, new_depth, consensus.id, new_cluster);
            new_cons.decompress();
            new_consensuses.push(new_cons);
        }
    }

    log::info!(
        "Consensus merging complete: {} -> {} consensuses",
        new_clusters.len(),
        new_consensuses.len()
    );


    let final_file = temp_dir.join("final_clusters_merged.tsv");
    write_clusters_tsv(&new_consensuses, twin_reads, &final_file, "final")
        .expect("Failed to write final_clusters_merged.tsv");

    // Write merged consensus sequences (before chimera filtering)
    let merged_fasta = temp_dir.join("merged_consensus_sequences.fasta");
    write_consensus_fasta(&new_consensuses, &merged_fasta, "merged")
        .expect("Failed to write merged_consensus_sequences.fasta");
    log::info!("Wrote {} merged consensus sequences to merged_consensus_sequences.fasta", new_consensuses.len());

    new_consensuses
}

/// Equivalence class: a set of ASVs that a group of reads maps to with equal quality
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct EquivalenceClass {
    asv_indices: Vec<usize>,
}

/// Refine ASV depths using EM algorithm on read-level mappings
/// Maps all reads back to final ASVs using k-mer/SNPmer-based approach
pub fn refine_asv_depths_with_em(
    twin_reads: &[TwinRead],
    consensuses: &mut Vec<ConsensusSequence>,
    kmer_info: &KmerGlobalInfo,
    args: &Cli,
    temp_dir: &PathBuf,
) {
    if consensuses.is_empty() {
        log::warn!("No consensuses to refine");
        return;
    }

    let mapping_file_writer = Mutex::new(std::io::BufWriter::new(
        std::fs::File::create(temp_dir.join("read_to_asv_mappings.tsv"))
            .expect("Failed to create read_to_asv_mappings.tsv"),
    ));

    log::info!("Refining ASV depths using EM algorithm with k-mer based read mapping");

    // Step 1: Load ASV sequences as TwinReads with SNPmers and minimizers
    let asv_fasta_path = temp_dir.join("final_asvs_for_em.fasta");
    write_consensus_fasta(&consensuses, &asv_fasta_path, "em_refinement")
        .expect("Failed to write ASVs for EM refinement");

    let asv_twin_reads = kmer_comp::twin_reads_from_fasta(&asv_fasta_path, kmer_info, args.kmer_size, args.c, args.blockmer_length, args.minimum_base_quality);
    log::info!("Loaded {} ASVs as TwinReads for k-mer comparison", asv_twin_reads.len());

    // Step 2: Build SNPmer index for all ASVs
    let k = args.kmer_size;
    let mask = !(3 << (k - 1)); // Mask for creating splitmers
    let mut asv_snpmer_index: FxHashMap<u64, Vec<(usize, Kmer48)>> = FxHashMap::default();

    for (asv_idx, asv_read) in asv_twin_reads.iter().enumerate() {
        let asv_snpmers = asv_read.snpmers_vec();
        for (_pos, kmer) in asv_snpmers {
            let splitmer = kmer.to_u64() & mask;
            asv_snpmer_index.entry(splitmer).or_insert_with(Vec::new).push((asv_idx, kmer));
        }
    }

    log::info!("Built SNPmer index with {} unique splitmers", asv_snpmer_index.len());
    log::info!("Mapping {} reads to {} ASVs using k-mer comparison", twin_reads.len(), consensuses.len());

    // Step 3: Map all reads to ASVs using k-mer comparison and collect equivalence classes
    let eq_classes = Mutex::new(HashMap::new());
    let filtered_reads_count = Mutex::new(0usize);
    let total_assigned_reads = Mutex::new(0usize);

    // Step 4: populate consensus seqs
    for cons in consensuses.iter_mut() {
        cons.unambig_best_read_map_count = Some(0);
        cons.ambig_read_map_count = Some(0);
        cons.num_map_leq_10nm = Some(0);
    }

    let unambig_read_map_count = Mutex::new(vec![0usize; consensuses.len()]);
    let ambig_read_map_count = Mutex::new(vec![0usize; consensuses.len()]);
    let num_map_leq_10nm = Mutex::new(vec![0usize; consensuses.len()]);

    twin_reads.par_iter().enumerate().for_each(|(_read_idx, twin_read)| {
        let read_snpmers = twin_read.snpmer_kmers();
        let read_minimizers: FxHashSet<Kmer48> = twin_read.minimizer_kmers().into_iter().cloned().collect();

        // Compare read against each ASV using SNPmers
        let candidate_stats = find_compatible_candidates(&asv_snpmer_index, &read_snpmers, k);

        // For each candidate ASV, calculate ratio: mismatched_snpmers / matched_minimizers
        let mut asv_scores: Vec<(usize, f64, usize, usize)> = Vec::new(); // (asv_idx, ratio, mismatches, minimizer_matches)

        for (asv_idx, (_matches, mismatches)) in candidate_stats {
            // Count matching minimizers between read and this ASV
            let asv_minimizers: FxHashSet<Kmer48> = asv_twin_reads[asv_idx].minimizer_kmers().into_iter().cloned().collect();
            let minimizer_matches = read_minimizers.intersection(&asv_minimizers).count();

            if minimizer_matches == 0 {
                continue; // Skip if no minimizer overlap
            }

            if (minimizer_matches as f64 / read_minimizers.len().min(asv_minimizers.len()) as f64)
                < (0.950f64).powi(k as i32) {
                continue; // Skip if less than 10% minimizer overlap
            }

            // Calculate ratio: mismatched_snpmers / matched_minimizers
            let ratio = mismatches as f64 / minimizer_matches as f64 / args.c as f64;

            asv_scores.push((asv_idx, ratio, mismatches, minimizer_matches));
        }

        // Find minimum ratio (best matches)
        if asv_scores.is_empty() {
            *filtered_reads_count.lock().unwrap() += 1;
            return;
        }

        let min_ratio = asv_scores.iter().min_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap();

        let min_mismatches = min_ratio.2;
        let max_mini = min_ratio.3;
        let min_ratio = min_ratio.1;

        // Filter by ratio threshold (0.005) and keep all ASVs with minimum ratio
        let threshold = 0.0050;
        let mut best_asv_indices: Vec<(usize, usize)> = asv_scores.iter()
            .filter(|(_, ratio, _, _)| *ratio <= threshold)
            .map(|(asv_idx, _, mismatches, _)| (*asv_idx, *mismatches))
            .collect();

        if best_asv_indices.is_empty() {
            *filtered_reads_count.lock().unwrap() += 1;
            return;
        }

        // Among best ASVs, keep those with lowest mismatches
        best_asv_indices.sort_by(|a, b| a.1.cmp(&b.1)); // Sort by mismatches
        let lowest_mismatches = best_asv_indices[0].1;
        best_asv_indices.retain(|(_, mismatches)| *mismatches == lowest_mismatches);

        // For each best ASV, count number of kmer matches for tie-breaking
        let mut best_alns: Vec<(usize, i32, usize)> = Vec::new(); // (asv_idx, nm, mismatches)
        let seq_u8 : Vec<u8> = twin_read.dna_seq.iter().map(|x| x.to_char().to_ascii_uppercase() as u8).collect();
        let aligner = Aligner::builder()
            .lrhq()
            .with_index_threads(1) // Use 1 thread per aligner since we parallelize over consensuses
            .with_cigar()
            .with_seq(&seq_u8)
            .expect("Failed to create aligner");
                
        for (asv_idx, mismatches) in best_asv_indices.into_iter() {
            let asv_tr = &asv_twin_reads[asv_idx];
            let seq_u8_asv : Vec<u8> = asv_tr.dna_seq.iter().map(|x| x.to_char().to_ascii_uppercase() as u8).collect();
            let alignment_result = aligner.map(&seq_u8_asv, true, false, None, None, None).unwrap();
            if alignment_result.is_empty() {
                continue;
            }
            best_alns.push((asv_idx, alignment_result[0].alignment.as_ref().unwrap().nm, mismatches));
        }

        // Sort best alignments by number of mismatches (ascending)
        best_alns.sort_by(|a, b| a.1.cmp(&b.1));

        let best_nm = best_alns.first().map(|x| x.1).unwrap_or(i32::MAX);
        let mut best_asv_aln_indices: Vec<usize> = best_alns.iter()
            .filter(|(_, nm, _)| *nm == best_nm)
            .map(|(asv_idx, _, _)| *asv_idx)
            .collect::<Vec<usize>>();

        // Map to file
        {
            let mut mapping_file_writer = mapping_file_writer.lock().unwrap();
            for (asv_idx, mini_matches, mismatches) in best_alns.iter().take(5) {
                writeln!(
                    mapping_file_writer,
                    "{}\tasv:{}\t{}\t{}",
                    twin_read.id,
                    consensuses[*asv_idx].id,
                   // asv_idx,
                    mismatches,
                    mini_matches
                ).expect("Failed to write to read_to_asv_mappings.tsv");
            }
        }

        // Only keep reads that have at least one good mapping
        if !best_asv_aln_indices.is_empty() {
            // Sort to ensure consistent equivalence class keys
            best_asv_aln_indices.sort();

            let eq_class = EquivalenceClass {
                asv_indices:best_asv_aln_indices,
            };

            if eq_class.asv_indices.len() == 1{
                let asv_idx = eq_class.asv_indices[0];
                let mut unambig_counts = unambig_read_map_count.lock().unwrap();
                unambig_counts[asv_idx] += 1;
            }
            else{
                for &asv_idx in eq_class.asv_indices.iter(){
                    let mut ambig_counts = ambig_read_map_count.lock().unwrap();
                    ambig_counts[asv_idx] += 1;
                }
            }

            if best_nm <= 10 {
                for &asv_idx in eq_class.asv_indices.iter(){
                    let mut leq10nm_counts = num_map_leq_10nm.lock().unwrap();
                    leq10nm_counts[asv_idx] += 1;
                }
            }

            // Add to equivalence class with count
            let mut eq_map = eq_classes.lock().unwrap();
            *eq_map.entry(eq_class).or_insert(0) += 1;
            *total_assigned_reads.lock().unwrap() += 1;
        } else {
            log::trace!("Read filtered out due to high ratio mappings: ratio {:.4} mismatch {} mini {}", min_ratio, min_mismatches, max_mini);
            *filtered_reads_count.lock().unwrap() += 1;
        }
    });

    let eq_classes = eq_classes.into_inner().unwrap();
    let filtered_reads = filtered_reads_count.into_inner().unwrap();
    let total_assigned = total_assigned_reads.into_inner().unwrap();

    log::info!("Filtered {} reads with ratio > 0.005", filtered_reads);
    log::info!("Total assigned reads: {}", total_assigned);
    log::info!("Total percentage assigned: {:.2}%", (total_assigned as f64 / (total_assigned + filtered_reads) as f64) * 100.0);
    log::info!("Number of unique equivalence classes: {}", eq_classes.len());

    let unambig_read_map_count = unambig_read_map_count.into_inner().unwrap();
    let ambig_read_map_count = ambig_read_map_count.into_inner().unwrap();
    let num_map_leq_10nm = num_map_leq_10nm.into_inner().unwrap();

    for i in 0..consensuses.len() {
        consensuses[i].unambig_best_read_map_count = Some(unambig_read_map_count[i]);
        consensuses[i].ambig_read_map_count = Some(ambig_read_map_count[i]);
        consensuses[i].num_map_leq_10nm = Some(num_map_leq_10nm[i]);
    }

    //debug equiv classes

    for (eq_class, count) in &eq_classes {
        log::trace!("Equivalence class: ASVs {:?}, Count {}", eq_class.asv_indices, count);
    }

    if eq_classes.is_empty() {
        log::warn!("No reads mapped well to ASVs. Keeping original depths.");
        return;
    }

    // Step 4: Run EM algorithm on equivalence classes
    let num_asvs = consensuses.len();
    let mut asv_abundances = vec![1.0 / num_asvs as f64; num_asvs];
    let convergence_threshold = 0.01 / total_assigned as f64;

    log::info!("Running EM algorithm with convergence threshold: {:.6e}", convergence_threshold);

    let mut iteration = 0;
    const MAX_ITERATIONS: usize = 10000;

    loop {
        iteration += 1;
        let mut new_asv_abundances = vec![0.0; num_asvs];

        // E-step + M-step: distribute reads proportionally based on current abundances
        for (eq_class, count) in &eq_classes {
            let denominator: f64 = eq_class.asv_indices.iter()
                .map(|&asv_idx| asv_abundances[asv_idx])
                .sum();

            if denominator > 0.0 {
                for &asv_idx in &eq_class.asv_indices {
                    let contribution = (*count as f64) * asv_abundances[asv_idx] / denominator;
                    new_asv_abundances[asv_idx] += contribution;
                }
            }
        }

        // Normalize by total assigned reads
        let total_counts: f64 = new_asv_abundances.iter().sum();
        if total_counts > 0.0 {
            for abundance in new_asv_abundances.iter_mut() {
                *abundance /= total_assigned as f64;
            }
        }

        // Check convergence
        let max_change = asv_abundances.iter()
            .zip(new_asv_abundances.iter())
            .map(|(old, new)| (old - new).abs())
            .fold(0.0, f64::max);

        asv_abundances = new_asv_abundances;

        if max_change < convergence_threshold || iteration >= MAX_ITERATIONS {
            log::info!("EM converged after {} iterations (max change: {:.6e})", iteration, max_change);
            break;
        }

        if iteration % 10 == 0 {
            log::debug!("EM iteration {}: max change {:.6e}", iteration, max_change);
        }
    }

    // Step 5: Update consensus depths based on EM abundances
    log::info!("Updating ASV depths based on EM abundances");
    for (asv_idx, consensus) in consensuses.iter_mut().enumerate() {
        let em_abundance = asv_abundances[asv_idx];
        let new_depth = (em_abundance * total_assigned as f64).round() as usize;
        if new_depth != 0 && consensus.depth / new_depth > 10 {
            log::debug!("ASV {} possible removal due to coverage drop of 10x: original depth {}, new depth {}",
                asv_idx, consensus.depth, new_depth);
            consensus.depth = new_depth;
            // consensus.depth = 0; // Force removal
        }
        else{
            log::debug!("ASV {}: Original depth = {}, EM abundance = {:.6}, New depth = {}",
                asv_idx, consensus.depth, em_abundance, new_depth);
            consensus.depth = new_depth;
        }
    }

    // Filter out ASVs with zero depth after EM
    let original_count = consensuses.len();
    consensuses.retain(|c| c.depth > 0);
    let filtered_count = original_count - consensuses.len();

    if filtered_count > 0 {
        log::info!("Filtered {} ASVs with zero depth after EM refinement", filtered_count);
    }

    log::info!("EM refinement complete: {} ASVs remaining", consensuses.len());
}




// Check for within-ASV heterogeneity by analyzing 50bp blocks of aligned reads
// Main entry point that processes all ASVs in parallel
// pub fn check_asv_heterogeneity(
//     twin_reads: &[TwinRead],
//     consensuses: &mut Vec<ConsensusSequence>,
//     _args: &Cli,
// ) {
//     let max_reads_per_asv = 500;
//     let block_size = 50;
//     let min_count_fraction = 0.15; // 15% of max count
//     let min_edit_distance = 1;
//     let edge_buffer = 1; // Skip first and last block
//     let min_haplotype_fraction = 0.10; // 10% minimum abundance for new variants

//     log::info!("Checking for within-ASV heterogeneity using 50bp blocks");

//     // Mutex for collecting new consensus sequences from parallel iterations
//     let new_consensuses_mutex: Mutex<Vec<ConsensusSequence>> = Mutex::new(Vec::new());

//     // Process each consensus in parallel
//     consensuses.par_iter().enumerate().for_each(|(asv_idx, consensus)| {
//         let new_variants = check_single_asv_heterogeneity(
//             asv_idx,
//             consensus,
//             twin_reads,
//             max_reads_per_asv,
//             block_size,
//             min_count_fraction,
//             min_edit_distance,
//             edge_buffer,
//             min_haplotype_fraction,
//         );

//         // Collect any new consensuses
//         if !new_variants.is_empty() {
//             let mut new_consensuses = new_consensuses_mutex.lock().unwrap();
//             new_consensuses.extend(new_variants);
//         }
//     });

//     // Append new variant consensuses to the main list
//     let new_consensuses = new_consensuses_mutex.into_inner().unwrap();
//     let num_new_variants = new_consensuses.len();

//     if num_new_variants > 0 {
//         log::info!("Generated {} new low-frequency variant consensuses from heterogeneity analysis", num_new_variants);
//         consensuses.extend(new_consensuses);
//     }

//     log::info!("Heterogeneity check complete: {} total ASVs", consensuses.len());
// }
