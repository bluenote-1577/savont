use crate::asv_cluster::find_compatible_candidates;
use crate::cli::ClusterArgs as Cli;
use crate::constants::*;
use crate::kmer_comp;
use crate::{seeding, types::*, utils};
use bio_seq::prelude::*;
use fxhash::{FxHashMap, FxHashSet};
use minimap2::Aligner;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Mutex;

type BaseObservation = (u8, u8, u8);

/// A decoded event. Only the packed word is stored at the reference position.
#[derive(Debug, PartialEq, Eq)]
enum PileupBase<'a> {
    /// A matched/mismatched base with quality score and homopolymer length
    /// (base, quality, hp_length from this read)
    Base(u8, u8, u8),
    /// A deletion at this reference position
    Deletion,
    /// The inserted bases, qualities and homopolymer lengths for one CIGAR event.
    Insertion(&'a [BaseObservation]),
}

#[derive(Default)]
struct Insertions {
    bases: Vec<BaseObservation>,
    ends: Vec<usize>,
}

impl Insertions {
    fn push(&mut self, bases: impl IntoIterator<Item = BaseObservation>) -> usize {
        let id = self.ends.len();
        self.bases.extend(bases);
        self.ends.push(self.bases.len());
        id
    }

    fn get(&self, id: usize) -> &[BaseObservation] {
        let start = if id == 0 { 0 } else { self.ends[id - 1] };
        &self.bases[start..self.ends[id]]
    }
}

/// All positions and their shared insertion storage for one consensus.
pub struct ClusterPileup {
    positions: Vec<Pileup>,
    insertions: Insertions,
}

/// Pileup information at a single reference position.
struct Pileup {
    ref_pos: usize,
    ref_base: u8,
    ref_hp_length: u8,
    events: Vec<u32>,
    wide_insertions: usize,
    alt_posterior: Option<f64>,
}

impl Pileup {
    const DELETION: u32 = 1 << 24;
    const INSERTION: u32 = 1 << 31;
    const WIDE_INSERTION: u32 = u32::MAX;

    /// Preallocate one event per aligned read to avoid repeated growth at each position.
    fn with_capacity(ref_pos: usize, ref_base: u8, ref_hp_length: u8, cap: usize) -> Self {
        Self {
            ref_pos,
            ref_base,
            ref_hp_length,
            events: Vec::with_capacity(cap),
            wide_insertions: 0,
            alt_posterior: None,
        }
    }

    fn add_base(&mut self, base: u8, quality: u8, hp_length: u8) {
        self.events
            .push(u32::from_le_bytes([base, quality, hp_length, 0]));
    }

    fn add_deletion(&mut self) {
        self.events.push(Self::DELETION);
    }

    fn add_insertion(&mut self, id: usize) {
        if id < (Self::WIDE_INSERTION ^ Self::INSERTION) as usize {
            self.events.push(Self::INSERTION | id as u32);
        } else {
            // Escape to a full-width index without limiting the number of events.
            let id = id as u64;
            self.events
                .extend([Self::WIDE_INSERTION, id as u32, (id >> 32) as u32]);
            self.wide_insertions += 1;
        }
    }

    fn observations<'a>(&'a self, insertions: &'a Insertions) -> impl Iterator<Item = PileupBase<'a>> {
        self.decoded_events().map(|(word, id)| match id {
            Some(id) => PileupBase::Insertion(insertions.get(id)),
            None if word == Self::DELETION => PileupBase::Deletion,
            None => {
                let [base, quality, hp_length, _] = word.to_le_bytes();
                PileupBase::Base(base, quality, hp_length)
            }
        })
    }

    fn decoded_events(&self) -> impl Iterator<Item = (u32, Option<usize>)> + '_ {
        let mut events = self.events.iter().copied();
        std::iter::from_fn(move || {
            let word = events.next()?;
            let id = if word == Self::WIDE_INSERTION {
                let low = events.next().unwrap() as u64;
                let high = events.next().unwrap() as u64;
                Some((low | (high << 32)) as usize)
            } else if word & Self::INSERTION != 0 {
                Some((word ^ Self::INSERTION) as usize)
            } else {
                None
            };
            Some((word, id))
        })
    }

    fn depth(&self) -> usize {
        self.events.len() - 2 * self.wide_insertions
    }

    fn depth_nodeletion(&self) -> (usize, usize, usize) {
        let mut insertions = 0;
        let mut deletions = 0;
        for (word, id) in self.decoded_events() {
            insertions += usize::from(id.is_some());
            deletions += usize::from(word == Self::DELETION);
        }
        (self.depth() - deletions, insertions, deletions)
    }
}

#[cfg(test)]
mod pileup_tests {
    use super::*;

    #[test]
    fn packed_observations_preserve_insertion_payloads_and_order() {
        let mut insertions = Insertions::default();
        let mut first = Pileup::with_capacity(10, b'A', 1, 4);
        let mut second = Pileup::with_capacity(11, b'C', 1, 1);
        let inserted_base = (b'G', 75, 255);
        let long_insertion = vec![inserted_base; 500];
        let short_insertion = [(b'T', 60, 3), (b'A', 69, 17)];

        first.add_base(b'N', 255, 255);
        first.add_insertion(insertions.push(long_insertion.iter().copied()));
        second.add_insertion(insertions.push([]));
        first.add_deletion();
        first.add_insertion(insertions.push(short_insertion));

        assert_eq!(first.events.len(), 4);
        assert_eq!(first.depth_nodeletion(), (3, 2, 1));
        assert_eq!(
            first.observations(&insertions).collect::<Vec<_>>(),
            vec![
                PileupBase::Base(b'N', 255, 255),
                PileupBase::Insertion(&long_insertion),
                PileupBase::Deletion,
                PileupBase::Insertion(&short_insertion),
            ]
        );
        assert_eq!(
            second.observations(&insertions).next(),
            Some(PileupBase::Insertion(&[]))
        );
    }

    #[test]
    fn wide_insertions_count_as_single_events() {
        let mut pileup = Pileup::with_capacity(0, b'A', 1, 4);
        let last_compact = (Pileup::INSERTION - 2) as usize;
        let first_wide = last_compact + 1;
        let high_id = usize::try_from((1_u64 << 32) + 11).unwrap_or(usize::MAX);
        pileup.add_insertion(last_compact);
        pileup.add_insertion(first_wide);
        pileup.add_insertion(high_id);
        pileup.add_deletion();

        assert_eq!(pileup.depth(), 4);
        assert_eq!(pileup.depth_nodeletion(), (3, 3, 1));
        assert_eq!(
            pileup.decoded_events().map(|(_, id)| id).collect::<Vec<_>>(),
            vec![Some(last_compact), Some(first_wide), Some(high_id), None]
        );
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
                        if query_seq[query_pos] != target_seq[target_pos]
                            && (query_seq[query_pos] != b'N' && target_seq[target_pos] != b'N')
                        {
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
                        if len < 10 {
                            error_count += 1;
                        } else {
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
                        if len < 10 {
                            error_count += 1;
                        } else {
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
    _cluster_idx: usize,
    no_band: bool,
) -> Vec<u8> {
    use spoars::align::{AlignmentEngine, AlignmentType, BandConfig, Scoring, SimdEngine};
    use spoars::graph::Graph;

    if sequences.is_empty() {
        return Vec::new();
    }

    let scoring = Scoring::new(3, -8, -6, -6, 0, 0).expect("invalid SPOA scoring");

    // Band = 0.1 * seq_len + |ref_len - seq_len|.  Use mean length as ref_len and the
    // cluster's max deviation as the fixed base so one engine covers all sequences.
    let ref_len = sequences.iter().map(|s| s.len()).sum::<usize>() / sequences.len();
    let max_deviation = sequences
        .iter()
        .map(|s| (ref_len as isize - s.len() as isize).abs() as u32)
        .max()
        .unwrap_or(0);

    let mut engine = if no_band {
        SimdEngine::new(AlignmentType::Overlap, scoring)
    } else {
        SimdEngine::banded(
            AlignmentType::Overlap,
            scoring,
            BandConfig {
                base: max_deviation,
                frac: 0.1,
            },
        )
    };
    let mut graph = Graph::new();

    for i in 0..sequences.len() {
        let u32_qual: Vec<u32> = qualities[i].iter().map(|&q| q as u32).collect();
        let (alignment, _score) = engine.align(&sequences[i], &graph);
        graph
            .add_alignment(&alignment, &sequences[i], &u32_qual)
            .expect("add_alignment failed");
    }

    graph.generate_consensus().into_bytes()
}

pub fn align_and_consensus(
    twin_reads: &[TwinRead],
    clusters: Vec<Vec<usize>>,
    args: &Cli,
    output_dir: &PathBuf,
) -> Vec<ConsensusSequence> {
    let max_seqs_consensus = 75;
    let accuracy: [f64; 16] = std::array::from_fn(|i| {
        let q = i as u8 * 3 + 33;
        1.0 - 10.0f64.powf(-((q - 33) as f64) / 10.0)
    });

    // Log which POA implementation is being used
    log::info!("Generating consensus sequences from SNPmer clusters...");

    let consensus_seqs = Mutex::new(Vec::new());
    // Implementation of alignment and consensus generation
    clusters.par_iter().enumerate().for_each(|(cluster_idx, cluster)| {
        // Compute the same ranking keys from the compact reads, and expand only the
        // selected seed and reads used for POA.
        let mut lengths_and_i: Vec<_> = cluster.iter().enumerate()
            .map(|(i, &read_id)| (twin_reads[read_id].dna_seq.len(), i)).collect();
        lengths_and_i.sort_by_key(|k| k.0);
        let mut avg_qual_and_i: Vec<_> = cluster.iter().enumerate()
            .map(|(i, &read_id)| {
                let read = &twin_reads[read_id];
                let q = if let Some(quality) = &read.qual_seq {
                    if quality.is_empty() { 1.0 } else {
                        quality.iter().map(|x| accuracy[x as usize]).sum::<f64>() / quality.len() as f64
                    }
                } else if read.dna_seq.is_empty() { 1.0 } else { 0.0 };
                (q, i)
            }).collect();
        avg_qual_and_i.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        let largest_sequence_index = lengths_and_i[(lengths_and_i.len() as f64 * 0.9) as usize].1;
        let range = avg_qual_and_i[0..max_seqs_consensus.min(avg_qual_and_i.len())].iter().map(|(_, i)| *i);

        let expand = |i: usize| {
            let read = &twin_reads[cluster[i]];
            let seq: Vec<u8> = read.dna_seq.iter()
                .map(|x| x.to_char().to_ascii_uppercase() as u8).collect();
            let mut qualities: Vec<u8> = if let Some(quality) = &read.qual_seq {
                quality.iter().flat_map(|x| std::iter::repeat_n((x as u8) * 3 + 33, QUALITY_SEQ_BIN)).collect()
            } else { vec![33; seq.len()] };
            if qualities.len() > seq.len() {
                qualities.truncate(seq.len());
            } else if qualities.len() < seq.len() {
                let last = *qualities.last().unwrap();
                qualities.resize(seq.len(), last);
            }
            (seq, qualities)
        };
        let (seed, seed_quality) = expand(largest_sequence_index);
        let aligner = Aligner::builder().map_ont().with_index_threads(args.threads).with_cigar().with_seq(&seed).expect("Failed to create aligner");
        let mut mappings = Vec::new();
        for i in range {
            if i == largest_sequence_index { continue; }
            let (seq, qual) = expand(i);
            let Ok(alignment) = aligner.map(&seq, true, false, None, None, None) else {
                log::debug!("No alignment found for read {} in cluster {}", i, cluster_idx);
                continue;
            };
            mappings.push((i, alignment, seq, qual));
        }
        mappings.sort_by_key(|k| k.0);
        let mut aligned_sequences = vec![seed];
        let mut aligned_qualities = vec![seed_quality];
        for (i, mappings, mut seq, mut qual) in mappings {
            let Some(best_mapping) = mappings.first() else {
                log::debug!("No alignment found for read {} in cluster {}", i, cluster_idx);
                continue;
            };
            let cigar_str = &best_mapping.alignment.as_ref().unwrap().cigar_str.as_ref().unwrap();
            log::trace!("Read {} len {}: CIGAR: {}, Query start: {}, Query end: {}, Target start: {}, Target end: {}, Cluster IDX: {}",
                i, best_mapping.query_len.unwrap(), cigar_str, best_mapping.query_start, best_mapping.query_end, best_mapping.target_start, best_mapping.target_end, cluster_idx);
            if best_mapping.strand == minimap2::Strand::Reverse {
                seq = utils::reverse_complement(&seq);
                qual.reverse();
            }
            aligned_sequences.push(seq);
            aligned_qualities.push(qual);
            if aligned_sequences.len() > max_seqs_consensus { break; }
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
        let hpc_consensus = generate_consensus_poa(&hpc_aligned_sequences, &hpc_aligned_qualities, coverage_threshold, cluster_idx, args.no_band);

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

        log::trace!("Completed alignment for cluster of size {}", cluster.len());
    });

    let mut consensus_seqs = consensus_seqs.into_inner().unwrap();
    consensus_seqs.sort_by(|a, b| b.3.cmp(&a.3).then_with(|| a.0.cmp(&b.0)));
    let consensus_seqs: Vec<ConsensusSequence> = consensus_seqs
        .into_iter()
        .map(|(id, seq, hp_lens, depth, cluster)| {
            ConsensusSequence::new(seq, hp_lens, depth, id, cluster)
        })
        .collect();

    // Write HPC consensus sequences to file
    let consensus_path = output_dir.join("consensus_sequences.fasta");
    write_consensus_fasta(&consensus_seqs, &consensus_path, "initial")
        .expect("Failed to write consensus_sequences.fasta");
    log::info!(
        "Wrote {} consensus sequences to consensus_sequences.fasta",
        consensus_seqs.len()
    );

    consensus_seqs
}

/// Keep the original conservative batch weight so packing reduces memory use per batch.
fn pileup_batch_weight(consensus: &ConsensusSequence) -> usize {
    let positions = consensus.sequence.len();
    let depth = consensus.cluster.len().min(MAX_SEQS_CONSENSUS);
    positions * depth * 32 + positions * 56
}

/// Groups all cluster indices into batches whose combined pileups stay within a memory
/// budget, so peak usage does not grow with the number of clusters.
///
/// Each batch is still processed in parallel across all threads, so this costs no
/// parallelism as long as a batch holds at least as many clusters as there are threads;
/// that floor is enforced even if it means exceeding the budget (a single cluster's
/// pileup is unavoidable, and the floor keeps all threads busy).
pub fn plan_pileup_batches(consensuses: &[ConsensusSequence], args: &Cli) -> Vec<Vec<usize>> {
    let all: Vec<usize> = (0..consensuses.len()).collect();
    plan_pileup_batches_for(consensuses, &all, args)
}

/// Same as `plan_pileup_batches`, for a subset of clusters (kept in the given order).
pub fn plan_pileup_batches_for(
    consensuses: &[ConsensusSequence],
    indices: &[usize],
    args: &Cli,
) -> Vec<Vec<usize>> {
    const PILEUP_BUDGET_BYTES: usize = 2 * 1024 * 1024 * 1024;
    let min_batch = args.threads.max(1);

    let mut batches: Vec<Vec<usize>> = Vec::new();
    let mut current: Vec<usize> = Vec::new();
    let mut current_bytes = 0usize;

    for &idx in indices {
        let cost = pileup_batch_weight(&consensuses[idx]);
        if !current.is_empty()
            && current.len() >= min_batch
            && current_bytes + cost > PILEUP_BUDGET_BYTES
        {
            batches.push(std::mem::take(&mut current));
            current_bytes = 0;
        }
        current.push(idx);
        current_bytes += cost;
    }
    if !current.is_empty() {
        batches.push(current);
    }
    batches
}

/// Builds pileups for the clusters named in `cluster_indices` (global indices into
/// `consensuses`). The returned Vec is ordered to match `cluster_indices`, not
/// `consensuses`. Callers process clusters in batches to limit peak memory.
pub fn generate_consensus_pileups(
    twin_reads: &[TwinRead],
    consensuses: &mut [ConsensusSequence],
    args: &Cli,
    cluster_indices: &[usize],
) -> Vec<ClusterPileup> {
    let max_seqs_consensus = MAX_SEQS_CONSENSUS;

    let pileups = Mutex::new(Vec::new());

    // Immutable reborrow so the parallel section can read while `consensuses` stays
    // `&mut` for the hp_length writeback below.
    let consensuses_ref: &[ConsensusSequence] = consensuses;

    // Process the requested consensuses and their reads in parallel.
    // `local_idx` indexes `cluster_indices`; `cluster_idx` is the global consensus index.
    cluster_indices
        .par_iter()
        .enumerate()
        .for_each(|(local_idx, &cluster_idx)| {
            let consensus = &consensuses_ref[cluster_idx];
            let cluster = &consensus.cluster;
            let consensus_seq = &consensus.sequence;
            let consensus_hp_lengths = &consensus.hp_lengths;

            // Align reads from this cluster back to consensus
            let reads_to_align = cluster.len().min(max_seqs_consensus);

            // Initialize pileup for this consensus with placeholder ref_hp_length (will be updated later)
            let mut cluster_pileup: Vec<Pileup> = consensus_seq
                .iter()
                .enumerate()
                .map(|(pos, &base)| {
                    Pileup::with_capacity(pos, base, consensus_hp_lengths[pos], reads_to_align)
                })
                .collect();
            let mut insertions = Insertions::default();

            // Create aligner with this consensus as reference
            let aligner = Aligner::builder()
                .map_ont()
                .with_cigar()
                .with_seq(consensus_seq)
                .expect("Failed to create aligner");

            for i in 0..reads_to_align {
                let read_idx = cluster[i];
                let twin_read = &twin_reads[read_idx];

                // Get sequence and quality
                let seq_u8: Vec<u8> = twin_read
                    .dna_seq
                    .iter()
                    .map(|x| x.to_char().to_ascii_uppercase() as u8)
                    .collect();

                let qual_u8 = if let Some(qual_seq) = &twin_read.qual_seq {
                    qual_seq.iter().map(|x| (x as u8) * 3 + 33).collect()
                } else {
                    vec![33; twin_read.dna_seq.len()]
                };

                // Bin qualities similar to align_and_consensus
                let bin_size = QUALITY_SEQ_BIN;
                let mut query_quals_u8: Vec<u8> =
                    qual_u8.iter().flat_map(|x| vec![*x; bin_size]).collect();

                // Adjust quality length to match sequence length
                if query_quals_u8.len() > seq_u8.len() {
                    query_quals_u8.truncate(seq_u8.len());
                } else if query_quals_u8.len() < seq_u8.len() {
                    let last_qual = query_quals_u8[query_quals_u8.len() - 1];
                    query_quals_u8.extend(vec![last_qual; seq_u8.len() - query_quals_u8.len()]);
                }

                // HPC compress the read before aligning to HPC consensus
                let (hpc_seq, hpc_qual, hp_lens) = utils::homopolymer_compress_with_quality(
                    &seq_u8,
                    &query_quals_u8,
                    args.use_hpc,
                );

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
                                } else {
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
                                                if ref_pos + j < cluster_pileup.len()
                                                    && query_pos + j < mapped_seq.len()
                                                {
                                                    let base = mapped_seq[query_pos + j];
                                                    let qual = mapped_qual[query_pos + j];
                                                    let hp_len = mapped_hp_lens[query_pos + j];
                                                    cluster_pileup[ref_pos + j]
                                                        .add_base(base, qual, hp_len);
                                                }
                                            }
                                            ref_pos += len;
                                            query_pos += len;
                                        }
                                        1 => {
                                            // Insertion in read - associate with previous ref position
                                            if ref_pos > 0
                                                && ref_pos - 1 < cluster_pileup.len()
                                                && query_pos + len <= mapped_seq.len()
                                            {
                                                let end = query_pos + len.min(MAX_INSERTION_LENGTH);
                                                let id = insertions.push((query_pos..end).map(|j| {
                                                    (mapped_seq[j], mapped_qual[j], mapped_hp_lens[j])
                                                }));
                                                cluster_pileup[ref_pos - 1].add_insertion(id);
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
                                            log::warn!(
                                                "Unexpected CIGAR operation in pileup: {}",
                                                op as char
                                            );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            log::trace!(
                "Generated pileup for consensus {} with {} positions",
                cluster_idx,
                cluster_pileup.len()
            );
            pileups.lock().unwrap().push((
                local_idx,
                ClusterPileup {
                    positions: cluster_pileup,
                    insertions,
                },
            ));
        });

    let mut pileups = pileups.into_inner().unwrap();
    pileups.sort_by_key(|k| k.0);
    let mut pileups: Vec<ClusterPileup> = pileups.into_iter().map(|(_, pileup)| pileup).collect();

    // Calculate modal HP lengths from aligned reads for each position
    for cluster in &mut pileups {
        for pileup in &mut cluster.positions {
            // Collect all HP lengths from aligned bases at this position
            let mut hp_lengths: Vec<u8> = Vec::new();
            for base_entry in pileup.observations(&cluster.insertions) {
                if let PileupBase::Base(_base, _qual, hp_len) = base_entry {
                    hp_lengths.push(hp_len);
                }
            }

            // Calculate modal HP length
            if !hp_lengths.is_empty() {
                // Count occurrences of each HP length
                let mut counts: std::collections::HashMap<u8, usize> =
                    std::collections::HashMap::new();
                for &hp_len in &hp_lengths {
                    *counts.entry(hp_len).or_insert(0) += 1;
                }

                // Find the most common HP length
                //let modal_hp_length = counts.iter()
                //    .max_by_key(|(_, count)| *count)
                //    .map(|(hp_len, _)| *hp_len)
                //    .unwrap_or(1);
                let median_hp_length = {
                    let mut sorted_hp_lengths = hp_lengths.clone();
                    sorted_hp_lengths.sort_unstable();
                    let mid = sorted_hp_lengths.len() / 2;
                    if sorted_hp_lengths.len() % 2 == 0 {
                        ((sorted_hp_lengths[mid - 1] as u16 + sorted_hp_lengths[mid] as u16) / 2)
                            as u8
                    } else {
                        sorted_hp_lengths[mid]
                    }
                };

                pileup.ref_hp_length = median_hp_length;
            } else {
                // No bases aligned, keep placeholder value
                pileup.ref_hp_length = 1;
            }
        }
    }

    // debug print out

    if log::log_enabled!(log::Level::Trace) {
        for (i, cluster) in pileups.iter().enumerate() {
            log::trace!("Pileup for consensus {}:", i);
            for pos in &cluster.positions {
                let bases_str: Vec<String> = pos
                    .observations(&cluster.insertions)
                    .map(|b| match b {
                        PileupBase::Base(base, qual, hp_len) => {
                            format!("{}(q={},hp={})", base as char, qual, hp_len)
                        }
                        PileupBase::Deletion => String::from("D"),
                        PileupBase::Insertion(data) => {
                            let bases: Vec<u8> = data.iter().map(|(b, _, _)| *b).collect();
                            format!("I({})", String::from_utf8_lossy(&bases))
                        }
                    })
                    .collect();
                log::trace!(
                    "Pos {}: Ref base: {}, Ref HP len: {}, Depth: {}, Bases: {}",
                    pos.ref_pos,
                    pos.ref_base as char,
                    pos.ref_hp_length,
                    pos.depth(),
                    bases_str.join(", ")
                );
            }
        }
    }

    // Update consensus HP lengths from modal values calculated from pileups.
    // Indexed through cluster_indices since `pileups` only covers the requested subset.
    // This is idempotent: ref_hp_length is derived from the aligned reads (or set to 1
    // when nothing aligns), never from the incoming hp_lengths, so regenerating a
    // cluster's pileup in a later pass yields the same result.
    for (local_idx, &cluster_idx) in cluster_indices.iter().enumerate() {
        let modal_hp_lengths: Vec<u8> = pileups[local_idx]
            .positions
            .iter()
            .map(|p| p.ref_hp_length)
            .collect();
        consensuses[cluster_idx].hp_lengths = modal_hp_lengths;
    }

    pileups
}

/// Estimate error rate as a function of quality score from pileup data
/// Uses top N clusters and filters positions with <5% error rate
/// The top `top_frac` of clusters by depth, as global indices into `consensuses`.
/// Error-rate estimation only needs these, so the caller builds pileups for just
/// this subset rather than for every cluster.
pub fn top_cluster_indices(consensuses: &[ConsensusSequence], top_frac: f64) -> Vec<usize> {
    let mut cluster_depths: Vec<(usize, usize)> = consensuses
        .iter()
        .enumerate()
        .map(|(idx, cons)| (idx, cons.depth))
        .collect();
    cluster_depths.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    cluster_depths
        .iter()
        .take((top_frac * cluster_depths.len() as f64).round() as usize)
        .map(|(idx, _)| *idx)
        .collect()
}

/// Per-quality `(errors, total)` base counts, accumulated over pileups.
pub type QualityStats = HashMap<u8, (usize, usize)>;

/// Adds one batch of pileups to `quality_stats`. The counts are plain sums, so batches can
/// be accumulated in any order with the same result. `pileups` must already be restricted
/// to the clusters of interest (see `top_cluster_indices`).
pub fn accumulate_quality_stats(quality_stats: &mut QualityStats, pileups: &[ClusterPileup]) {
    let prior_count = 1;

    for cluster in pileups {
        for pileup in &cluster.positions {
            // Calculate error rate at this position
            let mut total_bases = 0;
            let mut error_bases = 0;

            for base_entry in pileup.observations(&cluster.insertions) {
                match base_entry {
                    PileupBase::Base(base, _qual, _hp_len) => {
                        total_bases += 1;
                        if base != pileup.ref_base {
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
                    for base_entry in pileup.observations(&cluster.insertions) {
                        if let PileupBase::Base(base, qual, _hp_len) = base_entry {
                            let entry = quality_stats
                                .entry(qual)
                                .or_insert((prior_count, prior_count));
                            entry.1 += 1; // total count
                            if base != pileup.ref_base {
                                entry.0 += 1; // error count
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Per-quality error rates from the accumulated counts of the top `num_clusters` clusters.
pub fn estimate_quality_error_rates(
    quality_stats: QualityStats,
    num_clusters: usize,
) -> HashMap<u8, f64> {
    log::info!(
        "Analyzing quality error rates from top {} clusters",
        num_clusters
    );

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
    log::debug!(
        "Quality Error Rate Histogram (from {} high-confidence positions)",
        total_bases_analyzed
    );
    log::debug!(
        "Overall error rate: {:.4}% ({}/{})",
        overall_error_rate * 100.0,
        total_errors,
        total_bases_analyzed
    );
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

    return quality_stats
        .iter()
        .map(|(&q, &(e, t))| {
            let rate = if t > 0 { e as f64 / t as f64 } else { 0.0 };
            (q, rate)
        })
        .collect();
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

    for (_, consensus) in consensuses.iter().enumerate() {
        let cluster = &consensus.cluster;
        if cluster.is_empty() {
            continue;
        }

        let representative = cluster[0];
        writeln!(
            writer,
            "{}_cluster_{}\tsize_{}\trepresentative_{}\tmembers\n{}",
            prefix,
            consensus.id,
            cluster.len(),
            representative,
            cluster
                .iter()
                .map(|x| format!(
                    "{} {}",
                    &twin_reads[*x].id,
                    &twin_reads[*x].est_id.unwrap_or(100.)
                ))
                .collect::<Vec<_>>()
                .join("\n")
        )?;
    }

    Ok(())
}

/// Write final_assignments.tsv: one row per read in twin_reads (assigned or unassigned).
/// Reads that have no valid EM alignment appear with asv="unassigned".
/// `debug_to_final` maps internal debug_id → final_asv_index (depth-sorted rank).
/// ASV labels match feature-table.tsv: `final_asv_N_depth_M` (or `..._depth_A-B-C` pooled).
pub fn write_final_assignments(
    twin_reads: &[crate::types::TwinRead],
    consensuses: &[ConsensusSequence],
    mappings_path: &std::path::Path,
    debug_to_final: &std::collections::HashMap<usize, usize>,
    output_path: &std::path::Path,
) -> std::io::Result<()> {
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader, BufWriter, Write};

    // Build final_idx → (asv_label, asv_len) from consensuses (already depth-sorted, ids reassigned)
    let idx_to_info: HashMap<usize, (String, usize)> = consensuses
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let depth_str = if c.per_sample_depths.is_empty() {
                (c.depth + c.appended_depth).to_string()
            } else {
                c.per_sample_depths
                    .iter()
                    .map(|d| d.to_string())
                    .collect::<Vec<_>>()
                    .join("-")
            };
            let label = format!("final_asv_{}_depth_{}", i, depth_str);
            let asv_len: usize = c.hp_lengths.iter().map(|&x| x as usize).sum();
            (i, (label, asv_len))
        })
        .collect();

    // Parse read_to_asv_mappings.tsv: keep only the first (best) row per read.
    // k-mer/EM path columns: read_id, debug_id:N, snpmer_mm, minimap2_nm
    // minimap2 path columns: read_id, debug_id:N, minimap2_nm
    //
    // The `debug_id:` field is located by prefix rather than by a fixed column index.
    // Read names can themselves contain tabs (dorado writes BAM tags like `qs:f:32.2`
    // tab-separated in the FASTQ header), which shifts every later column and would
    // otherwise make every read look unassigned.
    let mut read_to_asv: HashMap<String, (usize, i32)> = HashMap::new(); // read_id → (final_idx, nm)
    if mappings_path.exists() {
        let reader = BufReader::new(std::fs::File::open(mappings_path)?);
        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();
            let Some(tag_pos) = parts.iter().position(|p| p.starts_with("debug_id:")) else {
                continue;
            };
            if tag_pos == 0 || parts.len() <= tag_pos + 1 {
                continue;
            }
            // Everything before the tag is the read name; trailing fields are numeric.
            let read_id = parts[..tag_pos].join("\t");
            if read_to_asv.contains_key(&read_id) {
                continue;
            } // first row wins
            let Some(debug_id) = parts[tag_pos]
                .strip_prefix("debug_id:")
                .and_then(|s| s.parse::<usize>().ok())
            else {
                continue;
            };
            // k-mer/EM path appends (snpmer_mm, minimap2_nm); minimap2 path appends (nm).
            // Either way the minimap2 NM is the last field.
            let nm: i32 = parts[parts.len() - 1].parse().unwrap_or(0);
            if let Some(&final_idx) = debug_to_final.get(&debug_id) {
                read_to_asv.insert(read_id, (final_idx, nm));
            }
        }
    }

    let mut writer = BufWriter::new(std::fs::File::create(output_path)?);
    writeln!(
        writer,
        "read_id\tasv\talignment_identity\test_read_identity"
    )?;

    for read in twin_reads {
        let est = match read.est_id {
            Some(e) => format!("{:.4}", e),
            None => "NA".to_string(),
        };
        // `base_id` is the read name up to the first whitespace. The full `id` keeps any
        // trailing FASTQ header comment, which for dorado output contains tab-separated
        // BAM tags and would break this TSV's columns.
        match read_to_asv.get(&read.base_id) {
            Some(&(final_idx, nm)) => {
                let (label, asv_len) = idx_to_info
                    .get(&final_idx)
                    .map(|(l, len)| (l.as_str(), *len))
                    .unwrap_or(("unassigned", 1));
                let identity = if asv_len > 0 {
                    let nm_clamped = (nm as usize).min(asv_len);
                    (asv_len - nm_clamped) as f64 / asv_len as f64 * 100.0
                } else {
                    0.0
                };
                writeln!(
                    writer,
                    "{}\t{}\t{:.4}\t{}",
                    read.base_id, label, identity, est
                )?;
            }
            None => {
                writeln!(writer, "{}\tunassigned\tNA\t{}", read.base_id, est)?;
            }
        }
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
        let start_non = consensus_seq
            .iter()
            .enumerate()
            .find(|&(_i, &b)| b != b'N')
            .map(|(i, _)| i)
            .unwrap_or(0);
        let end_non = consensus_seq
            .iter()
            .enumerate()
            .rfind(|&(_i, &b)| b != b'N')
            .map(|(i, _)| i + 1)
            .unwrap_or(consensus_seq.len());
        let depth_field = if consensus.per_sample_depths.is_empty() {
            (consensus.depth + consensus.appended_depth).to_string()
        } else {
            consensus
                .per_sample_depths
                .iter()
                .map(|d| d.to_string())
                .collect::<Vec<_>>()
                .join("-")
        };
        let header = format!(">{}_asv_{}_depth_{} debug_id:{} chimera_score:{} unambiguous_read_assignments:{} ambig_read_assignments:{} num_align_leq_10_mismatches:{}",
            prefix, i, depth_field, consensus.id, consensus.chimera_score.unwrap_or(0), consensus.unambig_best_read_map_count.unwrap_or(0),
            consensus.ambig_read_map_count.unwrap_or(0), consensus.num_map_leq_10nm.unwrap_or(0));
        writeln!(writer, "{}", header)?;

        // Use decompressed sequence if available, otherwise use HPC sequence
        let sequence = String::from_utf8_lossy(&consensus_seq[start_non..end_non]);

        writeln!(writer, "{}", sequence)?;
    }

    Ok(())
}

/// Polish consensus sequences using Bayesian inference with quality-aware error rates
/// Trims low coverage ends and calculates posterior probabilities for each base
/// Analyzes one batch of pileups and applies the results to the corresponding
/// consensuses. `pileups[i]` describes the cluster `cluster_indices[i]`.
/// Call `collect_low_quality_consensuses` once after all batches are analyzed.
pub fn analyze_pileup_consensuses(
    mut pileups: Vec<ClusterPileup>,
    cluster_indices: &[usize],
    consensuses: &mut Vec<ConsensusSequence>,
    quality_error_map: &HashMap<u8, f64>,
    args: &Cli,
) {
    let bad_length_threshold = 100;
    let min_coverage_abs = (args.min_cluster_size * 3 / 4).max(2) as usize;
    let deletion_insertion_quality = 48u8; // Fixed quality for indels

    // Get error rate for deletions/insertions
    let indel_error_rate = quality_error_map
        .get(&deletion_insertion_quality)
        .copied()
        .unwrap_or(DEFAULT_ERR_RATE); // Default 2% if not in map

    // Select consensuses to debug: most abundant (highest depth), 10th percentile, 90th percentile
    let mut sorted_by_depth: Vec<(usize, usize)> = consensuses
        .iter()
        .enumerate()
        .map(|(idx, cons)| (idx, cons.depth))
        .collect();
    sorted_by_depth.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let debug_indices: Vec<usize> = vec![30];

    // Process each consensus
    for (local_idx, cluster) in pileups.iter_mut().enumerate() {
        let cluster_pileup = &mut cluster.positions;
        // Global consensus index, so logs and debug selection stay stable across batches
        let cluster_idx = cluster_indices[local_idx];
        let min_coverage =
            (cluster_pileup.iter().map(|p| p.depth()).max().unwrap_or(0) / 3).max(min_coverage_abs);
        if cluster_pileup.is_empty() {
            continue;
        }

        // 1. Trim low coverage ends
        let mut start_idx = 0;
        let mut end_idx = cluster_pileup.len();

        // Find first position with sufficient coverage
        let mut trimmed_depths = vec![];
        for (i, pileup) in cluster_pileup.iter().enumerate() {
            if pileup.depth() >= min_coverage {
                trimmed_depths.push(pileup.depth());
                start_idx = i;
                break;
            } else {
                trimmed_depths.push(pileup.depth());
            }
        }

        // Find last position with sufficient coverage
        for (i, pileup) in cluster_pileup.iter().enumerate().rev() {
            if pileup.depth() >= min_coverage {
                trimmed_depths.push(pileup.depth());
                end_idx = i + 1;
                break;
            } else {
                trimmed_depths.push(pileup.depth());
            }
        }

        if start_idx >= end_idx {
            log::warn!(
                "Consensus {} has no positions with sufficient coverage",
                cluster_idx
            );
            continue;
        }

        log::debug!("Consensus {}: Trimming from {}-{} to {}-{} with min depth {}. Trimmed depths at ends: {:?}",
            cluster_idx, 0, cluster_pileup.len(), start_idx, end_idx, min_coverage, trimmed_depths);

        // Move the retained positions in place and release the trimmed events.
        cluster_pileup.truncate(end_idx);
        drop(cluster_pileup.drain(..start_idx));

        // 2. Calculate posterior probabilities for each position
        let mut posterior_probs = Vec::new();

        for pileup in cluster_pileup.iter_mut() {
            let ref_base = pileup.ref_base;

            // Calculate log P(Z | ref = ref_base)
            let mut log_prob_ref = 0.0;

            // Calculate log P(Z | ref != ref_base)
            let mut log_prob_not_ref = 0.0;

            for base_entry in pileup.observations(&cluster.insertions) {
                match base_entry {
                    PileupBase::Base(obs_base, qual, _hp_len) => {
                        let error_rate = quality_error_map
                            .get(&qual)
                            .copied()
                            .unwrap_or(DEFAULT_ERR_RATE);
                        let accuracy = 1.0 - error_rate;

                        if obs_base == ref_base {
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
                        let first_qual = insertion_data
                            .first()
                            .map(|(_, q, _)| *q)
                            .unwrap_or(deletion_insertion_quality);
                        let error_rate = quality_error_map
                            .get(&first_qual)
                            .copied()
                            .unwrap_or(DEFAULT_ERR_RATE);
                        log_prob_not_ref += (1.0 - error_rate).ln();
                        log_prob_ref += error_rate.ln();

                        let error_rates: Vec<f64> = insertion_data
                            .iter()
                            .take(0)
                            .map(|(_, q, _)| {
                                quality_error_map
                                    .get(q)
                                    .copied()
                                    .unwrap_or(DEFAULT_ERR_RATE)
                            })
                            .collect();
                        log_prob_ref += error_rates.iter().map(|&er| er.ln()).sum::<f64>();
                        log_prob_not_ref +=
                            error_rates.iter().map(|&er| (1.0 - er).ln()).sum::<f64>();
                    }
                }
            }

            // Calculate posterior: P(ref | Z) = P(Z | ref) / (P(Z | ref) + P(Z | not ref))
            // In log space: log P(ref | Z) = log_prob_ref - log_sum_exp(log_prob_ref, log_prob_not_ref)
            let log_normalizer = log_sum_exp(log_prob_ref, log_prob_not_ref);
            let alt_posterior = log_prob_not_ref - log_normalizer;

            //if alt_posterior > -30.0 {
            let post_threshold = args
                .posterior_threshold_ln
                .min((args.min_cluster_size * 3) as f64);
            if alt_posterior > -post_threshold {
                log::debug!("Low posterior probability at consensus {}, covs: {:?},  position {}: alternate_posterior {:.6}, log_prob_ref {:.4}, log_prob_not_ref {:.4}, depth {}, range {}-{}",
                    cluster_idx, pileup.depth_nodeletion(),  pileup.ref_pos, alt_posterior, log_prob_ref, log_prob_not_ref, pileup.depth(), start_idx, end_idx);
                let ref_count = pileup
                    .observations(&cluster.insertions)
                    .filter(|b| matches!(b, PileupBase::Base(base, _, _) if *base == ref_base))
                    .count();
                let non_ref_base_count = pileup.depth() - ref_count;
                log::debug!(
                    "    Reference base: {}, Ref count: {}, Non-ref count: {}",
                    ref_base as char,
                    ref_count,
                    non_ref_base_count
                );
                //print top 20 bases
                let mut dbg_string = String::from("Bases: ");
                for base_entry in pileup.observations(&cluster.insertions).take(20) {
                    match base_entry {
                        PileupBase::Base(obs_base, qual, hp_len) => {
                            dbg_string.push_str(&format!(
                                " Base: {} (q={}, hp={}) ",
                                obs_base as char, qual, hp_len
                            ));
                        }
                        PileupBase::Deletion => {
                            dbg_string.push_str(" Deletion ");
                        }
                        PileupBase::Insertion(insertion_data) => {
                            let bases: Vec<u8> =
                                insertion_data.iter().map(|(b, _, _)| *b).collect();
                            dbg_string.push_str(&format!(
                                " Insertion: {} ",
                                String::from_utf8_lossy(&bases)
                            ));
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
        if debug_posterior {
            if debug_indices.contains(&cluster_idx) {
                log::debug!("=================================================================");
                log::debug!(
                    "Posterior probabilities for consensus {} (depth {})",
                    cluster_idx,
                    consensuses.get(cluster_idx).map(|c| c.depth).unwrap_or(0)
                );
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

    for (local_idx, cluster) in pileups.iter().enumerate() {
        // Global consensus index for this batch entry
        let i = cluster_indices[local_idx];
        let mut low_confidence_positions = vec![];
        for pileup in &cluster.positions {
            if let Some(_) = pileup.alt_posterior {
                low_confidence_positions.push(pileup.ref_pos);
            }
        }

        if cluster.positions.is_empty() {
            log::warn!("Consensus {} has empty pileup after polishing", i);
            continue;
        }
        let left_start = cluster.positions.first().map(|p| p.ref_pos).unwrap();
        let right_end = cluster.positions.last().map(|p| p.ref_pos).unwrap() + 1; // +1 because ref_pos is 0-based and end index is exclusive

        let start_polish = bad_length_threshold + left_start;
        let end_polish = right_end - bad_length_threshold;

        let low_conf_region_left = low_confidence_positions
            .iter()
            .filter(|&&pos| pos < start_polish)
            .map(|x| *x)
            .max()
            .unwrap_or(left_start);
        let low_conf_region_right = low_confidence_positions
            .iter()
            .filter(|&&pos| pos >= end_polish)
            .map(|x| *x)
            .min()
            .unwrap_or(right_end);

        let consensus = &mut consensuses[i];
        if low_conf_region_left > 0 {
            log::debug!(
                "Consensus {}: Masking low-confidence region at start up to position {}",
                i,
                low_conf_region_left
            );
            for pos in 0..low_conf_region_left {
                consensus.sequence[pos] = b'N';
            }
        }
        if low_conf_region_right < consensus.sequence.len() {
            log::debug!(
                "Consensus {}: Masking low-confidence region at end from position {} to {}",
                i,
                low_conf_region_right,
                consensus.sequence.len()
            );
            for pos in low_conf_region_right..consensus.sequence.len() {
                consensus.sequence[pos] = b'N';
            }
        }
        for pileup in &cluster.positions {
            if let Some(_) = pileup.alt_posterior {
                if args.mask_low_quality {
                    consensus.sequence[pileup.ref_pos] = b'N';
                }
                if pileup.ref_pos > low_conf_region_left && pileup.ref_pos < low_conf_region_right {
                    log::debug!(
                        "Consensus {}: Marking position {} as low quality, ends {}-{}",
                        i,
                        pileup.ref_pos,
                        low_conf_region_left,
                        low_conf_region_right
                    );
                    consensus.low_quality_positions.push(pileup.ref_pos);
                }
            }
        }
    }

}

/// Partitions `consensuses` on the low-quality criteria, writes the cluster TSVs, and
/// returns the low-quality consensuses (which are removed from `consensuses`).
/// Runs once after every pileup batch has been analyzed, since it needs the final
/// per-consensus quality annotations for all clusters.
pub fn collect_low_quality_consensuses(
    consensuses: &mut Vec<ConsensusSequence>,
    twin_reads: &[TwinRead],
    args: &Cli,
    temp_dir: &PathBuf,
) -> Vec<ConsensusSequence> {
    let low_quality_consensuses = consensuses
        .iter()
        .filter(|c| lq_criteria(c, args))
        .map(|c| c.clone())
        .collect::<Vec<_>>();
    let mut lq_ids: Vec<usize> = low_quality_consensuses.iter().map(|c| c.id).collect();
    lq_ids.sort_unstable();
    log::info!("Low quality consensus sequences: {:?}", lq_ids);
    let low_quality_cluster_file = temp_dir.join("low_quality_clusters.tsv");
    write_clusters_tsv(
        &low_quality_consensuses,
        twin_reads,
        &low_quality_cluster_file,
        "low_quality",
    )
    .expect("Failed to write low_quality_clusters.tsv");
    log::info!("Wrote low quality cluster information to low_quality_clusters.tsv");
    consensuses.retain(|c| !lq_criteria(c, args));

    // Write clusters before filtering out low quality consensuses
    let postfilter_file = temp_dir.join("clusters_after_quality_filter_stage4.tsv");
    write_clusters_tsv(consensuses, twin_reads, &postfilter_file, "prefilter")
        .expect("Failed to write clusters_after_quality_filter_stage4.tsv");
    log::info!(
        "Wrote cluster information before filtering to clusters_after_quality_filter_stage4.tsv"
    );

    log::info!("Polishing complete");

    return low_quality_consensuses;
}

fn lq_criteria(consensus: &ConsensusSequence, args: &Cli) -> bool {
    (consensus.low_quality_positions.len() > 0)
        && (consensus.depth
            / (consensus.low_quality_positions.len() * consensus.low_quality_positions.len())
            < args.n_depth_cutoff)
}

fn remove_similar_seqs_kmers(mut consensuses: Vec<ConsensusSequence>) -> Vec<ConsensusSequence> {
    let mut kmer_index = std::collections::HashMap::new();
    let mut filtered_consensuses: Vec<ConsensusSequence> = Vec::new();
    let mut consensus_id_to_minis = std::collections::HashMap::new();
    let adapter_buffer = 25;

    for (i, consensus) in consensuses.iter().enumerate() {
        if consensus.sequence.len() < 100 {
            continue;
        }
        let mut minimizers = vec![];
        let mut positions = vec![];
        seeding::minimizer_seeds_positions(
            &consensus.sequence[adapter_buffer..consensus.sequence.len() - adapter_buffer],
            &mut minimizers,
            &mut positions,
            10,
            21,
        );
        for &mini in minimizers.iter() {
            kmer_index.entry(mini).or_insert_with(Vec::new).push(i);
        }
        consensus_id_to_minis.insert(i, minimizers);
    }

    for enum_id in 0..consensuses.len() {
        let Some(minimizers) = consensus_id_to_minis.get(&enum_id) else {
            continue;
        };
        let mut possible_greater_ids = std::collections::HashSet::new();
        let mut first = true;
        for mini in minimizers {
            if first {
                if let Some(ids) = kmer_index.get(mini) {
                    for id in ids {
                        if consensuses[*id].depth / 2 > consensuses[enum_id].depth {
                            possible_greater_ids.insert(*id);
                        }
                    }
                }
            } else {
                if let Some(ids) = kmer_index.get(mini) {
                    let id_set = ids
                        .iter()
                        .cloned()
                        .collect::<std::collections::HashSet<usize>>();
                    possible_greater_ids = possible_greater_ids
                        .intersection(&id_set)
                        .cloned()
                        .collect();
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
    log::info!(
        "Reduced consensus sequences from {} to {} after k-mer based deduplication",
        prev_size,
        consensuses.len()
    );

    // Write polished consensus sequences to FASTA for indexing
    let output_fasta_path = temp_dir.join("polished_consensuses.fasta");
    write_consensus_fasta(&consensuses, &output_fasta_path, "polished")
        .expect("Failed to write polished_consensuses.fasta");
    log::info!(
        "Wrote {} polished consensus sequences to polished_consensuses.fasta",
        consensuses.len()
    );

    // Build aligner using the first consensus as reference
    let mut aligner = Aligner::builder()
        .lrhq()
        .with_index_threads(args.threads)
        .with_cigar()
        .with_index(output_fasta_path.to_str().unwrap(), None)
        .expect("Failed to create aligner");

    aligner.mapopt.set_no_diag();
    aligner.mapopt.best_n = 75;

    // Store mappings: (query_idx, target_idx, nm, target_depth)
    let mappings = Mutex::new(Vec::new());

    // First, merge low quality consensuses into high quality consensuses
    log::info!(
        "Merging {} low quality consensuses into high quality consensuses",
        low_qual_consensuses.len()
    );
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
    let mut low_qual_mappings = low_qual_mappings.into_inner().unwrap();
    low_qual_mappings.sort_unstable();
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
        let valid_targets: Vec<(usize, usize, usize)> = mappings
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

                (relative_depth < threshold) || (1.0 / relative_depth < threshold)
            })
            .map(|(_, t_idx, nm, t_depth)| (*t_idx, *nm, *t_depth))
            .collect();

        // If there are valid targets, choose the one with highest depth
        if !valid_targets.is_empty() {
            let mut query_to_ref_mappings = vec![];
            let mut ref_to_query_mappings = vec![];
            for (t_idx, nm, t_depth) in &valid_targets {
                if consensuses[*t_idx].depth == consensuses[query_idx].depth {
                    if *nm == 0 {
                        // Perfect match with identical depth, merge based on index to avoid circular merges
                        if query_idx > *t_idx {
                            merge_map.insert(query_idx, *t_idx);
                        }
                    }
                    continue; // Skip identical depth consensuses here to avoid circular merges
                } else if consensuses[*t_idx].depth > consensuses[query_idx].depth {
                    query_to_ref_mappings.push((*t_idx, *nm, *t_depth, query_idx));
                } else {
                    ref_to_query_mappings.push((
                        query_idx,
                        *nm,
                        consensuses[query_idx].depth,
                        *t_idx,
                    ));
                }
            }
            if query_to_ref_mappings.len() > 0 {
                query_to_ref_mappings.sort_by(|a, b| {
                    b.2.cmp(&a.2)
                        .then_with(|| a.1.cmp(&b.1))
                        .then_with(|| a.0.cmp(&b.0))
                });
                let best_target = query_to_ref_mappings[0].0;
                merge_map.insert(query_idx, best_target);
            }

            for (_, _, _, t_idx) in ref_to_query_mappings {
                if !merge_map.contains_key(&t_idx) {
                    merge_map.insert(t_idx, query_idx);
                }
            }
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
    for query_idx in 0..consensuses.len() {
        let Some(&target_idx) = merged_into.get(&query_idx) else {
            continue;
        };
        log::debug!(
            "Merging consensus {} (depth {}) into consensus {} (depth {})",
            consensuses[query_idx].id,
            consensuses[query_idx].depth,
            consensuses[target_idx].id,
            consensuses[target_idx].depth,
        );

        // Move all reads from query cluster to target cluster
        let reads_to_move = new_clusters[query_idx].clone();
        new_clusters[target_idx].extend(reads_to_move);
        new_clusters[query_idx].clear();
    }

    for cluster in &mut new_clusters {
        cluster.sort_unstable();
    }

    // Build new consensus sequences with updated depths
    let mut new_consensuses = Vec::new();

    for (idx, consensus) in consensuses.into_iter().enumerate() {
        if !new_clusters[idx].is_empty() {
            let new_depth = new_clusters[idx].len();
            let new_cluster = new_clusters[idx].clone();
            let mut new_cons = ConsensusSequence::new(
                consensus.sequence,
                consensus.hp_lengths,
                new_depth,
                consensus.id,
                new_cluster,
            );
            new_cons.decompress();
            new_consensuses.push(new_cons);
        }
    }

    log::info!(
        "Consensus merging complete: {} -> {} consensuses",
        new_clusters.len(),
        new_consensuses.len()
    );

    new_consensuses.sort_by(|a, b| {
        b.depth
            .cmp(&a.depth)
            .then_with(|| a.sequence.cmp(&b.sequence))
            .then_with(|| a.id.cmp(&b.id))
    });

    let final_file = temp_dir.join("final_clusters_merged_stage5.tsv");
    write_clusters_tsv(&new_consensuses, twin_reads, &final_file, "final")
        .expect("Failed to write final_clusters_merged_stage5.tsv");

    // Write merged consensus sequences (before chimera filtering)
    let merged_fasta = temp_dir.join("merged_consensus_sequences.fasta");
    write_consensus_fasta(&new_consensuses, &merged_fasta, "merged")
        .expect("Failed to write merged_consensus_sequences.fasta");
    log::info!(
        "Wrote {} merged consensus sequences to merged_consensus_sequences.fasta",
        new_consensuses.len()
    );

    new_consensuses
}

/// Equivalence class: a set of ASVs that a group of reads maps to with equal quality
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct EquivalenceClass {
    asv_indices: Vec<usize>,
}

fn sorted_equivalence_classes(
    classes: HashMap<EquivalenceClass, usize>,
) -> Vec<(EquivalenceClass, usize)> {
    let mut classes: Vec<_> = classes.into_iter().collect();
    classes.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    classes
}

#[cfg(test)]
mod determinism_tests {
    use super::*;

    #[test]
    fn equivalence_classes_have_canonical_iteration_order() {
        let mut classes = HashMap::new();
        classes.insert(
            EquivalenceClass {
                asv_indices: vec![3],
            },
            4,
        );
        classes.insert(
            EquivalenceClass {
                asv_indices: vec![0, 2],
            },
            7,
        );
        classes.insert(
            EquivalenceClass {
                asv_indices: vec![0, 1],
            },
            5,
        );

        let sorted = sorted_equivalence_classes(classes);
        let keys: Vec<Vec<usize>> = sorted
            .into_iter()
            .map(|(class, _)| class.asv_indices)
            .collect();

        assert_eq!(keys, vec![vec![0, 1], vec![0, 2], vec![3]]);
    }
}

/// Refine ASV depths using minimap2 all-vs-all mapping (used in low-polymorphism mode).
/// Skips the SNPmer index entirely; maps every read to the ASV FASTA with minimap2.
fn refine_asv_depths_with_minimap2(
    twin_reads: &[TwinRead],
    consensuses: &mut Vec<ConsensusSequence>,
    args: &Cli,
    temp_dir: &PathBuf,
) {
    if consensuses.is_empty() {
        log::warn!("No consensuses to refine (minimap2 path)");
        return;
    }

    let mapping_file_writer = Mutex::new(std::io::BufWriter::new(
        std::fs::File::create(temp_dir.join("read_to_asv_mappings.tsv"))
            .expect("Failed to create read_to_asv_mappings.tsv"),
    ));

    log::info!("Low-polymorphism mode: refining ASV depths via minimap2 all-to-all mapping");

    // Write ASV sequences to FASTA (reuse the same path as the normal EM path)
    let asv_fasta_path = temp_dir.join("final_asvs_for_em.fasta");
    write_consensus_fasta(consensuses, &asv_fasta_path, "em_refinement")
        .expect("Failed to write ASVs for EM refinement (minimap2 path)");

    // Build a multi-sequence minimap2 index over all ASV sequences.
    // target_id in mapping results directly corresponds to the ASV index (0-based order in FASTA).
    let aligner = Aligner::builder()
        .lrhq()
        .with_index_threads(args.threads)
        .with_cigar()
        .with_index(asv_fasta_path.to_str().unwrap(), None)
        .expect("Failed to build minimap2 index for ASVs");

    log::info!(
        "Mapping {} reads to {} ASVs via minimap2",
        twin_reads.len(),
        consensuses.len()
    );

    let eq_classes: Mutex<HashMap<EquivalenceClass, usize>> = Mutex::new(HashMap::new());
    let filtered_reads_count = Mutex::new(0usize);
    let total_assigned_reads = Mutex::new(0usize);

    for cons in consensuses.iter_mut() {
        cons.unambig_best_read_map_count = Some(0);
        cons.ambig_read_map_count = Some(0);
        cons.num_map_leq_10nm = Some(0);
    }

    let unambig_read_map_count = Mutex::new(vec![0usize; consensuses.len()]);
    let ambig_read_map_count = Mutex::new(vec![0usize; consensuses.len()]);
    let num_map_leq_10nm = Mutex::new(vec![0usize; consensuses.len()]);

    twin_reads.par_iter().for_each(|twin_read| {
        let seq: Vec<u8> = twin_read
            .dna_seq
            .iter()
            .map(|x| x.to_char().to_ascii_uppercase() as u8)
            .collect();
        let mappings = aligner
            .map(&seq, true, false, None, None, None)
            .unwrap_or_default();

        // Keep only primary/supplementary hits with mapq > 0
        let valid: Vec<_> = mappings
            .iter()
            .filter(|m| m.mapq > 0 && m.alignment.is_some())
            .collect();

        if valid.is_empty() {
            *filtered_reads_count.lock().unwrap() += 1;
            return;
        }

        // Find the best NM among valid hits
        let best_nm = valid
            .iter()
            .map(|m| m.alignment.as_ref().unwrap().nm)
            .min()
            .unwrap();

        // Collect all hits tied at the best NM
        let mut best_asv_indices: Vec<usize> = valid
            .iter()
            .filter(|m| m.alignment.as_ref().unwrap().nm == best_nm)
            .map(|m| m.target_id as usize)
            .collect();
        best_asv_indices.sort();
        best_asv_indices.dedup();

        {
            // Format outside the lock; this read's rows still go out contiguously.
            let mut row_buf = String::new();
            for &asv_idx in &best_asv_indices {
                use std::fmt::Write as _;
                let _ = writeln!(
                    row_buf,
                    "{}\tdebug_id:{}\t{}",
                    twin_read.base_id, consensuses[asv_idx].id, best_nm
                );
            }
            let mut writer = mapping_file_writer.lock().unwrap();
            writer
                .write_all(row_buf.as_bytes())
                .expect("Failed to write read_to_asv_mappings.tsv");
        }

        if best_asv_indices.len() == 1 {
            let asv_idx = best_asv_indices[0];
            unambig_read_map_count.lock().unwrap()[asv_idx] += 1;
        } else {
            for &asv_idx in &best_asv_indices {
                ambig_read_map_count.lock().unwrap()[asv_idx] += 1;
            }
        }

        if best_nm <= 10 {
            for &asv_idx in &best_asv_indices {
                num_map_leq_10nm.lock().unwrap()[asv_idx] += 1;
            }
        }

        let eq_class = EquivalenceClass {
            asv_indices: best_asv_indices,
        };
        *eq_classes.lock().unwrap().entry(eq_class).or_insert(0) += 1;
        *total_assigned_reads.lock().unwrap() += 1;
    });

    let eq_classes = sorted_equivalence_classes(eq_classes.into_inner().unwrap());
    let filtered_reads = filtered_reads_count.into_inner().unwrap();
    let total_assigned = total_assigned_reads.into_inner().unwrap();

    log::info!(
        "Filtered {} reads with no valid minimap2 mapping",
        filtered_reads
    );
    log::info!("Total assigned reads: {}", total_assigned);
    log::info!(
        "Total percentage assigned: {:.2}%",
        (total_assigned as f64 / (total_assigned + filtered_reads) as f64) * 100.0
    );

    let unambig_read_map_count = unambig_read_map_count.into_inner().unwrap();
    let ambig_read_map_count = ambig_read_map_count.into_inner().unwrap();
    let num_map_leq_10nm = num_map_leq_10nm.into_inner().unwrap();
    for i in 0..consensuses.len() {
        consensuses[i].unambig_best_read_map_count = Some(unambig_read_map_count[i]);
        consensuses[i].ambig_read_map_count = Some(ambig_read_map_count[i]);
        consensuses[i].num_map_leq_10nm = Some(num_map_leq_10nm[i]);
    }

    if eq_classes.is_empty() {
        log::warn!("No reads mapped to any ASV (minimap2 path). Keeping original depths.");
        return;
    }

    // EM algorithm — identical to the SNPmer path
    let num_asvs = consensuses.len();
    let mut asv_abundances = vec![1.0 / num_asvs as f64; num_asvs];
    let convergence_threshold = 0.01 / total_assigned as f64;

    log::info!(
        "Running EM algorithm with convergence threshold: {:.6e}",
        convergence_threshold
    );

    let mut iteration = 0;
    const MAX_ITERATIONS: usize = 10000;

    loop {
        iteration += 1;
        let mut new_asv_abundances = vec![0.0; num_asvs];

        for (eq_class, count) in &eq_classes {
            let denominator: f64 = eq_class
                .asv_indices
                .iter()
                .map(|&asv_idx| asv_abundances[asv_idx])
                .sum();

            if denominator > 0.0 {
                for &asv_idx in &eq_class.asv_indices {
                    let contribution = (*count as f64) * asv_abundances[asv_idx] / denominator;
                    new_asv_abundances[asv_idx] += contribution;
                }
            }
        }

        let total_counts: f64 = new_asv_abundances.iter().sum();
        if total_counts > 0.0 {
            for abundance in new_asv_abundances.iter_mut() {
                *abundance /= total_assigned as f64;
            }
        }

        let max_change = asv_abundances
            .iter()
            .zip(new_asv_abundances.iter())
            .map(|(old, new)| (old - new).abs())
            .fold(0.0, f64::max);

        asv_abundances = new_asv_abundances;

        if max_change < convergence_threshold || iteration >= MAX_ITERATIONS {
            log::info!(
                "EM converged after {} iterations (max change: {:.6e})",
                iteration,
                max_change
            );
            break;
        }
    }

    log::info!("Updating ASV depths based on EM abundances");
    for (asv_idx, consensus) in consensuses.iter_mut().enumerate() {
        let em_abundance = asv_abundances[asv_idx];
        let new_depth = (em_abundance * total_assigned as f64).round() as usize;
        if new_depth != 0 && consensus.depth / new_depth > 10 {
            log::debug!("ASV {} possible removal due to coverage drop of 10x: original depth {}, new depth {}",
                asv_idx, consensus.depth, new_depth);
        }
        consensus.depth = new_depth;
    }

    let original_count = consensuses.len();
    consensuses.retain(|c| c.depth > 0);
    let filtered_count = original_count - consensuses.len();
    if filtered_count > 0 {
        log::info!(
            "Filtered {} ASVs with zero depth after EM refinement (minimap2 path)",
            filtered_count
        );
    }
    log::info!(
        "EM refinement complete (minimap2 path): {} ASVs remaining",
        consensuses.len()
    );
}

/// The read's lowest-mismatch candidate ASVs among those passing the minimizer-overlap and
/// mismatch-ratio filters, in `candidate_stats` order, plus the lowest ratio seen (for
/// logging). Candidates are visited by increasing SNPmer mismatches: once a mismatch level
/// has a passing ASV, higher levels could never be kept (only the fewest-mismatch passers
/// are), so their minimizer-set intersections, the expensive part, are skipped.
fn lowest_mismatch_candidates(
    candidate_stats: FxHashMap<usize, (usize, usize)>,
    read_minimizers: &FxHashSet<Kmer48>,
    asv_minimizer_sets: &[FxHashSet<Kmer48>],
    k: usize,
    c: usize,
) -> (Vec<(usize, usize)>, Option<(f64, usize, usize)>) {
    const RATIO_THRESHOLD: f64 = 0.0050;
    let mut ordered: Vec<(usize, usize, usize)> = candidate_stats
        .into_iter()
        .enumerate()
        .map(|(pos, (asv_idx, (_matches, mismatches)))| (mismatches, pos, asv_idx))
        .collect();
    ordered.sort_unstable();
    let mut best = Vec::new();
    let mut min_ratio: Option<(f64, usize, usize)> = None;
    let mut i = 0;
    while i < ordered.len() && best.is_empty() {
        let level = ordered[i].0;
        while i < ordered.len() && ordered[i].0 == level {
            let (mismatches, _, asv_idx) = ordered[i];
            i += 1;
            let asv_minimizers = &asv_minimizer_sets[asv_idx];
            let minimizer_matches = read_minimizers.intersection(asv_minimizers).count();
            if minimizer_matches == 0 {
                continue;
            }
            if (minimizer_matches as f64
                / read_minimizers.len().min(asv_minimizers.len()) as f64)
                < (0.950f64).powi(k as i32)
            {
                continue;
            }
            let ratio = mismatches as f64 / minimizer_matches as f64 / c as f64;
            if min_ratio.map_or(true, |m| ratio < m.0) {
                min_ratio = Some((ratio, mismatches, minimizer_matches));
            }
            if ratio <= RATIO_THRESHOLD {
                best.push((asv_idx, mismatches));
            }
        }
    }
    (best, min_ratio)
}

/// Refine ASV depths using EM algorithm on read-level mappings
/// Maps all reads back to final ASVs using k-mer/SNPmer-based approach
/// With `record_hits`, also returns each read's result for Stage 7b (per-sample counts).
pub fn refine_asv_depths_with_em(
    twin_reads: &[TwinRead],
    consensuses: &mut Vec<ConsensusSequence>,
    kmer_info: &KmerGlobalInfo,
    args: &Cli,
    temp_dir: &PathBuf,
    record_hits: bool,
) -> Option<ReadAsvHits> {
    if args.low_polymorphism {
        refine_asv_depths_with_minimap2(twin_reads, consensuses, args, temp_dir);
        return None;
    }

    if consensuses.is_empty() {
        log::warn!("No consensuses to refine");
        return None;
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

    let asv_twin_reads = kmer_comp::twin_reads_from_fasta(
        &asv_fasta_path,
        kmer_info,
        args.kmer_size,
        args.c,
        args.blockmer_length,
        args.minimum_base_quality,
    );
    log::info!(
        "Loaded {} ASVs as TwinReads for k-mer comparison",
        asv_twin_reads.len()
    );

    // Step 2: Build SNPmer index for all ASVs
    let k = args.kmer_size;
    let mask = !(3 << (k - 1)); // Mask for creating splitmers
    let mut asv_snpmer_index: FxHashMap<u64, Vec<(usize, Kmer48)>> = FxHashMap::default();

    for (asv_idx, asv_read) in asv_twin_reads.iter().enumerate() {
        let asv_snpmers = asv_read.snpmers_vec();
        for (_pos, kmer) in asv_snpmers {
            let splitmer = kmer.to_u64() & mask;
            asv_snpmer_index
                .entry(splitmer)
                .or_insert_with(Vec::new)
                .push((asv_idx, kmer));
        }
    }

    log::info!(
        "Built SNPmer index with {} unique splitmers",
        asv_snpmer_index.len()
    );

    // Precompute the per-ASV minimizer sets and uppercased sequences once. Both were
    // previously rebuilt inside the per-read loop for every candidate ASV, i.e. millions
    // of times over a run, even though the ASV set is fixed here.
    let asv_minimizer_sets: Vec<FxHashSet<Kmer48>> = asv_twin_reads
        .iter()
        .map(|t| t.minimizer_kmers().into_iter().collect())
        .collect();
    let asv_seqs_u8: Vec<Vec<u8>> = asv_twin_reads
        .iter()
        .map(|t| {
            t.dna_seq
                .iter()
                .map(|x| x.to_char().to_ascii_uppercase() as u8)
                .collect()
        })
        .collect();

    log::info!(
        "Mapping {} reads to {} ASVs using k-mer comparison",
        twin_reads.len(),
        consensuses.len()
    );

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

    // Each read's lowest-mismatch ASVs (with NM ranks), interned, for Stage 7b.
    let read_set: Option<Vec<AtomicU32>> = record_hits
        .then(|| (0..twin_reads.len()).map(|_| AtomicU32::new(NO_HITS)).collect());
    let hit_sets: Mutex<(FxHashMap<Vec<(u32, u32)>, u32>, Vec<Vec<(u32, u32)>>)> =
        Mutex::new(Default::default());

    twin_reads
        .par_iter()
        .enumerate()
        .for_each(|(read_idx, twin_read)| {
            let read_snpmers = twin_read.snpmer_kmers();
            let read_minimizers: FxHashSet<Kmer48> =
                twin_read.minimizer_kmers().into_iter().collect();

            // Compare read against each ASV using SNPmers
            let candidate_stats = find_compatible_candidates(&asv_snpmer_index, &read_snpmers, k);
            let (best_asv_indices, min_ratio) = lowest_mismatch_candidates(
                candidate_stats,
                &read_minimizers,
                &asv_minimizer_sets,
                k,
                args.c,
            );
            if best_asv_indices.is_empty() {
                *filtered_reads_count.lock().unwrap() += 1;
                return;
            }
            let (min_ratio, min_mismatches, max_mini) = min_ratio.unwrap();

            let lowest: Vec<usize> = best_asv_indices.iter().map(|&(asv_idx, _)| asv_idx).collect();

            // For each best ASV, count number of kmer matches for tie-breaking
            let mut best_alns: Vec<(usize, i32, usize)> = Vec::new(); // (asv_idx, nm, mismatches)
            let seq_u8: Vec<u8> = twin_read
                .dna_seq
                .iter()
                .map(|x| x.to_char().to_ascii_uppercase() as u8)
                .collect();
            // This index covers one read and serves only the few candidate ASVs below,
            // so index construction is a large share of the cost here. See
            // SHORT_REF_BUCKET_BITS.
            let mut builder = Aligner::builder().lrhq();
            builder.idxopt.bucket_bits = SHORT_REF_BUCKET_BITS;
            let aligner = builder
                .with_cigar()
                .with_seq(&seq_u8)
                .expect("Failed to create aligner");

            for (asv_idx, mismatches) in best_asv_indices.into_iter() {
                let seq_u8_asv = &asv_seqs_u8[asv_idx];
                let alignment_result = aligner
                    .map(seq_u8_asv, true, false, None, None, None)
                    .unwrap();
                if alignment_result.is_empty() {
                    continue;
                }
                best_alns.push((
                    asv_idx,
                    alignment_result[0].alignment.as_ref().unwrap().nm,
                    mismatches,
                ));
            }

            // Sort best alignments by number of mismatches (ascending)
            best_alns.sort_by(|a, b| a.1.cmp(&b.1));

            if let Some(read_set) = &read_set {
                let mut nms: Vec<i32> = best_alns.iter().map(|x| x.1).collect();
                nms.sort_unstable();
                nms.dedup();
                let mut set: Vec<(u32, u32)> = lowest
                    .iter()
                    .map(|&asv_idx| {
                        let rank = best_alns.iter().find(|x| x.0 == asv_idx).map_or(UNALIGNED, |x| {
                            nms.binary_search(&x.1).unwrap() as u32
                        });
                        (asv_idx as u32, rank)
                    })
                    .collect();
                set.sort_unstable();
                let mut sets = hit_sets.lock().unwrap();
                let next = sets.1.len() as u32;
                let id = *sets.0.entry(set.clone()).or_insert(next);
                if id == next {
                    sets.1.push(set);
                }
                read_set[read_idx].store(id, Ordering::Relaxed);
            }

            let best_nm = best_alns.first().map(|x| x.1).unwrap_or(i32::MAX);
            let mut best_asv_aln_indices: Vec<usize> = best_alns
                .iter()
                .filter(|(_, nm, _)| *nm == best_nm)
                .map(|(asv_idx, _, _)| *asv_idx)
                .collect::<Vec<usize>>();

            // Map to file. Format outside the lock, then emit this read's rows in one
            // write: the rows for a read stay contiguous and in the same order, but the
            // lock is no longer held across five separate formatting calls.
            {
                let mut row_buf = String::new();
                for (asv_idx, mini_matches, mismatches) in best_alns.iter().take(5) {
                    use std::fmt::Write as _;
                    let _ = writeln!(
                        row_buf,
                        "{}\tdebug_id:{}\t{}\t{}",
                        twin_read.base_id, consensuses[*asv_idx].id, mismatches, mini_matches
                    );
                }
                let mut mapping_file_writer = mapping_file_writer.lock().unwrap();
                mapping_file_writer
                    .write_all(row_buf.as_bytes())
                    .expect("Failed to write to read_to_asv_mappings.tsv");
            }

            // Only keep reads that have at least one good mapping
            if !best_asv_aln_indices.is_empty() {
                // Sort to ensure consistent equivalence class keys
                best_asv_aln_indices.sort();

                let eq_class = EquivalenceClass {
                    asv_indices: best_asv_aln_indices,
                };

                if eq_class.asv_indices.len() == 1 {
                    let asv_idx = eq_class.asv_indices[0];
                    let mut unambig_counts = unambig_read_map_count.lock().unwrap();
                    unambig_counts[asv_idx] += 1;
                } else {
                    for &asv_idx in eq_class.asv_indices.iter() {
                        let mut ambig_counts = ambig_read_map_count.lock().unwrap();
                        ambig_counts[asv_idx] += 1;
                    }
                }

                if best_nm <= 10 {
                    for &asv_idx in eq_class.asv_indices.iter() {
                        let mut leq10nm_counts = num_map_leq_10nm.lock().unwrap();
                        leq10nm_counts[asv_idx] += 1;
                    }
                }

                // Add to equivalence class with count
                let mut eq_map = eq_classes.lock().unwrap();
                *eq_map.entry(eq_class).or_insert(0) += 1;
                *total_assigned_reads.lock().unwrap() += 1;
            } else {
                log::trace!(
                    "Read filtered out due to high ratio mappings: ratio {:.4} mismatch {} mini {}",
                    min_ratio,
                    min_mismatches,
                    max_mini
                );
                *filtered_reads_count.lock().unwrap() += 1;
            }
        });

    // Captured before EM drops zero-depth ASVs: the hits refer to this order.
    let hits = read_set.map(|read_set| ReadAsvHits {
        asv_ids: consensuses.iter().map(|c| c.id).collect(),
        read_set: read_set.into_iter().map(AtomicU32::into_inner).collect(),
        sets: hit_sets.into_inner().unwrap().1,
    });

    let eq_classes = sorted_equivalence_classes(eq_classes.into_inner().unwrap());
    let filtered_reads = filtered_reads_count.into_inner().unwrap();
    let total_assigned = total_assigned_reads.into_inner().unwrap();

    log::info!("Filtered {} reads with ratio > 0.005", filtered_reads);
    log::info!("Total assigned reads: {}", total_assigned);
    log::info!(
        "Total percentage assigned: {:.2}%",
        (total_assigned as f64 / (total_assigned + filtered_reads) as f64) * 100.0
    );
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
        log::trace!(
            "Equivalence class: ASVs {:?}, Count {}",
            eq_class.asv_indices,
            count
        );
    }

    if eq_classes.is_empty() {
        log::warn!("No reads mapped well to ASVs. Keeping original depths.");
        return hits;
    }

    // Step 4: Run EM algorithm on equivalence classes
    let num_asvs = consensuses.len();
    let mut asv_abundances = vec![1.0 / num_asvs as f64; num_asvs];
    let convergence_threshold = 0.01 / total_assigned as f64;

    log::info!(
        "Running EM algorithm with convergence threshold: {:.6e}",
        convergence_threshold
    );

    let mut iteration = 0;
    const MAX_ITERATIONS: usize = 10000;

    loop {
        iteration += 1;
        let mut new_asv_abundances = vec![0.0; num_asvs];

        // E-step + M-step: distribute reads proportionally based on current abundances
        for (eq_class, count) in &eq_classes {
            let denominator: f64 = eq_class
                .asv_indices
                .iter()
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
        let max_change = asv_abundances
            .iter()
            .zip(new_asv_abundances.iter())
            .map(|(old, new)| (old - new).abs())
            .fold(0.0, f64::max);

        asv_abundances = new_asv_abundances;

        if max_change < convergence_threshold || iteration >= MAX_ITERATIONS {
            log::info!(
                "EM converged after {} iterations (max change: {:.6e})",
                iteration,
                max_change
            );
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
        } else {
            log::debug!(
                "ASV {}: Original depth = {}, EM abundance = {:.6}, New depth = {}",
                asv_idx,
                consensus.depth,
                em_abundance,
                new_depth
            );
            consensus.depth = new_depth;
        }
    }

    // Filter out ASVs with zero depth after EM
    let original_count = consensuses.len();
    consensuses.retain(|c| c.depth > 0);
    let filtered_count = original_count - consensuses.len();

    if filtered_count > 0 {
        log::info!(
            "Filtered {} ASVs with zero depth after EM refinement",
            filtered_count
        );
    }

    log::info!(
        "EM refinement complete: {} ASVs remaining",
        consensuses.len()
    );
    hits
}

/// What Stage 7 found for each read, kept so Stage 7b can split counts by sample without
/// mapping every read again.
pub struct ReadAsvHits {
    /// Debug ids of the ASVs that the indices in `sets` refer to (Stage 7's order).
    asv_ids: Vec<usize>,
    /// Per read: index into `sets`, or `NO_HITS` if no ASV passed the mapping filters.
    read_set: Vec<u32>,
    /// Distinct per-read results: the lowest-mismatch ASVs, each with the rank of its
    /// alignment NM within the read (0 = best), or `UNALIGNED`; sorted by ASV.
    sets: Vec<Vec<(u32, u32)>>,
}

const NO_HITS: u32 = u32::MAX;
const UNALIGNED: u32 = u32::MAX;

/// Stage 7b's per-read mapping against the final ASVs (identical to Stage 7's per-read logic).
struct SampleMapper {
    k: usize,
    snpmer_index: FxHashMap<u64, Vec<(usize, Kmer48)>>,
    minimizer_sets: Vec<FxHashSet<Kmer48>>,
    seqs_u8: Vec<Vec<u8>>,
}

impl SampleMapper {
    fn new(kmer_info: &KmerGlobalInfo, args: &Cli, asv_fasta_path: &Path) -> Self {
        let asv_twin_reads = kmer_comp::twin_reads_from_fasta(
            asv_fasta_path,
            kmer_info,
            args.kmer_size,
            args.c,
            args.blockmer_length,
            args.minimum_base_quality,
        );
        let k = args.kmer_size;
        let mask = !(3u64 << (k - 1));
        let mut snpmer_index: FxHashMap<u64, Vec<(usize, Kmer48)>> = FxHashMap::default();
        for (asv_idx, asv_read) in asv_twin_reads.iter().enumerate() {
            for (_pos, kmer) in asv_read.snpmers_vec() {
                snpmer_index
                    .entry(kmer.to_u64() & mask)
                    .or_insert_with(Vec::new)
                    .push((asv_idx, kmer));
            }
        }
        let minimizer_sets = asv_twin_reads
            .iter()
            .map(|t| t.minimizer_kmers().into_iter().collect())
            .collect();
        let seqs_u8 = asv_twin_reads
            .iter()
            .map(|t| {
                t.dna_seq
                    .iter()
                    .map(|x| x.to_char().to_ascii_uppercase() as u8)
                    .collect()
            })
            .collect();
        SampleMapper { k, snpmer_index, minimizer_sets, seqs_u8 }
    }

    /// The read's equivalence class, or None if it is filtered.
    fn class_of(&self, twin_read: &TwinRead, args: &Cli) -> Option<EquivalenceClass> {
        let k = self.k;
        let read_snpmers = twin_read.snpmer_kmers();
        let read_minimizers: FxHashSet<Kmer48> = twin_read.minimizer_kmers().into_iter().collect();

        let candidate_stats = find_compatible_candidates(&self.snpmer_index, &read_snpmers, k);
        let (best_asv_indices, _) =
            lowest_mismatch_candidates(candidate_stats, &read_minimizers, &self.minimizer_sets, k, args.c);
        if best_asv_indices.is_empty() {
            return None;
        }

        let seq_u8: Vec<u8> = twin_read
            .dna_seq
            .iter()
            .map(|x| x.to_char().to_ascii_uppercase() as u8)
            .collect();
        // Per-read index used for only a few candidate ASVs; see SHORT_REF_BUCKET_BITS.
        let mut builder = Aligner::builder().lrhq();
        builder.idxopt.bucket_bits = SHORT_REF_BUCKET_BITS;
        let aligner = builder
            .with_cigar()
            .with_seq(&seq_u8)
            .expect("Failed to create per-sample aligner");

        let mut best_alns: Vec<(usize, i32)> = Vec::new();
        for (asv_idx, _) in best_asv_indices {
            let aln = aligner
                .map(&self.seqs_u8[asv_idx], true, false, None, None, None)
                .unwrap();
            if aln.is_empty() {
                continue;
            }
            best_alns.push((asv_idx, aln[0].alignment.as_ref().unwrap().nm));
        }
        let best_nm = best_alns.iter().map(|x| x.1).min()?;
        let mut asv_indices: Vec<usize> = best_alns
            .iter()
            .filter(|(_, nm)| *nm == best_nm)
            .map(|(idx, _)| *idx)
            .collect();
        asv_indices.sort();
        Some(EquivalenceClass { asv_indices })
    }
}

/// What a Stage 7 result set becomes on the final (post-EM) ASVs.
enum SetOutcome {
    Class(EquivalenceClass),
    Filtered,
    /// Every ASV in the set was dropped after EM; the read must be mapped again.
    Remap,
}

/// Per-sample ASV depths: EM on each sample's equivalence classes against the final ASVs.
///
/// Stage 7 already mapped every read against the same ASV sequences; the final ASVs are
/// those with nonzero EM depth. Each read's 7b result is therefore its Stage 7
/// lowest-mismatch set restricted to the surviving ASVs (with the same alignment NMs), and
/// only reads whose whole set was dropped need to be mapped again. Without `stage7_hits`,
/// every read is mapped.
pub fn compute_per_sample_depths(
    twin_reads: &[TwinRead],
    n_samples: usize,
    consensuses: &[ConsensusSequence],
    kmer_info: &KmerGlobalInfo,
    args: &Cli,
    asv_fasta_path: &Path,
    stage7_hits: Option<&ReadAsvHits>,
) -> Vec<Vec<usize>> {
    let n_asvs = consensuses.len();
    let mut result = vec![vec![0usize; n_samples]; n_asvs];

    if n_asvs == 0 || n_samples == 0 {
        return result;
    }

    if args.low_polymorphism {
        return compute_per_sample_depths_minimap2(
            twin_reads,
            n_samples,
            consensuses,
            args,
            asv_fasta_path,
        );
    }

    let outcomes: Vec<SetOutcome> = stage7_hits.map_or_else(Vec::new, |hits| {
        let final_index: FxHashMap<usize, usize> =
            consensuses.iter().enumerate().map(|(i, c)| (c.id, i)).collect();
        hits.sets
            .iter()
            .map(|set| {
                let surviving: Vec<(usize, u32)> = set
                    .iter()
                    .filter_map(|&(asv, rank)| {
                        final_index.get(&hits.asv_ids[asv as usize]).map(|&i| (i, rank))
                    })
                    .collect();
                if surviving.is_empty() {
                    return SetOutcome::Remap;
                }
                match surviving.iter().map(|x| x.1).filter(|&r| r != UNALIGNED).min() {
                    None => SetOutcome::Filtered,
                    Some(best) => {
                        let mut asv_indices: Vec<usize> =
                            surviving.iter().filter(|x| x.1 == best).map(|x| x.0).collect();
                        asv_indices.sort();
                        SetOutcome::Class(EquivalenceClass { asv_indices })
                    }
                }
            })
            .collect()
    });
    let mapper: std::sync::OnceLock<SampleMapper> = std::sync::OnceLock::new();
    let mut remapped_reads = 0usize;

    for sample_k in 0..n_samples {
        let sample_k_u32 = sample_k as u32;
        let reads: Vec<usize> = (0..twin_reads.len())
            .filter(|&r| twin_reads[r].file_idx == sample_k_u32)
            .collect();
        if reads.is_empty() {
            log::info!("Sample {}: no reads found", sample_k);
            continue;
        }

        let mut eq_classes: HashMap<EquivalenceClass, usize> = HashMap::new();
        let mut total_assigned = 0usize;
        let mut to_map: Vec<usize> = Vec::new();
        match stage7_hits {
            Some(hits) => {
                for &r in &reads {
                    let set = hits.read_set[r];
                    if set == NO_HITS {
                        continue; // no ASV passed in Stage 7, so none of the survivors can
                    }
                    match &outcomes[set as usize] {
                        SetOutcome::Class(class) => {
                            *eq_classes.entry(class.clone()).or_insert(0) += 1;
                            total_assigned += 1;
                        }
                        SetOutcome::Filtered => {}
                        SetOutcome::Remap => to_map.push(r),
                    }
                }
            }
            None => to_map = reads.clone(),
        }
        if !to_map.is_empty() {
            remapped_reads += to_map.len();
            let mapper = mapper.get_or_init(|| SampleMapper::new(kmer_info, args, asv_fasta_path));
            let classes: Vec<Option<EquivalenceClass>> = to_map
                .par_iter()
                .map(|&r| mapper.class_of(&twin_reads[r], args))
                .collect();
            for class in classes.into_iter().flatten() {
                *eq_classes.entry(class).or_insert(0) += 1;
                total_assigned += 1;
            }
        }

        let eq_classes = sorted_equivalence_classes(eq_classes);
        let filtered = reads.len() - total_assigned;

        log::info!(
            "Sample {}: {} reads assigned, {} filtered",
            sample_k,
            total_assigned,
            filtered
        );

        if eq_classes.is_empty() || total_assigned == 0 {
            continue;
        }

        // EM
        let mut asv_abundances = vec![1.0 / n_asvs as f64; n_asvs];
        let convergence_threshold = 0.01 / total_assigned as f64;
        let mut iteration = 0;
        const MAX_ITERATIONS: usize = 10000;
        loop {
            iteration += 1;
            let mut new_abundances = vec![0.0f64; n_asvs];
            for (eq_class, count) in &eq_classes {
                let denom: f64 = eq_class
                    .asv_indices
                    .iter()
                    .map(|&i| asv_abundances[i])
                    .sum();
                if denom > 0.0 {
                    for &asv_idx in &eq_class.asv_indices {
                        new_abundances[asv_idx] +=
                            (*count as f64) * asv_abundances[asv_idx] / denom;
                    }
                }
            }
            let total: f64 = new_abundances.iter().sum();
            if total > 0.0 {
                for a in new_abundances.iter_mut() {
                    *a /= total_assigned as f64;
                }
            }
            let max_change = asv_abundances
                .iter()
                .zip(new_abundances.iter())
                .map(|(o, n)| (o - n).abs())
                .fold(0.0, f64::max);
            asv_abundances = new_abundances;
            if max_change < convergence_threshold || iteration >= MAX_ITERATIONS {
                break;
            }
        }

        for asv_idx in 0..n_asvs {
            result[asv_idx][sample_k] =
                (asv_abundances[asv_idx] * total_assigned as f64).round() as usize;
        }
    }
    if stage7_hits.is_some() {
        log::debug!("Stage 7b re-mapped {} reads whose Stage 7 ASVs were all dropped", remapped_reads);
    }

    result
}

fn compute_per_sample_depths_minimap2(
    twin_reads: &[TwinRead],
    n_samples: usize,
    consensuses: &[ConsensusSequence],
    args: &Cli,
    asv_fasta_path: &Path,
) -> Vec<Vec<usize>> {
    let n_asvs = consensuses.len();
    let mut result = vec![vec![0usize; n_samples]; n_asvs];

    let aligner = Aligner::builder()
        .lrhq()
        .with_index_threads(args.threads)
        .with_cigar()
        .with_index(asv_fasta_path.to_str().unwrap(), None)
        .expect("Failed to build minimap2 index for per-sample depths");

    for sample_k in 0..n_samples {
        let sample_k_u32 = sample_k as u32;

        let eq_classes: Mutex<HashMap<EquivalenceClass, usize>> = Mutex::new(HashMap::new());
        let filtered_count = Mutex::new(0usize);
        let total_assigned = Mutex::new(0usize);

        twin_reads
            .par_iter()
            .filter(|r| r.file_idx == sample_k_u32)
            .for_each(|twin_read| {
                let seq: Vec<u8> = twin_read
                    .dna_seq
                    .iter()
                    .map(|x| x.to_char().to_ascii_uppercase() as u8)
                    .collect();
                let mappings = aligner
                    .map(&seq, true, false, None, None, None)
                    .unwrap_or_default();

                let valid: Vec<_> = mappings
                    .iter()
                    .filter(|m| m.mapq > 0 && m.alignment.is_some())
                    .collect();

                if valid.is_empty() {
                    *filtered_count.lock().unwrap() += 1;
                    return;
                }

                let best_nm = valid
                    .iter()
                    .map(|m| m.alignment.as_ref().unwrap().nm)
                    .min()
                    .unwrap();
                let mut best_asv_indices: Vec<usize> = valid
                    .iter()
                    .filter(|m| m.alignment.as_ref().unwrap().nm == best_nm)
                    .map(|m| m.target_id as usize)
                    .collect();
                best_asv_indices.sort();
                best_asv_indices.dedup();

                let eq_class = EquivalenceClass {
                    asv_indices: best_asv_indices,
                };
                *eq_classes.lock().unwrap().entry(eq_class).or_insert(0) += 1;
                *total_assigned.lock().unwrap() += 1;
            });

        let eq_classes = sorted_equivalence_classes(eq_classes.into_inner().unwrap());
        let total_assigned = total_assigned.into_inner().unwrap();
        let filtered = filtered_count.into_inner().unwrap();

        log::info!(
            "Sample {} (minimap2): {} reads assigned, {} filtered",
            sample_k,
            total_assigned,
            filtered
        );

        if eq_classes.is_empty() || total_assigned == 0 {
            continue;
        }

        let mut asv_abundances = vec![1.0 / n_asvs as f64; n_asvs];
        let convergence_threshold = 0.01 / total_assigned as f64;
        let mut iteration = 0;
        const MAX_ITERATIONS: usize = 10000;
        loop {
            iteration += 1;
            let mut new_abundances = vec![0.0f64; n_asvs];
            for (eq_class, count) in &eq_classes {
                let denom: f64 = eq_class
                    .asv_indices
                    .iter()
                    .map(|&i| asv_abundances[i])
                    .sum();
                if denom > 0.0 {
                    for &asv_idx in &eq_class.asv_indices {
                        new_abundances[asv_idx] +=
                            (*count as f64) * asv_abundances[asv_idx] / denom;
                    }
                }
            }
            let total: f64 = new_abundances.iter().sum();
            if total > 0.0 {
                for a in new_abundances.iter_mut() {
                    *a /= total_assigned as f64;
                }
            }
            let max_change = asv_abundances
                .iter()
                .zip(new_abundances.iter())
                .map(|(o, n)| (o - n).abs())
                .fold(0.0, f64::max);
            asv_abundances = new_abundances;
            if max_change < convergence_threshold || iteration >= MAX_ITERATIONS {
                break;
            }
        }

        for asv_idx in 0..n_asvs {
            result[asv_idx][sample_k] =
                (asv_abundances[asv_idx] * total_assigned as f64).round() as usize;
        }
    }

    result
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
