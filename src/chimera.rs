use crate::types::*;
use crate::cli::Cli;
use crate::utils;
use minimap2::{Aligner, Strand};
use rayon::prelude::*;
use std::sync::Mutex;
use std::collections::HashMap;

/// Represents a chimeric consensus sequence
#[derive(Debug, Clone)]
pub struct ChimeraInfo {
    pub query_idx: usize,
    pub left_parent_idx: usize,
    pub right_parent_idx: usize,
    pub left_match_len: usize,
    pub right_match_len: usize,
    pub query_len: usize,
    pub coverage_fraction: f64,
}

/// Information about the best left and right alignments for a query
#[derive(Debug, Clone)]
struct BestAlignments {
    best_left_ref: Option<usize>,
    best_left_len: usize,
    best_right_ref: Option<usize>,
    best_right_len: usize,
}

/// Detect chimeric consensus sequences
/// A consensus is chimeric if:
/// 1. It has lower depth than its potential parents
/// 2. Its left portion matches one parent perfectly (excluding indels)
/// 3. Its right portion matches a different parent perfectly (excluding indels)
/// 4. The two parents are <99% similar to each other
/// 5. Combined, the matches cover >=95% of the query
pub fn detect_chimeras(
    consensuses: &[ConsensusSequence],
    args: &Cli,
) -> Vec<ChimeraInfo> {
    if consensuses.is_empty() {
        return Vec::new();
    }

    log::info!("Starting chimera detection for {} consensuses", consensuses.len());

    // Calculate pairwise similarities between all consensuses (for the <99% check)
    let similarities = calculate_pairwise_similarities(consensuses, args);

    let chimeras = Mutex::new(Vec::new());

    // For each consensus, check if it's a chimera (using decompressed sequences)
    consensuses.par_iter().enumerate().for_each(|(query_idx, query_consensus)| {
        let query_seq = query_consensus.decompressed_sequence.as_ref()
            .expect("Consensus sequence must be decompressed before chimera detection");
        let query_depth = query_consensus.depth;
        let query_len = query_seq.len();

        // Find best left and right alignments
        let mut best_alignments = BestAlignments {
            best_left_ref: None,
            best_left_len: 0,
            best_right_ref: None,
            best_right_len: 0,
        };

        // Align query to each potential parent
        for (ref_idx, ref_consensus) in consensuses.iter().enumerate() {
            if ref_idx == query_idx {
                continue;
            }

            // Only consider higher-depth consensuses as parents
            if ref_consensus.depth <= query_depth * 3 {
                continue;
            }

            // Use decompressed sequence for alignment
            let ref_seq = ref_consensus.decompressed_sequence.as_ref()
                .expect("Consensus sequence must be decompressed before chimera detection");

            // Align query to this reference
            let aligner = Aligner::builder()
                .map_ont()
                .with_cigar()
                .with_seq(ref_seq)
                .expect("Failed to create aligner");

            if let Ok(mappings) = aligner.map(query_seq, true, false, None, None, None) {
                for mapping in mappings.iter() {
                    if let Some(ref alignment) = mapping.alignment {
                        if let Some(ref cigar) = alignment.cigar {
                            // Handle reverse complement and adjust coordinates
                            let rc;
                            let (final_query_seq, query_start, query_end) = if mapping.strand == Strand::Reverse {
                                let qstart = query_len - mapping.query_end as usize;
                                let qend = query_len - mapping.query_start as usize;
                                rc = true;
                                (utils::reverse_complement(query_seq), qstart, qend)
                            } else {
                                rc = false;
                                (query_seq.clone(), mapping.query_start as usize, mapping.query_end as usize)
                            };

                            // Calculate perfect match lengths from both sides (on decompressed sequences)
                            let (left_match, right_match) = calculate_match_lengths(
                                cigar,
                                &final_query_seq,
                                ref_seq,
                                query_start,
                                query_end,
                                mapping.target_start as usize,
                                mapping.target_end as usize,
                                rc,
                            );

                            // Update best left alignment
                            if left_match > best_alignments.best_left_len {
                                best_alignments.best_left_len = left_match;
                                best_alignments.best_left_ref = Some(ref_idx);
                            }

                            // Update best right alignment
                            if right_match > best_alignments.best_right_len {
                                best_alignments.best_right_len = right_match;
                                best_alignments.best_right_ref = Some(ref_idx);
                            }

                            // if query_idx == 57{
                                // println!("Query 57 vs Ref {}: Left match: {}, Right match: {}, CS: {}", ref_idx, left_match, right_match, alignment.cs.as_ref().unwrap());
                            // }
                        }
                    }
                }
            }
        }

        // Check if this is a chimera
        if let (Some(left_ref), Some(right_ref)) = (best_alignments.best_left_ref, best_alignments.best_right_ref) {
            log::debug!(
                "Query {} best left ref: {:?}, best right ref: {:?}, Left: {}, Right: {}",
                consensuses[query_idx].id, best_alignments.best_left_ref, best_alignments.best_right_ref, best_alignments.best_left_len, best_alignments.best_right_len
            );
            // Must be two different parents
            if left_ref != right_ref {
                // Check if parents are <99% similar
                let parent_similarity = similarities.get(&(left_ref.min(right_ref), left_ref.max(right_ref)))
                    .copied()
                    .unwrap_or(0.0);

                if parent_similarity < 0.97 || (parent_similarity < 99.5 && (consensuses[left_ref].depth > query_depth * 10) && consensuses[right_ref].depth > query_depth * 10) {
                    // Check if coverage is >=95%
                    let total_match = best_alignments.best_left_len + best_alignments.best_right_len;
                    let coverage_fraction = total_match as f64 / query_len as f64;

                    if coverage_fraction >= 0.9 * parent_similarity && (coverage_fraction < 1.5 || (parent_similarity < 0.99 && coverage_fraction < 1.8)) {
                        log::info!(
                            "Detected chimera: consensus {} (depth {}) = left_parent {} + right_parent {} (coverage: {:.2}%, parent similarity: {:.2}%)",
                            consensuses[query_idx].id, query_depth, left_ref, right_ref, coverage_fraction * 100.0, parent_similarity * 100.0
                        );

                        chimeras.lock().unwrap().push(ChimeraInfo {
                            query_idx,
                            left_parent_idx: left_ref,
                            right_parent_idx: right_ref,
                            left_match_len: best_alignments.best_left_len,
                            right_match_len: best_alignments.best_right_len,
                            query_len,
                            coverage_fraction,
                        });
                    }
                    else{
                        log::debug!(
                            "Consensus {} failed coverage check: coverage {:.2}%, required {:.2}%. Parent similarity {:.2}%",
                            consensuses[query_idx].id, coverage_fraction * 100.0, 95.0, parent_similarity * 100.0
                        );
                    }
                }
                else{
                    log::debug!(
                        "Consensus {} failed parent similarity check: similarity {:.2}%, required <99%",
                        consensuses[query_idx].id, parent_similarity * 100.0
                    );
                }
            }
        }
    });

    let chimeras = chimeras.into_inner().unwrap();
    log::info!("Detected {} chimeric consensuses", chimeras.len());

    chimeras
}

/// Calculate left and right perfect match lengths from CIGAR by actually comparing bases
/// Returns (left_match_length, right_match_length)
/// Op codes: 0 = Match/Mismatch, 1 = Insertion, 2 = Deletion, 100 = Mismatch, 101 = Soft clip, 102 = Hard clip
fn calculate_match_lengths(
    cigar: &[(u32, u8)],
    query_seq: &[u8],
    target_seq: &[u8],
    query_start: usize,
    query_end: usize,
    target_start: usize,
    target_end: usize,
    rc: bool,
) -> (usize, usize) {
    
    // Track match positions in the query
    let mut left_max_perfect = 0;
    let mut right_max_perfect = 0;

    // Process CIGAR to find matching positions
    {
        let mut num_errs = 0;
        let mut query_pos = query_start;
        let mut target_pos = target_start;

        for &(length, op) in cigar {
            if num_errs > 1{
                break;
            }
            let len = length as usize;
            match op {
                0 => {
                    // Match or mismatch - check actual bases
                    for i in 0..len {
                        if query_pos + i < query_seq.len() && target_pos + i < target_seq.len() {
                            if query_seq[query_pos + i] == target_seq[target_pos + i] {
                                left_max_perfect += 1;
                            }
                            else{
                                num_errs += 1;
                                if num_errs > 1{
                                    break;
                                }
                            }
                        }
                    }
                    query_pos += len;
                    target_pos += len;
                }
                1 => {
                    // Insertion in query - skip query bases
                    query_pos += len;
                }
                2 => {
                    // Deletion in query - skip target bases
                    target_pos += len;
                }
                _ => {
                    log::warn!("Unexpected CIGAR operation: {}", op);
                }
            }
        }
    }

    
    {
        let mut query_pos_right = query_end;
        let mut target_pos_right = target_end;
        let mut num_errs = 0;

        // Process CIGAR to find matching positions
        for &(length, op) in cigar.iter().rev() {
            if num_errs > 1{
                break;
            }
            let len = length as usize;
            match op {
                0 => {
                    // Match or mismatch - check actual bases NEED -1 because the intervals are [start,end)
                    for i in 0..len {
                        if query_seq[query_pos_right - i - 1] == target_seq[target_pos_right - i - 1] {
                            right_max_perfect += 1;
                        }
                        else{
                            num_errs += 1;
                            if num_errs > 1{
                                break;
                            }
                        }
                    }
                    query_pos_right -= len;
                    target_pos_right -= len;
                }
                1 => {
                    // Insertion in query - skip query bases
                    query_pos_right -= len;
                }
                2 => {
                    // Deletion in query - skip target bases
                    target_pos_right -= len;
                }
                _ => {
                    log::warn!("Unexpected CIGAR operation: {}", op);
                }
            }
        }
    }

    if right_max_perfect < 100 || left_max_perfect >= right_max_perfect {
        right_max_perfect = 0;
    }

    if left_max_perfect < 100 || right_max_perfect >= left_max_perfect {
        left_max_perfect = 0;
    }

    if rc{
        (right_max_perfect, left_max_perfect)
    }
    else{
        (left_max_perfect, right_max_perfect)
    }
}

/// Calculate pairwise similarities between all consensuses using NM / alignment_length
/// Returns a HashMap of (idx1, idx2) -> similarity
fn calculate_pairwise_similarities(
    consensuses: &[ConsensusSequence],
    _args: &Cli,
) -> HashMap<(usize, usize), f64> {
    let similarities = Mutex::new(HashMap::new());

    log::info!("Calculating pairwise similarities for {} consensuses (using decompressed sequences)", consensuses.len());

    consensuses.par_iter().enumerate().for_each(|(i, cons_i)| {
        let seq_i = cons_i.decompressed_sequence.as_ref()
            .expect("Consensus sequence must be decompressed before chimera detection");

        for (j, cons_j) in consensuses.iter().enumerate() {
            if i >= j {
                continue; // Only calculate once per pair
            }

            let seq_j = cons_j.decompressed_sequence.as_ref()
                .expect("Consensus sequence must be decompressed before chimera detection");

            // Align cons_i to cons_j (using decompressed sequences)
            let aligner = Aligner::builder()
                .map_ont()
                .with_cigar()
                .with_seq(seq_j)
                .expect("Failed to create aligner");

            if let Ok(mappings) = aligner.map(seq_i, false, false, None, None, None) {
                if let Some(best_mapping) = mappings.first() {
                    if let Some(ref alignment) = best_mapping.alignment {
                        // Calculate identity as 1 - (NM / alignment_length)
                        let alignment_len = (best_mapping.query_end - best_mapping.query_start) as f64;
                        let nm = alignment.nm as f64;

                        let identity = if alignment_len > 0.0 {
                            1.0 - (nm / alignment_len)
                        } else {
                            0.0
                        };

                        similarities.lock().unwrap().insert((i, j), identity);
                    }
                }
            }
        }
    });

    similarities.into_inner().unwrap()
}

/// Remove chimeric consensuses from the list
pub fn filter_chimeras(
    consensuses: Vec<ConsensusSequence>,
    chimeras: &[ChimeraInfo],
) -> Vec<ConsensusSequence> {
    let chimera_indices: std::collections::HashSet<usize> = chimeras.iter()
        .map(|c| c.query_idx)
        .collect();

    let original_count = consensuses.len();

    let filtered: Vec<ConsensusSequence> = consensuses.into_iter()
        .enumerate()
        .filter(|(idx, _)| !chimera_indices.contains(idx))
        .map(|(_, cons)| cons)
        .collect();

    log::info!("Filtered {} chimeric consensuses, {} remaining",
        original_count - filtered.len(), filtered.len());

    filtered
}
