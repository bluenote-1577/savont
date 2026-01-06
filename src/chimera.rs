use crate::types::*;
use crate::cli::ClusterArgs as Cli;
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
    best_left_refs: Vec<usize>,
    best_left_lens: Vec<usize>,
    best_right_refs: Vec<usize>,
    best_right_lens: Vec<usize>,
}

/// Detect chimeric consensus sequences
/// A consensus is chimeric if:
/// 1. It has lower depth than its potential parents
/// 2. Its left portion matches one parent perfectly (excluding indels)
/// 3. Its right portion matches a different parent perfectly (excluding indels)
/// 4. The two parents are <99% similar to each other
/// 5. Combined, the matches cover >=95% of the query
pub fn detect_chimeras(
    consensuses: &mut [ConsensusSequence],
    args: &Cli,
) -> Vec<ChimeraInfo> {
    if consensuses.is_empty() {
        return Vec::new();
    }

    log::info!("Starting chimera detection for {} consensuses", consensuses.len());

    // Calculate pairwise similarities between all consensuses (for the <99% check)
    let similarities = calculate_pairwise_similarities(consensuses, args);
    let chimera_scores : Mutex<HashMap<usize, f64, _>> = Mutex::new(HashMap::new());

    let chimeras = Mutex::new(Vec::new());
    let min_match_length = args.chimera_detect_length.unwrap_or(args.min_read_length / 10);
    if min_match_length < 10{
        log::warn!("Chimera detection match length is set to a very low value of {} < 10. This may lead to false positives.", min_match_length);
    }

    // For each consensus, check if it's a chimera (using decompressed sequences)
    consensuses.par_iter().enumerate().for_each(|(query_idx, query_consensus)| {
        let query_seq = query_consensus.decompressed_sequence.as_ref()
            .expect("Consensus sequence must be decompressed before chimera detection");
        let query_depth = query_consensus.depth;
        let query_len = query_seq.len();

        // Find best left and right alignments
        let mut best_alignments = BestAlignments {
            best_left_refs: Vec::new(),
            best_left_lens: Vec::new(),
            best_right_refs: Vec::new(),
            best_right_lens: Vec::new(),
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
                                args,
                            );


                            if let Some(left_match) = left_match{
                                best_alignments.best_left_lens.push(left_match);
                                best_alignments.best_left_refs.push(ref_idx);
                            }
                            if let Some(right_match) = right_match{
                                best_alignments.best_right_lens.push(right_match);
                                best_alignments.best_right_refs.push(ref_idx);
                            }
                            // if query_idx == 57{
                                // println!("Query 57 vs Ref {}: Left match: {}, Right match: {}, CS: {}", ref_idx, left_match, right_match, alignment.cs.as_ref().unwrap());
                            // }
                        }
                    }
                }
            }
        }

        let mut min_chimera_score = 0.0 as f64;
        for (&left_ref, &left_len) in best_alignments.best_left_refs.iter().zip(best_alignments.best_left_lens.iter()) {
            let similarity_score = similarities.get(&(left_ref.min(query_idx), left_ref.max(query_idx)))
                .copied()
                .unwrap_or(1.0);
            if similarity_score < 0.85 && left_len < 500{
                continue;
            }
            let log_score = (similarity_score).ln() * left_len as f64;
            min_chimera_score = min_chimera_score.min(log_score);
        }
        for (&right_ref, &right_len) in best_alignments.best_right_refs.iter().zip(best_alignments.best_right_lens.iter()) {
            let similarity_score = similarities.get(&(right_ref.min(query_idx), right_ref.max(query_idx)))
                .copied()
                .unwrap_or(1.0);
            if similarity_score < 0.85 && right_len < 500{
                continue;
            }
            let log_score = (similarity_score).ln() * right_len as f64;
            min_chimera_score = min_chimera_score.min(log_score);
        }
        chimera_scores.lock().unwrap().insert(query_idx, min_chimera_score);

        // Check if this is a chimera
        for (&left_ref, &left_len) in best_alignments.best_left_refs.iter().zip(best_alignments.best_left_lens.iter()) {
            for (&right_ref, &right_len) in best_alignments.best_right_refs.iter().zip(best_alignments.best_right_lens.iter()) {
                log::trace!(
                    "Query {} alignment: Left ref {}, left len {}, Right ref {}, right len {}",
                    consensuses[query_idx].id, consensuses[left_ref].id, left_len,
                    consensuses[right_ref].id, right_len
                );
                // Must be two different parents
                if left_ref != right_ref {
                    // Check if parents are <99% similar
                    let parent_similarity = similarities.get(&(left_ref.min(right_ref), left_ref.max(right_ref)))
                        .copied()
                        .unwrap_or(0.0);

                    if parent_similarity < 0.97 || (parent_similarity < 0.995 && (consensuses[left_ref].depth > query_depth * 10) && consensuses[right_ref].depth > query_depth * 10) {
                        // Check if coverage is >=95%
                        let total_match = left_len + right_len;
                        let coverage_fraction = total_match as f64 / query_len as f64;

                        if coverage_fraction >= (0.9 * parent_similarity).min(0.8) && (coverage_fraction < 1.5 || (parent_similarity < 0.99 && coverage_fraction < 1.8)) {
                            log::debug!(
                                "Detected chimera: consensus {} (depth {}) = left_parent {} + right_parent {} (coverage: {:.2}%, parent similarity: {:.2}%)",
                                consensuses[query_idx].id, query_depth, consensuses[left_ref].id, consensuses[right_ref].id, coverage_fraction * 100.0, parent_similarity * 100.0
                            );

                            chimeras.lock().unwrap().push(ChimeraInfo {
                                query_idx,
                                left_parent_idx: left_ref,
                                right_parent_idx: right_ref,
                                left_match_len: left_len,
                                right_match_len: right_len,
                                query_len,
                                coverage_fraction,
                            });
                            break;
                        }
                        else{
                            log::trace!(
                                "Consensus {} failed coverage check: coverage {:.2}%, required {:.2}%. Parent similarity {:.2}%",
                                consensuses[query_idx].id, coverage_fraction * 100.0, 95.0, parent_similarity * 100.0
                            );
                        }
                    }
                    else{
                        log::trace!(
                            "Consensus {} failed parent similarity check: similarity {:.2}%, required <99%",
                            consensuses[query_idx].id, parent_similarity * 100.0
                        );
                    }
                }
            }
        }


        // Detection step 2: if a consensus has > 90% covered perfectly by any single parent, and > 20 mismatches
        // Call it as a chimera
        let match_ref = best_alignments.best_left_refs.iter().zip(best_alignments.best_left_lens.iter())
            .chain(best_alignments.best_right_refs.iter().zip(best_alignments.best_right_lens.iter()));

        for (reference, match_len) in match_ref {
            if *match_len >= (query_len - min_match_length) {
                let sim = similarities.get(&(*reference.min(&query_idx), *reference.max(&query_idx)))
                    .copied()
                    .unwrap_or(1.0);
                let total_mismatches = ((1.0 - sim) * query_len as f64) as usize;
                let ratio_depth = consensuses[*reference].depth as f64 / query_depth as f64;
                if ratio_depth < 3.0{
                    continue;
                }
                if total_mismatches as f64 > 20.0 / ratio_depth.log2(){
                    //Chimera
                    chimeras.lock().unwrap().push(ChimeraInfo {
                        query_idx,
                        left_parent_idx: *reference,
                        right_parent_idx: *reference,
                        left_match_len: *match_len,
                        right_match_len: 0,
                        query_len,
                        coverage_fraction: (*match_len as f64) / (query_len as f64),
                    });
                    log::debug!("Detected chimera by single-parent match: consensus {} (depth {}) = parent {} (match length {}, mismatches {})",
                        consensuses[query_idx].id, query_depth, consensuses[*reference].id, *match_len, total_mismatches);
                }
            }
        }
    });



    for i in 0..consensuses.len() {
        if !chimera_scores.lock().unwrap().contains_key(&i) {
            chimera_scores.lock().unwrap().insert(i, 0.0);
        }
        else{
            log::debug!("Consensus {} chimera score: {:.4}", consensuses[i].id, chimera_scores.lock().unwrap()[&i]);
            consensuses[i].chimera_score = Some(chimera_scores.lock().unwrap()[&i] as i64);
        }
    }

    let chimeras = chimeras.into_inner().unwrap();
    log::debug!("Detected {} chimeric combinations", chimeras.len());

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
    args: &Cli,
) -> (Option<usize>, Option<usize>) {
    
    // Track match positions in the query
    let mut left_max_perfect = 0;
    let mut right_max_perfect = 0;
    let pcr_slack = 15;

    // Process CIGAR to find matching positions
    {
        let mut num_errs = 0;
        let mut query_pos = query_start;
        let mut target_pos = target_start;

        for &(length, op) in cigar {
            if num_errs > args.chimera_allowable_errors{
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
                                // Elevated error rates near the edges: PCR primer mismatches and polishing issues...
                                num_errs += 1;
                                if num_errs > args.chimera_allowable_errors && query_pos + i >= pcr_slack{
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
            if num_errs > args.chimera_allowable_errors{
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
                            if num_errs > args.chimera_allowable_errors && query_pos_right - i + pcr_slack <= query_seq.len(){
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

    let mut right_max_perfect_opt = Some(right_max_perfect);
    let mut left_max_perfect_opt = Some(left_max_perfect);
    let min_match_length = args.chimera_detect_length.unwrap_or((args.min_read_length / 10).max(100));

    if right_max_perfect < min_match_length || left_max_perfect >= right_max_perfect {
        right_max_perfect_opt = None;
    }

    if left_max_perfect < min_match_length || right_max_perfect >= left_max_perfect {
        left_max_perfect_opt = None;
    }

    if rc{
        (right_max_perfect_opt, left_max_perfect_opt)
    }
    else{
        (left_max_perfect_opt, right_max_perfect_opt)
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
                .lrhq()
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

    let average_chimera_score: f64 = chimera_indices.iter()
        .filter_map(|&idx| consensuses.get(idx).and_then(|c| c.chimera_score.map(|s| s as f64)))
        .sum::<f64>() / chimera_indices.len() as f64;

    let std_chimera_score: f64 = (chimera_indices.iter()
        .filter_map(|&idx| consensuses.get(idx).and_then(|c| c.chimera_score.map(|s| s as f64)))
        .map(|score| (score - average_chimera_score).powi(2))
        .sum::<f64>() / chimera_indices.len() as f64).sqrt();

    let filtered: Vec<ConsensusSequence> = consensuses.into_iter()
        .enumerate()
        .filter(|(idx, _)| !chimera_indices.contains(idx))
        .map(|(_, cons)| cons)
        .collect();

    log::info!("Filtered {} chimeric consensuses with mean chimera score {:.2} and std {:.2}, {} remaining",
        original_count - filtered.len(), average_chimera_score, std_chimera_score, filtered.len());

    filtered
}
