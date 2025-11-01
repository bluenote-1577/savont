use crate::cli::Cli;
use std::sync::{Arc, Mutex};
use rayon::prelude::*;
use crate::types::*;
use fxhash::FxHashMap;
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;

/// Add a read to the inverted index
fn add_read_to_index(
    index: &mut FxHashMap<Kmer48, Vec<usize>>,
    read_id: usize,
    read_kmers: &[Kmer48],
) {
    for &kmer in read_kmers {
        index.entry(kmer).or_insert_with(Vec::new).push(read_id);
    }
}

/// Query a read against the index and calculate similarities
/// Returns Vec<(read_id, similarity)> for all candidates
fn query_read_against_index(
    index: &FxHashMap<Kmer48, Vec<usize>>,
    query_kmers: &[Kmer48],
    twin_reads: &[TwinRead],
    k: usize,
) -> Vec<(usize, f64)> {
    let mut candidates: FxHashMap<usize, usize> = FxHashMap::default();

    // Find all candidate reads via inverted index and count shared k-mers
    for kmer in query_kmers {
        if let Some(read_ids) = index.get(kmer) {
            for &read_id in read_ids {
                *candidates.entry(read_id).or_insert(0) += 1;
            }
        }
    }

    // Calculate similarity for each candidate
    let mut similarities = Vec::new();
    for (candidate_id, shared_count) in candidates {
        let candidate_kmer_count = twin_reads[candidate_id].minimizer_positions.len();
        let query_kmer_count = query_kmers.len();
        let max_count = candidate_kmer_count.max(query_kmer_count);

        if max_count == 0 {
            continue;
        }

        let ratio = shared_count as f64 / max_count as f64;
        let similarity = ratio.powf(1.0 / k as f64);

        similarities.push((candidate_id, similarity));
    }

    similarities.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    similarities
}

/// Function 1: Greedy sequential k-mer clustering
/// For each read: query against index. If no match > threshold, add to index.
/// Reads in the index are cluster representatives.
pub fn cluster_reads_by_kmers(
    twin_reads: &[TwinRead],
    args: &Cli,
    output_dir: &PathBuf,
) -> Vec<Vec<usize>> {
    log::info!("Starting greedy k-mer based clustering...");

    let k = args.kmer_size;
    let threshold = 0.95;

    // Inverted index: kmer -> Vec<read_id>
    let mut index: FxHashMap<Kmer48, Vec<usize>> = FxHashMap::default();

    // Cluster assignments: read_id -> representative_read_id
    let mut cluster_assignment: HashMap<usize, usize> = HashMap::new();

    // Representative reads (reads added to the index)
    let mut representatives: Vec<usize> = Vec::new();

    // Process reads sequentially
    for (read_id, read) in twin_reads.iter().enumerate() {
        let read_kmers = read.minimizer_kmers();

        // Query this read against the current index
        let similarities = query_read_against_index(&index, &read_kmers, twin_reads, k);

        // Check if any match exceeds threshold
        let best_match = similarities.first();

        if let Some(&(representative_id, similarity)) = best_match {
            if similarity > threshold {
                // Assign this read to the cluster of the representative
                cluster_assignment.insert(read_id, representative_id);
                continue; // Don't add this read to the index
            }
        }

        // No match found - this read becomes a new representative
        add_read_to_index(&mut index, read_id, &read_kmers);
        cluster_assignment.insert(read_id, read_id); // Represents itself
        representatives.push(read_id);

        if read_id % 10 == 0 && read_id > 0 {
            log::debug!(
                "Processed {} / {} reads. Current representatives: {}",
                read_id,
                twin_reads.len(),
                representatives.len()
            );
        }
    }

    log::info!(
        "Greedy clustering complete. {} cluster representatives found",
        representatives.len()
    );

    // Build clusters from assignments
    let mut clusters_map: HashMap<usize, Vec<usize>> = HashMap::new();
    for (read_id, rep_id) in cluster_assignment {
        clusters_map.entry(rep_id).or_insert_with(Vec::new).push(read_id);
    }

    let mut clusters: Vec<Vec<usize>> = clusters_map.into_values().collect();
    clusters.sort_by(|a, b| b.len().cmp(&a.len()));

    // Sort members within each cluster because lower IDs have better estimated accuracy
    for cluster in clusters.iter_mut(){
        cluster.sort_unstable();
    }

    // Remove small clusters
    clusters.retain(|cluster| cluster.len() >= args.min_cluster_size);

    // Write clusters to file
    let cluster_file = output_dir.join("kmer_clusters.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&cluster_file).unwrap());

    writeln!(writer, "cluster_id\tsize\trepresentative\tmembers").unwrap();

    for (cluster_id, cluster) in clusters.iter().enumerate() {
        let representative = cluster[0]; // First member is always the representative
        writeln!(
            writer,
            "cluster_{}\t{}\t{}\t{}",
            cluster_id,
            cluster.len(),
            representative,
            cluster.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
        ).unwrap();
    }

    log::info!(
        "K-mer clustering complete. {} clusters found. Largest cluster: {} reads",
        clusters.len(),
        clusters.first().map(|c| c.len()).unwrap_or(0)
    );
    log::info!("Wrote k-mer clusters to {}", cluster_file.display());

    clusters
}

/// Add a read's SNPmers to the SNPmer inverted index
fn add_read_snpmers_to_index(
    index: &mut FxHashMap<u64, Vec<(usize, Kmer48)>>,
    read_id: usize,
    read_snpmers: &[(u32, Kmer48)],
    k: usize,
) {
    let mask = !(3 << (k - 1));
    for &(_pos, kmer) in read_snpmers {
        let splitmer = kmer.to_u64() & mask;
        index.entry(splitmer).or_insert_with(Vec::new).push((read_id, kmer));
    }
}

/// Query a read's SNPmers against the SNPmer index
/// Returns (matches, mismatches) for each candidate read_id
/// Candidates with 0 mismatches are compatible
fn find_compatible_candidates(
    index: &FxHashMap<u64, Vec<(usize, Kmer48)>>,
    query_snpmers: &[(u32, Kmer48)],
    k: usize,
) -> FxHashMap<usize, (usize, usize)> {
    let mask = !(3 << (k - 1));

    // Track matches and mismatches for each candidate
    let mut candidate_stats: FxHashMap<usize, (usize, usize)> = FxHashMap::default();

    // Initialize all candidates that exist in the index
    for candidates in index.values() {
        for &(candidate_id, _) in candidates {
            candidate_stats.entry(candidate_id).or_insert((0, 0));
        }
    }

    // Check each query SNPmer
    for &(_pos, query_kmer) in query_snpmers {
        let query_splitmer = query_kmer.to_u64() & mask;

        if let Some(candidates) = index.get(&query_splitmer) {
            for &(candidate_id, candidate_kmer) in candidates {
                let stats = candidate_stats.get_mut(&candidate_id).unwrap();
                if query_kmer == candidate_kmer {
                    stats.0 += 1; // matches
                } else {
                    stats.1 += 1; // mismatches
                }
            }
        }
    }

    candidate_stats
}

/// Function 2: Greedy SNPmer clustering within each k-mer cluster
/// For each read in a k-mer cluster: query SNPmers against index.
/// Only add to index if NO SNPmer mismatches found.
pub fn cluster_reads_by_snpmers(
    twin_reads: &[TwinRead],
    kmer_clusters: &[Vec<usize>],
    args: &Cli,
    output_dir: &PathBuf,
) -> Vec<Vec<usize>> {
    log::info!("Starting greedy SNPmer-based clustering within k-mer clusters...");

    let k = args.kmer_size;

    // Shared data structures wrapped in Arc<Mutex<>> for thread safety
    let snpmer_cluster_assignment = Arc::new(Mutex::new(HashMap::new()));
    let local_clusters_map = Arc::new(Mutex::new(FxHashMap::default()));

    // Process each k-mer cluster independently in parallel
    kmer_clusters.par_iter().enumerate().for_each(|(kmer_cluster_id, kmer_cluster)| {
         // Skip empty k-mer clusters
        if kmer_cluster.len() < 1 {
            return;
        }

        // SNPmer inverted index for this k-mer cluster: splitmer -> Vec<(read_id, full_kmer)>
        let mut snpmer_index: FxHashMap<u64, Vec<(usize, Kmer48)>> = FxHashMap::default();

        // Local cluster assignments within this k-mer cluster
        let mut local_assignment: HashMap<usize, usize> = HashMap::new();

        let mut rep_size: HashMap<usize, usize> = HashMap::new();

        // Process reads in this k-mer cluster sequentially
        for &read_id in kmer_cluster {
            let read_snpmers = twin_reads[read_id].snpmers_vec();

            // Query this read's SNPmers against the current index to find compatible candidates
            let candidate_stats = find_compatible_candidates(&snpmer_index, &read_snpmers, k);

            // Filter to only compatible candidates (0 mismatches, >0 matches)
            let compatible_candidates: Vec<(usize, usize, usize)> = candidate_stats
                .iter()
                .filter(|(_, (matches, mismatches))| *mismatches == 0 && *matches > 0)
                .map(|(candidate_id, (matches, _))| (*candidate_id, *matches, rep_size[candidate_id]))
                .collect();

            if compatible_candidates.is_empty() {
                // No compatible candidates - this read becomes a new SNPmer representative
                add_read_snpmers_to_index(&mut snpmer_index, read_id, &read_snpmers, k);
                local_assignment.insert(read_id, read_id);
                rep_size.insert(read_id, 1);
            } else {
                // Choose the best compatible candidate
                // Sort by: fewest members (smallest cluster), then most matches
                let mut candidates_sorted: Vec<_> = compatible_candidates
                    .iter()
                    .map(|(cand_id, matches, size)| (*size, -(*matches as i64), *cand_id))
                    .collect();
                candidates_sorted.sort();

                let best_rep = candidates_sorted[0].2;
                local_assignment.insert(read_id, best_rep);
                *rep_size.entry(best_rep).or_insert(0) += 1;

                log::trace!(
                    "Assigned read {} to representative {} ({} matches)",
                    read_id, best_rep, candidate_stats[&best_rep].0
                );
            }
        }

        // Build SNPmer subclusters within this k-mer cluster
        let mut cluster_map = HashMap::new();
        for (read_id, rep_id) in &local_assignment {
            if read_id != rep_id {
                continue;
            }
            cluster_map.entry(*rep_id).or_insert_with(Vec::new).push(*rep_id);
        }
        for (read_id, rep_id) in local_assignment.iter() {
            if read_id == rep_id {
                continue;
            }
            cluster_map.entry(*rep_id).or_insert_with(Vec::new).push(*read_id);
        }

        let mut local_clusters: Vec<Vec<usize>> = cluster_map.into_values().collect();
        local_clusters.sort_by(|a, b| b.len().cmp(&a.len()));

        local_clusters.retain(|cluster| cluster.len() >= args.min_cluster_size);

        // Update shared data structures
        {
            let mut assignment = snpmer_cluster_assignment.lock().unwrap();
            for (read_id, rep_id) in local_assignment.iter() {
                assignment.insert(*read_id, *rep_id);
            }
        }

        {
            let mut clusters_map = local_clusters_map.lock().unwrap();
            clusters_map.entry(kmer_cluster_id).or_insert_with(Vec::new).extend(local_clusters);
        }


        if kmer_cluster_id % 100 == 0 && kmer_cluster_id > 0 {
            log::debug!(
                "Processed {} / {} k-mer clusters for SNPmer clustering",
                kmer_cluster_id,
                kmer_clusters.len()
            );
        }
    });

    // Extract data from Arc after parallel processing
    let local_clusters_map = Arc::try_unwrap(local_clusters_map)
        .unwrap()
        .into_inner()
        .unwrap();

    // Write SNPmer clusters to TSV file
    let cluster_file = output_dir.join("snpmer_clusters.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&cluster_file).unwrap());
    writeln!(writer, "kmer_cluster_id\tsnpmer_cluster_id\tsize\trepresentative\tmembers").unwrap();

    for (kmer_cluster_id, snpmer_clusters) in local_clusters_map.iter() {
        for (local_snpmer_id, snpmer_cluster) in snpmer_clusters.iter().enumerate() {
            if snpmer_cluster.is_empty() {
                continue;
            }
            let representative = snpmer_cluster[0];
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                kmer_cluster_id,
                local_snpmer_id,
                snpmer_cluster.len(),
                representative,
                snpmer_cluster.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
            ).unwrap();
        }
    }

    log::info!("Wrote SNPmer clusters to {}", cluster_file.display());

    let recluster = true;
    let local_clusters_all: Vec<Vec<usize>>;
    if recluster {
        local_clusters_all = recluster_using_consensus_reps(
            local_clusters_map,
            twin_reads,
            args,
        );
    } else {
        // Flatten all local clusters into a single list
        let mut all_clusters: Vec<Vec<usize>> = Vec::new();
        for (_kmer_cluster_id, snpmer_clusters) in local_clusters_map {
            for snpmer_cluster in snpmer_clusters {
                if !snpmer_cluster.is_empty() {
                    all_clusters.push(snpmer_cluster);
                }
            }
        }

        // Sort by size descending
        all_clusters.sort_by(|a, b| b.len().cmp(&a.len()));
        local_clusters_all = all_clusters;
    }

    log::info!(
        "SNPmer clustering complete. {} total SNPmer clusters from {} k-mer clusters",
        local_clusters_all.len(),
        kmer_clusters.len()
    );

    let final_file = output_dir.join("final_clusters.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&final_file).unwrap());

    for (snpmer_rep_id, cluster) in local_clusters_all.iter().enumerate() {
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

    return local_clusters_all;
}


/// Build consensus SNPmer representative for a cluster
fn build_consensus_snpmers(
    cluster: &[usize],
    twin_reads: &[TwinRead],
    k: usize,
) -> Vec<ConsensusSnpmer> {
    let mask = !(3 << (k - 1));

    // Map: splitmer -> Map: full_kmer -> (count, Vec<positions>)
    let mut splitmer_data: FxHashMap<u64, FxHashMap<Kmer48, (usize, Vec<u32>)>> = FxHashMap::default();

    // Step 1: Collect all SNPmers from all reads in the cluster with positions
    for &read_id in cluster {
        let snpmers = twin_reads[read_id].snpmers_vec();
        for &(pos, kmer) in &snpmers {
            let splitmer = kmer.to_u64() & mask;
            let entry = splitmer_data
                .entry(splitmer)
                .or_insert_with(FxHashMap::default)
                .entry(kmer)
                .or_insert((0, Vec::new()));
            entry.0 += 1;
            entry.1.push(pos);
        }
    }

    // Step 2: For each splitmer, find the most common full k-mer and compute median position
    let mut consensus_snpmers = Vec::new();
    for (splitmer, kmer_data) in splitmer_data {
        // Find the k-mer with maximum count
        if let Some((&best_kmer, (count, positions))) = kmer_data.iter().max_by_key(|(_, (count, _))| count) {
            if *count > (cluster.len()/10).max(1) {
                // Calculate median position
                let mut pos_sorted = positions.clone();
                pos_sorted.sort_unstable();
                let median_pos = if !pos_sorted.is_empty() {
                    pos_sorted[pos_sorted.len() / 2]
                } else {
                    0
                };
                consensus_snpmers.push(ConsensusSnpmer::new(median_pos, splitmer, best_kmer, *count as u32));
            }
        }
    }

    consensus_snpmers.sort_by_key(|cs| (cs.position, cs.splitmer));
    consensus_snpmers
}

/// Compare two consensus SNPmer sets and count matches and mismatches
/// Returns (matches, mismatches)
fn compare_consensus_snpmers(
    consensus1: &[ConsensusSnpmer],
    consensus2: &[ConsensusSnpmer],
) -> (usize, usize) {
    // Build index of consensus2's splitmer -> ConsensusSnpmer
    let mut splitmer_to_snpmer: FxHashMap<u64, &ConsensusSnpmer> = FxHashMap::default();
    for cs in consensus2 {
        splitmer_to_snpmer.insert(cs.splitmer, cs);
    }

    let mut matches = 0;
    let mut mismatches = 0;

    // Check each SNPmer in consensus1
    for cs1 in consensus1 {
        // If consensus2 has this splitmer position, check if k-mers match
        if let Some(cs2) = splitmer_to_snpmer.get(&cs1.splitmer) {
            if cs1.kmer == cs2.kmer {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
    }

    (matches, mismatches)
}

/// Check if two consensus SNPmer sets are concordant (no mismatches)
fn are_consensus_concordant(
    consensus1: &[ConsensusSnpmer],
    consensus2: &[ConsensusSnpmer],
) -> bool {
    let (matches, mismatches) = compare_consensus_snpmers(consensus1, consensus2);
    mismatches == 0 && matches >= consensus1.len().min(consensus2.len())
}

/// Reassign reads to their best matching cluster based on SNPmer comparison
/// Returns (reassigned_clusters, num_reassignments)
fn reassign_reads_to_best_cluster(
    current_clusters: Vec<Vec<usize>>,
    twin_reads: &[TwinRead],
    k: usize,
) -> (Vec<Vec<usize>>, usize) {
    // Build consensus for each cluster
    let cluster_consensus: Vec<Vec<ConsensusSnpmer>> = current_clusters
        .iter()
        .map(|cluster| build_consensus_snpmers(cluster, twin_reads, k))
        .collect();

    // Pre-build lookup maps for each cluster
    let splitmer_to_consensus_kmers: Vec<FxHashMap<u64, Kmer48>> = cluster_consensus
        .iter()
        .map(|consensus| {
            let mut map: FxHashMap<u64, Kmer48> = FxHashMap::default();
            for cs in consensus {
                map.insert(cs.splitmer, cs.kmer);
            }
            map
        })
        .collect();

    let mask = !(3 << (k - 1));

    // Create new cluster assignments
    let new_clusters: Mutex<Vec<Vec<usize>>> = Mutex::new(vec![Vec::new(); current_clusters.len()]);
    let num_reassignments = Mutex::new(0);

    // For each cluster, check each read
    //for (cluster_idx, cluster) in current_clusters.iter().enumerate() {
    current_clusters.par_iter().enumerate().for_each(|(cluster_idx, cluster)| {
        for &read_id in cluster {
            let read_snpmers = twin_reads[read_id].snpmers_vec();

            // Find the best matching cluster for this read
            let mut best_cluster = cluster_idx;
            let mut best_score = (usize::MAX, 0); // (mismatches, matches) - lower mismatches and higher matches is better

            for (candidate_idx, _) in cluster_consensus.iter().enumerate() {
                // Build index of consensus splitmer -> kmer
                let splitmer_to_kmer = &splitmer_to_consensus_kmers[candidate_idx];

                let mut matches = 0;
                let mut mismatches = 0;

                // Compare read's SNPmers against this consensus
                for &(_pos, kmer) in &read_snpmers {
                    let splitmer = kmer.to_u64() & mask;
                    if let Some(&consensus_kmer) = splitmer_to_kmer.get(&splitmer) {
                        if kmer == consensus_kmer {
                            matches += 1;
                        } else {
                            mismatches += 1;
                        }
                    }
                }

                // Update best if this is better (fewer mismatches, or same mismatches but more matches)
                let current_score = (mismatches, matches);
                if mismatches < best_score.0 || (mismatches == best_score.0 && matches > best_score.1) {
                    best_score = current_score;
                    best_cluster = candidate_idx;
                }
            }

            // Assign read to best cluster
            new_clusters.lock().unwrap()[best_cluster].push(read_id);

            if best_cluster != cluster_idx {
                *num_reassignments.lock().unwrap() += 1;
                log::trace!(
                    "Reassigned read {} from cluster {} to cluster {} (mismatches: {}, matches: {})",
                    read_id, cluster_idx, best_cluster, best_score.0, best_score.1
                );
            }
        }
    });

    // Remove empty clusters
    let mut new_clusters = new_clusters.into_inner().unwrap();
    new_clusters.retain(|cluster| !cluster.is_empty());
    let num_reassignments = num_reassignments.into_inner().unwrap();

    log::debug!("Reassignment complete: {} reads reassigned to better clusters", num_reassignments);

    (new_clusters, num_reassignments)
}

/// Perform one round of reclustering based on consensus representatives
/// Returns (merged_clusters, num_merges_performed)
fn recluster_one_round(
    current_clusters: Vec<Vec<usize>>,
    twin_reads: &[TwinRead],
    k: usize,
) -> (Vec<Vec<usize>>, usize) {
    // Build consensus representatives for each cluster ONCE before the loop
    let mut all_clusters: Vec<(Vec<usize>, Vec<ConsensusSnpmer>)> = Vec::new();

    for cluster in current_clusters {
        if cluster.is_empty() {
            continue;
        }

        let consensus = build_consensus_snpmers(&cluster, twin_reads, k);
        all_clusters.push((cluster, consensus));
    }

    // Sort clusters by size (descending) so we merge smaller into larger
    all_clusters.sort_by(|a, b| b.0.len().cmp(&a.0.len()));

    // Merge clusters with concordant consensus representatives
    let mut cluster_merged: Vec<bool> = vec![false; all_clusters.len()];
    let mut cluster_needs_consensus_rebuild: Vec<bool> = vec![false; all_clusters.len()];
    let mut merged_clusters: Vec<Vec<usize>> = Vec::new();
    let mut num_merges = 0;

    for i in 0..all_clusters.len() {
        if cluster_merged[i] {
            continue;
        }

        // If this cluster was merged into in a previous iteration, rebuild its consensus
        if cluster_needs_consensus_rebuild[i] {
            all_clusters[i].1 = build_consensus_snpmers(&all_clusters[i].0, twin_reads, k);
            cluster_needs_consensus_rebuild[i] = false;
        }

        // Try to merge smaller clusters into this one
        for j in i + 1..all_clusters.len() {
            if cluster_merged[j] {
                continue;
            }

            // Check if consensus representatives are concordant (bidirectional)
            // Use the pre-built consensus from the beginning of the function
            let concordant = {
                let consensus_i = &all_clusters[i].1;
                let consensus_j = &all_clusters[j].1;
                are_consensus_concordant(consensus_i, consensus_j) &&
                are_consensus_concordant(consensus_j, consensus_i)
            };

            if concordant {
                // Merge smaller cluster into larger cluster
                let old_size = all_clusters[i].0.len();
                let cluster_j_size = all_clusters[j].0.len();

                // Clone cluster j's reads first to avoid borrow checker issues
                let cluster_j_reads = all_clusters[j].0.clone();

                // Extend cluster i with cluster j's reads
                all_clusters[i].0.extend(cluster_j_reads);

                // Mark that cluster i needs consensus rebuild (will be done before next comparison)
                cluster_needs_consensus_rebuild[i] = true;

                cluster_merged[j] = true;
                num_merges += 1;
                log::debug!(
                    "Merged cluster {} (size {}) into cluster {} (old size: {}, new size: {})",
                    j, cluster_j_size, i, old_size, all_clusters[i].0.len()
                );
            } else {
                log::trace!(
                    "Clusters {} and {} not concordant, not merging",
                    i, j
                );
            }
        }

        // Rebuild consensus one final time if any merges happened for cluster i
        if cluster_needs_consensus_rebuild[i] {
            all_clusters[i].1 = build_consensus_snpmers(&all_clusters[i].0, twin_reads, k);
        }

        merged_clusters.push(all_clusters[i].0.clone());
    }

    // Sort by size descending
    merged_clusters.sort_by(|a, b| b.len().cmp(&a.len()));

    (merged_clusters, num_merges)
}

pub fn recluster_using_consensus_reps(
    clusters: FxHashMap<usize, Vec<Vec<usize>>>,
    twin_reads: &[TwinRead],
    args: &Cli,
) -> Vec<Vec<usize>>{
    log::info!("Starting iterative reclustering using consensus representatives...");

    let k = args.kmer_size;

    // Preserve the hierarchical structure: FxHashMap<kmer_cluster_id, Vec<snpmer_clusters>>
    let mut current_clusters = clusters;

    let total_groups = current_clusters.len();
    let total_initial_clusters: usize = current_clusters.values().map(|v| v.len()).sum();
    log::info!("Starting with {} k-mer groups containing {} total SNPmer clusters", total_groups, total_initial_clusters);

    // Step 2: Iteratively recluster until convergence
    let mut iteration = 0;
    loop {
        iteration += 1;
        let total_merges = Mutex::new(0);
        let total_reassignments = Mutex::new(0);

        let new_clusters: Mutex<FxHashMap<usize, Vec<Vec<usize>>>> = Mutex::new(FxHashMap::default());


        // Process each k-mer group independently
        //for (kmer_cluster_id, snpmer_clusters) in current_clusters {
        current_clusters.into_par_iter().for_each(|(kmer_cluster_id, snpmer_clusters)| {
            log::debug!("Processing k-mer group {} with {} SNPmer clusters", kmer_cluster_id, snpmer_clusters.len());

            // Step 2a: Merge clusters within this group
            let (merged_clusters, num_merges) = recluster_one_round(
                snpmer_clusters,
                twin_reads,
                k
            );

            *total_merges.lock().unwrap() += num_merges;

            // Step 2b: Reassign reads to best matching clusters within this group
            let (reassigned_clusters, num_reassignments) = reassign_reads_to_best_cluster(
                merged_clusters,
                twin_reads,
                k
            );

            *total_reassignments.lock().unwrap() += num_reassignments;

            // Store the updated clusters for this k-mer group
            if !reassigned_clusters.is_empty() {
                new_clusters.lock().unwrap().insert(kmer_cluster_id, reassigned_clusters);
            }
        });

        let new_clusters = new_clusters.into_inner().unwrap();
        let total_merges = total_merges.into_inner().unwrap();
        let total_reassignments = total_reassignments.into_inner().unwrap();

        log::info!(
            "Iteration {}: {} total merges, {} total reassignments across {} k-mer groups",
            iteration,
            total_merges,
            total_reassignments,
            new_clusters.len()
        );

        current_clusters = new_clusters;

        // Check for convergence (no merges and no reassignments across all groups)
        if total_merges == 0 {
            log::info!("Convergence reached after {} iterations", iteration);
            break;
        }
    }

    // Step 3: Flatten the hierarchical structure for final output and debugging
    let mut final_clusters: Vec<Vec<usize>> = Vec::new();
    for (_kmer_cluster_id, snpmer_clusters) in &current_clusters {
        for cluster in snpmer_clusters {
            if !cluster.is_empty() {
                final_clusters.push(cluster.clone());
            }
        }
    }

    // Sort by size descending
    final_clusters.sort_by(|a, b| b.len().cmp(&a.len()));
    final_clusters.retain(|cluster| cluster.len() >= args.min_cluster_size);

    log::info!("Final result: {} total clusters across {} k-mer groups", final_clusters.len(), current_clusters.len());

    // Step 4: Debugging - Compare all final clusters and output mismatch statistics
    log::info!("=== Debugging: Pairwise cluster comparison ===");

    // Build consensus for each final cluster
    let final_consensus: Vec<Vec<ConsensusSnpmer>> = final_clusters
        .iter()
        .map(|cluster| build_consensus_snpmers(cluster, twin_reads, k))
        .collect();

    // Compare all pairs of clusters
    for i in 0..final_clusters.len() {
        let rep_i = final_clusters[i][0]; // First member as representative

        for j in (i + 1)..final_clusters.len() {
            let rep_j = final_clusters[j][0];

            // Count matches and mismatches
            let (matches, mismatches) = compare_consensus_snpmers(
                &final_consensus[i],
                &final_consensus[j]
            );

            if matches > 0 || mismatches > 0 {
                log::debug!(
                    "Cluster {} (rep: {}, size: {}) vs Cluster {} (rep: {}, size: {}): {} matches, {} mismatches",
                    i, rep_i, final_clusters[i].len(),
                    j, rep_j, final_clusters[j].len(),
                    matches, mismatches
                );
            }
        }
    }

    log::info!("=== End debugging ===");

    // Step 5: Return flattened clusters
    final_clusters
}
