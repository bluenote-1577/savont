use crate::cli::ClusterArgs as Cli;
use std::sync::{Arc, Mutex};
use rayon::prelude::*;
use crate::types::*;
use fxhash::FxHashMap;
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;

/// Polymorphic marker type for clustering
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PolyMarkerType {
    Snpmer,
    Blockmer,
}

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
        let min_count = candidate_kmer_count.max(query_kmer_count);

        if min_count == 0 {
            continue;
        }

        let ratio = shared_count as f64 / min_count as f64;
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
    let threshold = 0.950;

    // MinHash LSH parameters
    let use_bucketed = true;  // Flag to enable/disable bucketed approach
    let num_hash_tables = 20; // Number of hash functions
    let bucket_size = 3;      // Number of minimizers per bucket
    let top_n_candidates = 10; // Number of top candidates to verify

    // Inverted index: kmer -> Vec<read_id>
    let mut index: FxHashMap<Kmer48, Vec<usize>> = FxHashMap::default();

    // Bucket index for LSH
    let mut bucket_index: BucketIndex = vec![FxHashMap::default(); num_hash_tables];

    // Cluster assignments: read_id -> representative_read_id
    let cluster_assignment: Mutex<HashMap<usize, usize>> = Mutex::new(HashMap::new());

    // Representative reads (reads added to the index)
    let mut representatives: Vec<usize> = Vec::new();

    // Process reads sequentially
    for (read_id, read) in twin_reads.iter().enumerate() {
        let read_kmers = read.minimizer_kmers();

        let best_rep = if use_bucketed {
            // Use bucketed LSH approach
            let bucket_hits = query_read_against_bucket_index(&bucket_index, &read_kmers, bucket_size);

            if bucket_hits.is_empty() {
                None
            } else {
                // Get top candidates sorted by number of bucket hits
                let mut candidates: Vec<(usize, usize)> = bucket_hits.into_iter().collect();
                candidates.sort_by(|a, b| b.1.cmp(&a.1)); // Sort by hits descending
                //println!("Top candidates for read {}: {:?}", read_id, &candidates[..candidates.len().min(5)]);

                // Find the maximum number of hits
                let max_hits = candidates[0].1;

                // Collect all candidates with max hits, plus top_n_candidates
                let mut candidates_to_check = Vec::new();
                for (cand_id, hits) in candidates.iter() {
                    if *hits == max_hits || candidates_to_check.len() < top_n_candidates {
                        candidates_to_check.push(*cand_id);
                    } else {
                        break;
                    }
                }

                // Verify candidates with full similarity check
                let mut best_similarity = 0.0;
                let mut best_candidate = None;

                let read_kmer_set: std::collections::HashSet<_> = read_kmers.iter().cloned().collect();

                for &cand_id in &candidates_to_check {
                    let rep_kmers = twin_reads[cand_id].minimizer_kmers();
                    let mut count = 0;

                    for kmer in read_kmer_set.iter() {
                        if rep_kmers.contains(kmer) {
                            count += 1;
                        }
                    }

                    let ratio = count as f64 / read_kmer_set.len().max(rep_kmers.len()) as f64;
                    let similarity = ratio.powf(1.0 / k as f64);

                    if similarity > best_similarity {
                        best_similarity = similarity;
                        best_candidate = Some(cand_id);
                    }
                }

                if best_similarity > threshold {
                    best_candidate
                } else {
                    None
                }
            }
        } else {
            // Use standard index-based approach for small representative sets
            let similarities = query_read_against_index(&index, &read_kmers, twin_reads, k);

            if let Some(&(representative_id, similarity)) = similarities.first() {
                if similarity > threshold {
                    Some(representative_id)
                } else {
                    None
                }
            } else {
                None
            }
        };

        if let Some(rep_id) = best_rep {
            // Assign this read to the cluster of the representative
            cluster_assignment.lock().unwrap().insert(read_id, rep_id);
        } else {
            // No match found - this read becomes a new representative
            add_read_to_index(&mut index, read_id, &read_kmers);
            if use_bucketed {
                add_read_to_bucket_index(&mut bucket_index, read_id, &read_kmers, bucket_size);
            }
            cluster_assignment.lock().unwrap().insert(read_id, read_id); // Represents itself
            representatives.push(read_id);
        }

        if read_id % 10000 == 0 && read_id > 0 {
            log::info!(
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

    let cluster_assignment = cluster_assignment.into_inner().unwrap();

    // Build clusters from assignments
    let mut clusters_map: HashMap<usize, Vec<usize>> = HashMap::new();
    for (read_id, rep_id) in cluster_assignment {
        clusters_map.entry(rep_id).or_insert_with(Vec::new).push(read_id);
    }

    let mut clusters: Vec<Vec<usize>> = clusters_map.into_values().collect();
    clusters.sort_by(|a, b| b.len().cmp(&a.len()));

    // Sort members within each cluster because lower IDs have better estimated accuracy
    for cluster in clusters.iter_mut(){
        cluster.sort();
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

/// Bucket index for MinHash LSH clustering
/// Each hash table maps bucket signature -> list of representative read IDs
type BucketIndex = Vec<FxHashMap<u64, Vec<usize>>>;

/// Create bucket signature from minimizers using a specific hash seed
/// Takes the bottom `bucket_size` minimizers after hashing and combines them into a single u64
fn create_bucket_signature(minimizers: &[Kmer48], hash_seed: u64, bucket_size: usize) -> Option<u64> {
    if minimizers.len() < bucket_size {
        return None;
    }

    use fxhash::FxHasher64;
    use std::hash::Hasher;

    // Hash each minimizer with the seed and sort by hash value
    let mut hashed: Vec<(u64, Kmer48)> = minimizers
        .iter()
        .map(|&kmer| {
            let mut hasher = FxHasher64::default();
            hasher.write_u64(hash_seed);
            hasher.write_u64(kmer.to_u64());
            (hasher.finish(), kmer)
        })
        .collect();

    hashed.sort_by_key(|x| x.0);

    // Take bottom bucket_size minimizers and combine into signature
    let mut signature: u64 = 0;
    for i in 0..bucket_size {
        // XOR the k-mer values together
        signature ^= hashed[i].1.to_u64().wrapping_mul(i as u64 + 1);
    }

    Some(signature)
}

/// Add a representative read to the bucket index
fn add_read_to_bucket_index(
    bucket_index: &mut BucketIndex,
    read_id: usize,
    minimizers: &[Kmer48],
    bucket_size: usize,
) {
    let _num_tables = bucket_index.len();

    for (table_idx, table) in bucket_index.iter_mut().enumerate() {
        let hash_seed = table_idx as u64;
        if let Some(signature) = create_bucket_signature(minimizers, hash_seed, bucket_size) {
            table.entry(signature).or_insert_with(Vec::new).push(read_id);
        }
    }
}

/// Query a read against the bucket index to find candidate representatives
/// Returns a map of candidate_id -> number of bucket hits
fn query_read_against_bucket_index(
    bucket_index: &BucketIndex,
    minimizers: &[Kmer48],
    bucket_size: usize,
) -> FxHashMap<usize, usize> {
    use rayon::prelude::*;

    let num_tables = bucket_index.len();

    // Query each table in parallel
    let hits_per_table: Vec<FxHashMap<usize, usize>> = (0..num_tables)
        .into_par_iter()
        .map(|table_idx| {
            let hash_seed = table_idx as u64;
            let mut local_hits: FxHashMap<usize, usize> = FxHashMap::default();

            if let Some(signature) = create_bucket_signature(minimizers, hash_seed, bucket_size) {
                if let Some(candidates) = bucket_index[table_idx].get(&signature) {
                    for &candidate_id in candidates {
                        *local_hits.entry(candidate_id).or_insert(0) += 1;
                    }
                }
            }

            local_hits
        })
        .collect();

    // Merge hits from all tables
    let mut total_hits: FxHashMap<usize, usize> = FxHashMap::default();
    for table_hits in hits_per_table {
        for (candidate_id, count) in table_hits {
            *total_hits.entry(candidate_id).or_insert(0) += count;
        }
    }

    total_hits
}

/// Add a read's SNPmers to the SNPmer inverted index
fn add_read_snpmers_to_index(
    index: &mut FxHashMap<u64, Vec<(usize, Kmer48)>>,
    read_id: usize,
    read_snpmers: &[Kmer48],
    k: usize,
) {
    let mask = !(3 << (k - 1));
    for &kmer in read_snpmers {
        let splitmer = kmer.to_u64() & mask;
        index.entry(splitmer).or_insert_with(Vec::new).push((read_id, kmer));
    }
}

/// Query a read's SNPmers against the SNPmer index
/// Returns (matches, mismatches) for each candidate read_id
/// Candidates with 0 mismatches are compatible
pub fn find_compatible_candidates(
    index: &FxHashMap<u64, Vec<(usize, Kmer48)>>,
    query_snpmers: &[Kmer48],
    k: usize,
) -> FxHashMap<usize, (usize, usize)> {
    let mask = !(3 << (k - 1));

    // Track matches and mismatches for each candidate
    let mut candidate_stats: FxHashMap<usize, (usize, usize)> = FxHashMap::default();

    // Check each query SNPmer
    for &query_kmer in query_snpmers {
        let query_splitmer = query_kmer.to_u64() & mask;

        if let Some(candidates) = index.get(&query_splitmer) {
            for &(candidate_id, candidate_kmer) in candidates {
                let stats = candidate_stats.entry(candidate_id).or_insert((0, 0));
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

/// Assign a read to a representative cluster
fn assign_read_to_representative(
    read_id: usize,
    rep_id: usize,
    local_assignment: &mut HashMap<usize, usize>,
    rep_size: &mut HashMap<usize, usize>,
) {
    local_assignment.insert(read_id, rep_id);
    *rep_size.entry(rep_id).or_insert(0) += 1;
}

/// Create a new representative for a read
fn create_new_representative(
    read_id: usize,
    read_snpmers: &[Kmer48],
    k: usize,
    representatives: &mut Vec<usize>,
    snpmer_index: &mut FxHashMap<u64, Vec<(usize, Kmer48)>>,
    local_assignment: &mut HashMap<usize, usize>,
    rep_size: &mut HashMap<usize, usize>,
) {
    representatives.push(read_id);
    add_read_snpmers_to_index(snpmer_index, read_id, read_snpmers, k);
    local_assignment.insert(read_id, read_id);
    rep_size.insert(read_id, 1);
}

/// Find best representative using iterative parallel search
fn find_best_representative_iterative(
    read_id: usize,
    splitmer_to_kmer: &FxHashMap<u64, Kmer48>,
    representatives: &[usize],
    twin_reads: &[TwinRead],
    k: usize,
    args: &Cli,
) -> Option<usize> {
    let mask = !(3 << (k - 1));

    representatives.par_iter().with_max_len(1).find_any(|&&rep_id| {
        let rep_snpmers = twin_reads[rep_id].snpmer_kmers();

        // Check SNPmer compatibility
        let mut matches = 0;
        let mut mismatches = 0;

        for &rep_kmer in rep_snpmers.iter() {
            let rep_splitmer = rep_kmer.to_u64() & mask;
            if let Some(&query_kmer) = splitmer_to_kmer.get(&rep_splitmer) {
                if query_kmer == rep_kmer {
                    matches += 1;
                } else {
                    mismatches += 1;
                }
            }
        }

        // Compatible if no mismatches and at least one match
        if mismatches == 0 && matches > 0{
            // Additional blockmer validation if enabled and representative found
            if args.use_blockmers {
                let blockmer_comp = compare_blockmers(
                    &twin_reads[read_id],
                    &twin_reads[rep_id],
                    k,
                    args.blockmer_length,
                );
                if blockmer_comp.1 > args.blockmer_length {
                    false
                } else {
                    true
                }
            } else {
                true
            }
        }
        else{
            false
        }
    }).copied()
}

/// Find best representative using index-based search with blockmer validation
fn find_best_representative_indexed(
    read_id: usize,
    read_snpmers: &[Kmer48],
    snpmer_index: &FxHashMap<u64, Vec<(usize, Kmer48)>>,
    rep_size: &HashMap<usize, usize>,
    twin_reads: &[TwinRead],
    k: usize,
    blockmer_length: usize,
    use_blockmers: bool,
) -> Option<usize> {
    // Query this read's SNPmers against the index
    let candidate_stats = find_compatible_candidates(snpmer_index, read_snpmers, k);

    // Filter to only compatible candidates (0 mismatches, >0 matches)
    let compatible_candidates: Vec<(usize, usize, usize)> = candidate_stats
        .iter()
        .filter(|(_, (matches, mismatches))| *mismatches == 0 && *matches > 0)
        .map(|(candidate_id, (matches, _))| (*candidate_id, *matches, rep_size[candidate_id]))
        .collect();

    if compatible_candidates.is_empty() {
        return None;
    }

    // Sort by: most matches, then fewest members (smallest cluster)
    let mut candidates_sorted: Vec<_> = compatible_candidates
        .iter()
        .map(|(cand_id, matches, size)| (-(*matches as i64), *size, *cand_id))
        .collect();
    candidates_sorted.sort();

    // Validate with blockmers if enabled
    if use_blockmers {
        validate_candidates_with_blockmers(
            read_id,
            &candidates_sorted,
            twin_reads,
            k,
            blockmer_length,
        )
    } else {
        Some(candidates_sorted[0].2)
    }
}

/// Validate candidates using blockmer comparison
fn validate_candidates_with_blockmers(
    read_id: usize,
    candidates_sorted: &[(i64, usize, usize)],
    twin_reads: &[TwinRead],
    k: usize,
    blockmer_length: usize,
) -> Option<usize> {
    let mut blockmer_cands = Vec::new();

    for candidate in candidates_sorted {
        let blockmer_comp = compare_blockmers(
            &twin_reads[read_id],
            &twin_reads[candidate.2],
            k,
            blockmer_length,
        );

        log::trace!(
            "Read {} compatible with candidate {} ({} matches, cluster size {}). Blockmer match/mismatch = {} / {}",
            twin_reads[read_id].base_id,
            twin_reads[candidate.2].base_id,
            -candidate.0,
            candidate.1,
            blockmer_comp.0,
            blockmer_comp.1
        );

        blockmer_cands.push((candidate.2, blockmer_comp.0, blockmer_comp.1));
    }

    // Sort by fewest mismatches, then most matches
    blockmer_cands.sort_by(|a, b| a.2.cmp(&b.2).then(b.1.cmp(&a.1)));

    // Check if best candidate passes blockmer validation
    if blockmer_cands[0].2 > 1 {
        log::trace!(
            "Read {} has no fully concordant candidates; creating new SNPmer representative",
            twin_reads[read_id].base_id,
        );
        None
    } else {
        Some(blockmer_cands[0].0)
    }
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

    // Determine whether to use iterative mode (parallel search over representatives)

    // Process each k-mer cluster independently in parallel
    kmer_clusters.par_iter().with_max_len(1).enumerate().for_each(|(kmer_cluster_id, kmer_cluster)| {
         // Skip empty k-mer clusters
        if kmer_cluster.len() < 1 {
            return;
        }

        // SNPmer inverted index for this k-mer cluster: splitmer -> Vec<(read_id, full_kmer)>
        let mut snpmer_index: FxHashMap<u64, Vec<(usize, Kmer48)>> = FxHashMap::default();

        // Local cluster assignments within this k-mer cluster
        let mut local_assignment: HashMap<usize, usize> = HashMap::new();

        // Track representative reads for iterative mode
        let mut representatives: Vec<usize> = Vec::new();

        let mut rep_size: HashMap<usize, usize> = HashMap::new();
        let mut count = 0;
        
        for &read_id in kmer_cluster {
            let read_snpmers = twin_reads[read_id].snpmer_kmers();

            // Determine which search strategy to use
            let use_iterative = representatives.len() > 1000;
            //let use_iterative = false;

            // Find best matching representative
            let best_rep = if use_iterative {
                let splitmer_to_kmer: FxHashMap<u64, Kmer48> = read_snpmers
                    .iter()
                    .map(|&kmer| (kmer.to_u64() & !(3 << (k - 1)), kmer))
                    .collect();

                find_best_representative_iterative(
                    read_id,
                    &splitmer_to_kmer,
                    &representatives,
                    twin_reads,
                    k,
                    args,
                )
                
            } else {
                find_best_representative_indexed(
                    read_id,
                    read_snpmers,
                    &snpmer_index,
                    &rep_size,
                    twin_reads,
                    k,
                    args.blockmer_length,
                    args.use_blockmers,
                )
            };

            // Assign to representative or create new one
            if let Some(rep_id) = best_rep {
                assign_read_to_representative(read_id, rep_id, &mut local_assignment, &mut rep_size);
            } else {
                create_new_representative(
                    read_id,
                    read_snpmers,
                    k,
                    &mut representatives,
                    &mut snpmer_index,
                    &mut local_assignment,
                    &mut rep_size,
                );
            }

            count += 1;
            if count % 10_000 == 0 {
                log::info!("Processed {} / {} reads for SNPmer clustering in k-mer cluster {} with {} reps",  count, kmer_cluster.len(), kmer_cluster_id, representatives.len());
            }

            // if count % 1_000_000 == 0 {
            //     //log::info!("Processed {} reads for SNPmer clustering", count);

            //     // Perform intermediate reclustering to merge fragmented clusters
            //     log::info!("Count {} / {} for cluster {}: performing intermediate reclustering within k-mer cluster {}", count, kmer_cluster.len(), kmer_cluster_id, kmer_cluster_id);

            //     // Step 1: Build current clusters from local_assignment
            //     let mut cluster_map = HashMap::new();
            //     for (read_id, rep_id) in &local_assignment {
            //         if cluster_map.contains_key(rep_id) == false {
            //             cluster_map.insert(*rep_id, vec![rep_id.clone()]);
            //         }
            //         else{
            //             cluster_map.entry(*rep_id).or_insert_with(Vec::new).push(*read_id);
            //         }
            //     }

            //     let mut current_clusters: Vec<Vec<usize>> = cluster_map.into_values().collect();
            //     current_clusters.sort_by(|a, b| b.len().cmp(&a.len()));

            //     log::trace!("Before intermediate recluster: {} clusters", current_clusters.len());

            //     // Step 2: Perform one round of reclustering using top 100 reads for consensus
            //     let (merged_clusters, num_merges) = recluster_one_round_top_n(
            //         current_clusters,
            //         twin_reads,
            //         k,
            //         args.blockmer_length,
            //         PolyMarkerType::Snpmer,
            //         Some(200), // Use only top 100 reads for fast consensus building
            //     );

            //     log::debug!("Intermediate recluster for kmer group {} : merged {} clusters, now have {} clusters",
            //         kmer_cluster_id, num_merges, merged_clusters.len());

            //     // Step 3: Rebuild state from merged clusters
            //     // Clear the index and rebuild from scratch
            //     snpmer_index.clear();
            //     local_assignment.clear();
            //     rep_size.clear();

            //     // Rebuild from merged clusters
            //     for cluster in &merged_clusters {
            //         if cluster.is_empty() {
            //             continue;
            //         }

            //         let rep_id = cluster[0]; // First read becomes the representative

            //         // Add representative's SNPmers to index
            //         let rep_snpmers = twin_reads[rep_id].snpmer_kmers();
            //         add_read_snpmers_to_index(&mut snpmer_index, rep_id, &rep_snpmers, k);

            //         // Assign all reads in cluster to this representative
            //         for &read_id in cluster {
            //             local_assignment.insert(read_id, rep_id);
            //         }

            //         rep_size.insert(rep_id, cluster.len());
            //     }

            //     log::trace!("Rebuilt index with {} representatives", snpmer_index.len());
            // }
        }

        // Add representatives
        let mut cluster_map = HashMap::new();
        for (read_id, rep_id) in &local_assignment {
            if read_id != rep_id {
                continue;
            }
            cluster_map.entry(*rep_id).or_insert_with(Vec::new).push(*rep_id);
        }
        // Add non representatives
        for (read_id, rep_id) in local_assignment.iter() {
            if read_id == rep_id {
                continue;
            }
            cluster_map.entry(*rep_id).or_insert_with(Vec::new).push(*read_id);
        }

        let mut local_clusters: Vec<Vec<usize>> = cluster_map.into_values().collect();
        local_clusters.sort_by(|a, b| b.len().cmp(&a.len()));

        // Remove small clusters
        log::debug!("Before size filtering: {} SNPmer clusters in k-mer cluster {}", local_clusters.len(), kmer_cluster_id);
        local_clusters.retain(|cluster| cluster.len() >= args.min_cluster_size);
        log::debug!("After size filtering: {} SNPmer clusters in k-mer cluster {}", local_clusters.len(), kmer_cluster_id);

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
            log::info!(
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
            cluster.iter().map(|x| format!("{} {}", &twin_reads[*x].id, &twin_reads[*x].est_id.unwrap_or(100.))).collect::<Vec<_>>().join("\n")
        ).unwrap();
    }

    return local_clusters_all;
}

fn compare_blockmers(
    twin_read1: &TwinRead,
    twin_read2: &TwinRead,
    _k: usize,
    l: usize
) -> (usize, usize) {
    let blockmers1 = twin_read1.blockmers_vec();
    let blockmers2 = twin_read2.blockmers_vec();

    let mut matches = 0;
    let mut mismatches = 0;

    let mut map2: FxHashMap<u64, u64> = FxHashMap::default();
    for &(_pos, kmer) in &blockmers2 {
        let anchor = kmer >> l * 2;
        map2.insert(anchor, kmer);
    }

    for &(_pos, kmer1) in &blockmers1 {
        let anchor = kmer1 >> l * 2;
        if let Some(&kmer2) = map2.get(&anchor) {
            if kmer1 == kmer2 {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
    }

    (matches, mismatches)
}


/// Build consensus SNPmer representative for a cluster
fn build_consensus_snpmers(
    cluster: &[usize],
    twin_reads: &[TwinRead],
    k: usize,
) -> Vec<ConsensusSnpmer> {
    build_consensus_snpmers_top_n(cluster, twin_reads, k, None)
}

/// Build consensus SNPmer representative for a cluster using only top N reads
fn build_consensus_snpmers_top_n(
    cluster: &[usize],
    twin_reads: &[TwinRead],
    k: usize,
    top_n: Option<usize>,
) -> Vec<ConsensusSnpmer> {
    let mask = !(3 << (k - 1));

    // Map: splitmer -> Map: full_kmer -> (count, Vec<positions>)
    let mut splitmer_data: FxHashMap<u64, FxHashMap<Kmer48, (usize, Vec<u32>)>> = FxHashMap::default();

    // Determine how many reads to use for consensus building
    let reads_to_use = if let Some(n) = top_n {
        cluster.len().min(n)
    } else {
        cluster.len()
    };

    // Step 1: Collect SNPmers from top N reads in the cluster
    for &read_id in &cluster[..reads_to_use] {
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
            if *count >= (cluster.len()/6).max(1) {
                // Calculate median position
                let mut pos_sorted = positions.clone();
                pos_sorted.sort();
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

fn build_consensus_blockmers(
    cluster: &[usize],
    twin_reads: &[TwinRead],
    k: usize,
    l: usize,
) -> Vec<ConsensusPoly> {
    build_consensus_blockmers_top_n(cluster, twin_reads, k, l, None)
}

fn build_consensus_blockmers_top_n(
    cluster: &[usize],
    twin_reads: &[TwinRead],
    _k: usize,
    l: usize,
    top_n: Option<usize>,
) -> Vec<ConsensusPoly> {
    // For blockmers: splitmer is the anchor k-mer (shift right by 2*l to remove suffix)
    // full k-mer is the entire (k+l)-mer

    // Map: splitmer (anchor) -> Map: full_kmer -> (count, Vec<positions>)
    let mut splitmer_data: FxHashMap<u64, FxHashMap<Kmer48, (usize, Vec<u32>)>> = FxHashMap::default();

    // Determine how many reads to use for consensus building
    let reads_to_use = if let Some(n) = top_n {
        cluster.len().min(n)
    } else {
        cluster.len()
    };

    // Step 1: Collect blockmers from top N reads in the cluster
    for &read_id in &cluster[..reads_to_use] {
        let blockmers = twin_reads[read_id].blockmers_vec();
        for &(pos, kmer_u64) in &blockmers {
            // Extract anchor k-mer by shifting right by 2*l
            let splitmer = kmer_u64 >> (2 * l);
            let kmer = Kmer48::from_u64(kmer_u64);

            let entry = splitmer_data
                .entry(splitmer)
                .or_insert_with(FxHashMap::default)
                .entry(kmer)
                .or_insert((0, Vec::new()));
            entry.0 += 1;
            entry.1.push(pos);
        }
    }

    // Step 2: For each splitmer (anchor), find the most common full blockmer and compute median position
    let mut consensus_blockmers = Vec::new();
    for (splitmer, kmer_data) in splitmer_data {
        // Find the blockmer with maximum count
        if let Some((&best_kmer, (count, positions))) = kmer_data.iter().max_by_key(|(_, (count, _))| count) {
            if *count >= (cluster.len()/6).max(1) {
                // Calculate median position
                let mut pos_sorted = positions.clone();
                pos_sorted.sort();
                let median_pos = if !pos_sorted.is_empty() {
                    pos_sorted[pos_sorted.len() / 2]
                } else {
                    0
                };
                consensus_blockmers.push(ConsensusPoly::new(median_pos, splitmer, best_kmer, *count as u32));
            }
        }
    }

    consensus_blockmers.sort_by_key(|cs| (cs.position, cs.splitmer));
    consensus_blockmers
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
    mismatches == 0 && matches >= consensus1.len().min(consensus2.len().max(2))
}

/// Reassign reads to their best matching cluster based on SNPmer comparison
/// Returns (reassigned_clusters, num_reassignments)
fn reassign_reads_to_best_cluster(
    current_clusters: Vec<Vec<usize>>,
    twin_reads: &[TwinRead],
    k: usize,
    l: usize,
    marker_type: PolyMarkerType,
    args: &Cli,
) -> (Vec<Vec<usize>>, usize) {
    // Build consensus for each cluster based on marker type
    let cluster_consensus: Vec<Vec<ConsensusPoly>> = current_clusters
        .iter()
        .map(|cluster| {
            match marker_type {
                PolyMarkerType::Snpmer => build_consensus_snpmers(cluster, twin_reads, k),
                PolyMarkerType::Blockmer => build_consensus_blockmers(cluster, twin_reads, k, l),
            }
        })
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

    // Create new cluster assignments
    let new_clusters: Mutex<Vec<Vec<usize>>> = Mutex::new(vec![Vec::new(); current_clusters.len()]);
    let num_reassignments = Mutex::new(0);

    // For each cluster, check each read
    current_clusters.par_iter().enumerate().for_each(|(cluster_idx, cluster)| {
        let mut new_clusters_local = vec![Vec::new(); current_clusters.len()];
        for &read_id in cluster {
            // Get read markers based on marker type
            let read_markers: Vec<(u32, Kmer48)> = match marker_type {
                PolyMarkerType::Snpmer => twin_reads[read_id].snpmers_vec(),
                PolyMarkerType::Blockmer => {
                    twin_reads[read_id].blockmers_vec()
                        .into_iter()
                        .map(|(pos, kmer_u64)| (pos, Kmer48::from_u64(kmer_u64)))
                        .collect()
                }
            };

            // Find the best matching cluster for this read
            let mut best_cluster = cluster_idx;
            let mut best_score = (usize::MAX, 0); // (mismatches, matches) - lower mismatches and higher matches is better

            for (candidate_idx, _) in cluster_consensus.iter().enumerate() {
                // Build index of consensus splitmer -> kmer
                let splitmer_to_kmer = &splitmer_to_consensus_kmers[candidate_idx];

                let mut matches = 0;
                let mut mismatches = 0;

                // Compare read's markers against this consensus
                for &(_pos, kmer) in &read_markers {
                    // Extract splitmer based on marker type
                    let splitmer = match marker_type {
                        PolyMarkerType::Snpmer => {
                            let mask = !(3 << (k - 1));
                            kmer.to_u64() & mask
                        }
                        PolyMarkerType::Blockmer => {
                            // Anchor k-mer is top k bases, shift right by 2*l
                            kmer.to_u64() >> (2 * l)
                        }
                    };

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
            //new_clusters.lock().unwrap()[best_cluster].push(read_id);
            new_clusters_local[best_cluster].push(read_id);

            if best_cluster != cluster_idx {
                *num_reassignments.lock().unwrap() += 1;
                log::trace!(
                    "Reassigned read {} from cluster {} to cluster {} (mismatches: {}, matches: {})",
                    read_id, cluster_idx, best_cluster, best_score.0, best_score.1
                );
            }
        }

        // Merge local new clusters into global new clusters
        let mut global_clusters = new_clusters.lock().unwrap();
        for (i, cluster) in new_clusters_local.into_iter().enumerate() {
            global_clusters[i].extend(cluster);
        }
    });

    // Remove empty clusters
    let mut new_clusters = new_clusters.into_inner().unwrap();
    new_clusters.retain(|cluster| !cluster.is_empty() && cluster.len() >= args.min_cluster_size);
    let num_reassignments = num_reassignments.into_inner().unwrap();

    log::trace!("Reassignment complete: {} reads reassigned to better clusters", num_reassignments);

    (new_clusters, num_reassignments)
}

/// Perform one round of reclustering based on consensus representatives
/// Returns (merged_clusters, num_merges_performed)
fn recluster_one_round(
    current_clusters: Vec<Vec<usize>>,
    twin_reads: &[TwinRead],
    k: usize,
    l: usize,
    marker_type: PolyMarkerType,
) -> (Vec<Vec<usize>>, usize) {
    recluster_one_round_top_n(current_clusters, twin_reads, k, l, marker_type, None)
}

/// Perform one round of reclustering using only top N reads for consensus building
/// Returns (merged_clusters, num_merges_performed)
fn recluster_one_round_top_n(
    current_clusters: Vec<Vec<usize>>,
    twin_reads: &[TwinRead],
    k: usize,
    l: usize,
    marker_type: PolyMarkerType,
    top_n: Option<usize>,
) -> (Vec<Vec<usize>>, usize) {
    // Build consensus representatives for each cluster ONCE before the loop
    let mut all_clusters: Vec<(Vec<usize>, Vec<ConsensusPoly>)> = Vec::new();

    for cluster in current_clusters {
        if cluster.is_empty() {
            continue;
        }

        let consensus = match marker_type {
            PolyMarkerType::Snpmer => build_consensus_snpmers_top_n(&cluster, twin_reads, k, top_n),
            PolyMarkerType::Blockmer => build_consensus_blockmers_top_n(&cluster, twin_reads, k, l, top_n),
        };
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
            all_clusters[i].1 = match marker_type {
                PolyMarkerType::Snpmer => build_consensus_snpmers_top_n(&all_clusters[i].0, twin_reads, k, top_n),
                PolyMarkerType::Blockmer => build_consensus_blockmers_top_n(&all_clusters[i].0, twin_reads, k, l, top_n),
            };
            cluster_needs_consensus_rebuild[i] = false;
        }

        // Try to merge smaller clusters into this one
        for j in i + 1..all_clusters.len() {
            if cluster_merged[j] {
                continue;
            }

            // Check if consensus representatives are concordant (bidirectional)
            // Use the pre-built consensus from the beginning of the function
            let consensus_i = &all_clusters[i].1;
            let consensus_j = &all_clusters[j].1;
            let mut concordant = {
                are_consensus_concordant(consensus_i, consensus_j) &&
                are_consensus_concordant(consensus_j, consensus_i)
            };

            let (matches, mismatches) = compare_consensus_snpmers(consensus_i, consensus_j);
            let max_len = all_clusters[i].0.len().max(all_clusters[j].0.len());
            let min_len = all_clusters[i].0.len().min(all_clusters[j].0.len());
            let _missing = consensus_i.len().min(consensus_j.len()) - matches;
            if mismatches == 0  && 
            (matches as f64) > (consensus_i.len().min(consensus_j.len()) as f64) * 0.975 && 
            max_len / min_len > 50 {
                concordant = true;
                //println!("Note: Clusters {} and {} potential merging due to size disparity vs mismatches: max_len {}, min_len {}, matches {}, mismatches {}, cons len 1 {}, cons len 2 {}", 
                //i, j, max_len, min_len, matches, mismatches, consensus_i.len(), consensus_j.len());
            }
            
            if mismatches == 0 && max_len / min_len > 500 && min_len <= 2 {
                concordant = true;
                //println!("Note: Clusters {} and {} potential merging due to extreme size disparity vs mismatches: max_len {}, min_len {}, matches {}, mismatches {}, cons len 1 {}, cons len 2 {}", 
                //i, j, max_len, min_len, matches, mismatches, consensus_i.len(), consensus_j.len());
            }

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
                log::trace!(
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
            all_clusters[i].1 = match marker_type {
                PolyMarkerType::Snpmer => build_consensus_snpmers_top_n(&all_clusters[i].0, twin_reads, k, top_n),
                PolyMarkerType::Blockmer => build_consensus_blockmers_top_n(&all_clusters[i].0, twin_reads, k, l, top_n),
            };
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
    let marker_type = if args.use_blockmers {
        PolyMarkerType::Blockmer
    } else {
        PolyMarkerType::Snpmer
    };

    // Preserve the hierarchical structure: FxHashMap<kmer_cluster_id, Vec<snpmer_clusters>>
    let mut current_clusters = clusters;

    let total_groups = current_clusters.len();
    let total_initial_clusters: usize = current_clusters.values().map(|v| v.len()).sum();
    log::info!("Starting with {} k-mer groups containing {} total SNPmer clusters", total_groups, total_initial_clusters);

    // Step 2: Iteratively recluster until convergence
    let mut iteration = 0;
    loop {
        if iteration >= args.max_iterations_recluster {
            log::info!("Reached maximum reclustering iterations ({})", args.max_iterations_recluster);
            break;
        }
        iteration += 1;
        let total_merges = Mutex::new(0);
        let total_reassignments = Mutex::new(0);

        let new_clusters: Mutex<FxHashMap<usize, Vec<Vec<usize>>>> = Mutex::new(FxHashMap::default());


        // Process each k-mer group independently
        //for (kmer_cluster_id, snpmer_clusters) in current_clusters {
        current_clusters.into_par_iter().for_each(|(kmer_cluster_id, snpmer_clusters)| {
            log::trace!("Processing k-mer group {} with {} SNPmer clusters", kmer_cluster_id, snpmer_clusters.len());

            // Step 2a: Merge clusters within this group
            let (merged_clusters, num_merges) = recluster_one_round(
                snpmer_clusters,
                twin_reads,
                k,
                args.blockmer_length,
                marker_type, // TODO: Make this configurable
            );

            *total_merges.lock().unwrap() += num_merges;
            log::trace!("Processing k-mer group {} --- merged", kmer_cluster_id);

            // Step 2b: Reassign reads to best matching clusters within this group
            let reassign = true;
            if reassign{
                let (reassigned_clusters, num_reassignments) = reassign_reads_to_best_cluster(
                    merged_clusters,
                    twin_reads,
                    k,
                    args.blockmer_length,
                    marker_type, // TODO: Make this configurable
                    args,
                );

                *total_reassignments.lock().unwrap() += num_reassignments;

                // Store the updated clusters for this k-mer group
                if !reassigned_clusters.is_empty() {
                    new_clusters.lock().unwrap().insert(kmer_cluster_id, reassigned_clusters);
                }
            }
            else{
                // Store the merged clusters without reassignment
                if !merged_clusters.is_empty() {
                    new_clusters.lock().unwrap().insert(kmer_cluster_id, merged_clusters);
                }
            }
            log::trace!("Processing k-mer group {} --- done", kmer_cluster_id);
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
    if log::log_enabled!(log::Level::Trace) {
        log::info!("Performing detailed pairwise cluster comparison for debugging...");
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
    }

    // Step 5: Return flattened clusters
    final_clusters
}
