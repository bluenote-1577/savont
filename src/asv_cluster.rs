use crate::cli::ClusterArgs as Cli;
use crate::constants::LSH_NUM_TABLES;
use crate::types::*;
use fxhash::{FxHashMap, FxHashSet};
use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

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
        // Unfiltered count, matching the minimizers this index is built from. (Only
        // reachable when use_bucketed is false, which it never is.)
        let candidate_kmer_count = twin_reads[candidate_id].minimizer_count();
        let query_kmer_count = query_kmers.len();
        let min_count = candidate_kmer_count.max(query_kmer_count);

        if min_count == 0 {
            continue;
        }

        let ratio = shared_count as f64 / min_count as f64;
        let similarity = ratio.powf(1.0 / k as f64);

        similarities.push((candidate_id, similarity));
    }

    similarities.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap().then_with(|| a.0.cmp(&b.0)));
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
    let threshold = args.primary_clustering_threshold;

    // MinHash LSH parameters
    let use_bucketed = true; // Flag to enable/disable bucketed approach
    let top_n_candidates = 10; // Number of top candidates to verify

    // Inverted index: kmer -> Vec<read_id>
    let mut index: FxHashMap<Kmer48, Vec<usize>> = FxHashMap::default();

    // Bucket index for LSH
    let mut bucket_index: BucketIndex = vec![FxHashMap::default(); LSH_NUM_TABLES];

    // Cluster assignments: read_id -> representative_read_id.
    // This loop is sequential, so no lock is needed.
    let mut cluster_assignment: HashMap<usize, usize> = HashMap::new();

    // Representative reads (reads added to the index)
    let mut representatives: Vec<usize> = Vec::new();

    // Minimizer set (and minimizer count) per representative, kept from when the read
    // became a representative. The verification step below used to call `Vec::contains` (a
    // linear scan) once per read k-mer per candidate; probing a prebuilt set instead turns
    // that inner loop from O(|read| * |rep|) comparisons into O(|read|) hash lookups.
    let mut rep_kmer_sets: FxHashMap<usize, (FxHashSet<Kmer48>, usize)> = FxHashMap::default();

    // The loop below is sequential, so everything per-read that does not depend on earlier
    // reads (recomputing the minimizers from the sequence, their LSH signatures and k-mer
    // set) is done for a bounded batch in parallel first; reads are then assigned in their
    // original order. (Building the set inside the loop instead measured ~20% slower.)
    const SIGNATURE_BATCH_SIZE: usize = 4096;
    let read_batches = twin_reads.chunks(SIGNATURE_BATCH_SIZE).flat_map(|batch| {
        batch
            .par_iter()
            .map(|read| {
                let kmers = read.minimizer_kmers();
                let signatures = lsh_signatures(&kmers);
                let kmer_set: FxHashSet<Kmer48> = kmers.iter().cloned().collect();
                (kmers, kmer_set, signatures)
            })
            .collect::<Vec<_>>()
    });
    for (read_id, (read_kmers, read_kmer_set, signatures)) in read_batches.enumerate() {
        let signatures = &signatures[..];

        let best_rep = if use_bucketed {
            // Use bucketed LSH approach
            let bucket_hits = query_read_against_bucket_index(&bucket_index, signatures);

            if bucket_hits.is_empty() {
                None
            } else {
                // Get top candidates sorted by number of bucket hits
                let mut candidates: Vec<(usize, usize)> = bucket_hits.into_iter().collect();
                // Sort by hits descending. This runs once per read inside a sequential
                // loop, so rayon's task overhead outweighs any gain on these short lists.
                // `(hits, id)` is a total order (ids are unique), so unstable sorting is
                // deterministic and matches the previous parallel sort exactly.
                candidates.sort_unstable_by(|a, b| (b.1, b.0).cmp(&(a.1, a.0)));
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

                for &cand_id in &candidates_to_check {
                    // Counts distinct read k-mers present in the representative, exactly as
                    // the previous `rep_kmers.contains(..)` scan did. The denominator keeps
                    // using the minimizer count (which may include duplicates), not the set size.
                    let (rep_set, rep_kmer_count) = &rep_kmer_sets[&cand_id];
                    let count = read_kmer_set.iter().filter(|k| rep_set.contains(k)).count();

                    let ratio = count as f64 / read_kmer_set.len().max(*rep_kmer_count) as f64;
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
            cluster_assignment.insert(read_id, rep_id);
        } else {
            // No match found - this read becomes a new representative
            if use_bucketed {
                add_read_to_bucket_index(&mut bucket_index, read_id, signatures);
            } else {
                add_read_to_index(&mut index, read_id, &read_kmers);
            }
            cluster_assignment.insert(read_id, read_id); // Represents itself
            representatives.push(read_id);
            rep_kmer_sets.insert(read_id, (read_kmer_set, read_kmers.len()));
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

    // (sequential loop above, so this is already a plain map)

    // Build clusters from assignments
    let mut clusters_map: HashMap<usize, Vec<usize>> = HashMap::new();
    for (read_id, rep_id) in cluster_assignment {
        clusters_map
            .entry(rep_id)
            .or_insert_with(Vec::new)
            .push(read_id);
    }

    let mut clusters: Vec<Vec<usize>> = clusters_map.into_values().collect();
    // HashMap iteration randomizes member order. Canonicalize members before
    // using the first member as the equal-size cluster tie-break.
    for cluster in &mut clusters {
        cluster.sort_unstable();
    }
    clusters.sort_by(|a, b| {
        b.len()
            .cmp(&a.len())
            .then_with(|| a.first().cmp(&b.first()))
    });

    // Remove small clusters
    clusters.retain(|cluster| cluster.len() >= args.min_cluster_size);

    // Write clusters to file
    let cluster_file = output_dir.join("kmer_clusters_stage2.tsv");
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
            cluster
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        )
        .unwrap();
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
fn _create_bucket_signature(
    minimizers: &[Kmer48],
    hash_seed: u64,
    bucket_size: usize,
) -> Option<u64> {
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

/// Add a representative read to the bucket index using its precomputed LSH signatures
fn add_read_to_bucket_index(
    bucket_index: &mut BucketIndex,
    read_id: usize,
    signatures: &[Option<u64>],
) {
    for (table_idx, table) in bucket_index.iter_mut().enumerate() {
        if let Some(Some(signature)) = signatures.get(table_idx) {
            table
                .entry(*signature)
                .or_insert_with(Vec::new)
                .push(read_id);
        }
    }
}

/// Query a read against the bucket index using its precomputed LSH signatures.
/// Returns a map of candidate_id -> number of bucket hits.
fn query_read_against_bucket_index(
    bucket_index: &BucketIndex,
    signatures: &[Option<u64>],
) -> FxHashMap<usize, usize> {
    let num_tables = bucket_index.len();

    // Query each table in parallel using precomputed signatures
    let hits_per_table: Vec<FxHashMap<usize, usize>> = (0..num_tables)
        .into_iter()
        .map(|table_idx| {
            let mut local_hits: FxHashMap<usize, usize> = FxHashMap::default();

            if let Some(Some(signature)) = signatures.get(table_idx) {
                if let Some(candidates) = bucket_index[table_idx].get(signature) {
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
        index
            .entry(splitmer)
            .or_insert_with(Vec::new)
            .push((read_id, kmer));
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

/// Find best representative using iterative parallel search.
/// `rep_snpmers[i]` holds the SNPmers of `representatives[i]`.
fn find_best_representative_iterative(
    read_id: usize,
    splitmer_to_kmer: &FxHashMap<u64, Kmer48>,
    representatives: &[usize],
    rep_snpmers: &[Vec<Kmer48>],
    twin_reads: &[TwinRead],
    k: usize,
    args: &Cli,
) -> Option<usize> {
    let mask = !(3 << (k - 1));

    // find_first preserves representative order while retaining parallel search.
    representatives
        .par_iter()
        .zip(rep_snpmers.par_iter())
        .with_max_len(1)
        .find_first(|&(&rep_id, rep_snpmers)| {
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
            if mismatches == 0 && matches > 0 {
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
            } else {
                false
            }
        })
        .map(|(&rep_id, _)| rep_id)
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
    blockmer_cands.sort_by(|a, b| a.2.cmp(&b.2).then(b.1.cmp(&a.1)).then(a.0.cmp(&b.0)));

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
    // In low-polymorphism mode the SNPmer signal is driven by sequencing errors rather than
    // real biological variation. Skip SNPmer-based sub-clustering and pass k-mer clusters
    // through directly so SNPmer-less reads are not discarded.
    if args.low_polymorphism {
        log::info!(
            "Low-polymorphism mode: skipping SNPmer clustering, passing k-mer clusters through."
        );
        let mut clusters: Vec<Vec<usize>> = kmer_clusters
            .iter()
            .filter(|c| c.len() >= args.min_cluster_size)
            .cloned()
            .collect();
        clusters.sort_by(|a, b| {
            b.len()
                .cmp(&a.len())
                .then_with(|| a.first().cmp(&b.first()))
        });
        return clusters;
    }

    log::info!("Starting greedy SNPmer-based clustering within k-mer clusters...");

    let k = args.kmer_size;

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
        // SNPmers of each representative, aligned with `representatives`. SNPmers are
        // recomputed from the sequence on every `snpmer_kmers()` call, and the iterative search
        // compares each read against every representative, so keep one copy per representative.
        let mut rep_snpmers: Vec<Vec<Kmer48>> = Vec::new();

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
                    &rep_snpmers,
                    twin_reads,
                    k,
                    args,
                )
                
            } else {
                find_best_representative_indexed(
                    read_id,
                    &read_snpmers,
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
                    &read_snpmers,
                    k,
                    &mut representatives,
                    &mut snpmer_index,
                    &mut local_assignment,
                    &mut rep_size,
                );
                rep_snpmers.push(read_snpmers);
            }

            count += 1;
            if count % 10_000 == 0 {
                log::info!("Processed {} / {} reads for SNPmer clustering in k-mer cluster {} with {} reps",  count, kmer_cluster.len(), kmer_cluster_id, representatives.len());
            }
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

        let mut local_clusters: Vec<Vec<usize>> = cluster_map.into_values()
            .map(|mut c| { c.sort_unstable(); c })
            .collect();
        local_clusters.sort_by(|a, b| b.len().cmp(&a.len())
            .then_with(|| a.first().cmp(&b.first())));

        // Remove small clusters
        log::debug!("Before size filtering: {} SNPmer clusters in k-mer cluster {}", local_clusters.len(), kmer_cluster_id);
        local_clusters.retain(|cluster| cluster.len() >= args.min_cluster_size);
        log::debug!("After size filtering: {} SNPmer clusters in k-mer cluster {}", local_clusters.len(), kmer_cluster_id);

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
    let cluster_file = output_dir.join("snpmer_clusters_before_reclust2.5.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&cluster_file).unwrap());
    writeln!(
        writer,
        "kmer_cluster_id\tsnpmer_cluster_id\tsize\trepresentative\tmembers"
    )
    .unwrap();

    let mut kmer_cluster_ids: Vec<usize> = local_clusters_map.keys().copied().collect();
    kmer_cluster_ids.sort_unstable();
    for kmer_cluster_id in kmer_cluster_ids {
        let snpmer_clusters = &local_clusters_map[&kmer_cluster_id];
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
                snpmer_cluster
                    .iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )
            .unwrap();
        }
    }

    log::info!("Wrote SNPmer clusters to {}", cluster_file.display());

    let recluster = true;
    let local_clusters_all: Vec<Vec<usize>>;
    if recluster {
        local_clusters_all = recluster_using_consensus_reps(local_clusters_map, twin_reads, args);
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

    let final_file = output_dir.join("final_snpmer_clusters_stage3.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&final_file).unwrap());

    for (snpmer_rep_id, cluster) in local_clusters_all.iter().enumerate() {
        let representative = cluster[0];
        writeln!(
            writer,
            "final_cluster_{}\tsize_{}\trepresentative_{}\tmembers\n{}",
            snpmer_rep_id,
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
        )
        .unwrap();
    }

    return local_clusters_all;
}

fn compare_blockmers(
    twin_read1: &TwinRead,
    twin_read2: &TwinRead,
    _k: usize,
    l: usize,
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
type Marker = (u64, Kmer48);

struct MarkerCluster {
    reads: Vec<usize>,
    consensus: Option<Vec<Marker>>,
}

fn for_each_marker(
    read: &TwinRead,
    k: usize,
    l: usize,
    marker_type: PolyMarkerType,
    mut visit: impl FnMut(Marker),
) {
    match marker_type {
        PolyMarkerType::Snpmer => {
            let mask = !(3 << (k - 1));
            for (_, kmer) in read.snpmers() {
                visit((kmer.to_u64() & mask, kmer));
            }
        }
        PolyMarkerType::Blockmer => {
            for (_, kmer) in read.blockmers() {
                visit((kmer >> (2 * l), Kmer48::from_u64(kmer)));
            }
        }
    }
}

fn build_marker_consensus(
    cluster: &[usize],
    twin_reads: &[TwinRead],
    k: usize,
    l: usize,
    marker_type: PolyMarkerType,
) -> Vec<Marker> {
    let mut counts: FxHashMap<Marker, usize> = FxHashMap::default();
    for &read_id in cluster {
        for_each_marker(&twin_reads[read_id], k, l, marker_type, |marker| {
            *counts.entry(marker).or_insert(0) += 1;
        });
    }

    let min_support = (cluster.len() / 6).max(1);
    let mut supported_alleles: Vec<_> = counts
        .into_iter()
        .filter(|&(_, count)| count >= min_support)
        .collect();
    // Put the most frequent allele first within each splitmer, breaking ties by k-mer.
    supported_alleles.sort_unstable_by_key(|&((splitmer, kmer), count)| {
        (splitmer, Reverse(count), kmer.to_u64())
    });
    supported_alleles.dedup_by_key(|((splitmer, _), _)| *splitmer);
    supported_alleles
        .into_iter()
        .map(|(marker, _)| marker)
        .collect()
}

/// Compare markers sorted by splitmer. The consensus has one allele per splitmer;
/// the query can repeat a splitmer, so score every occurrence before advancing the consensus.
#[inline]
fn compare_markers(
    query: &[Marker],
    mut consensus: &[Marker],
    max_mismatches: usize,
) -> (usize, usize) {
    let (mut matches, mut mismatches) = (0, 0);
    for &(query_splitmer, query_kmer) in query {
        loop {
            let Some((&(ref_splitmer, ref_kmer), rest)) = consensus.split_first() else {
                return (matches, mismatches);
            };
            if ref_splitmer < query_splitmer {
                consensus = rest;
                continue;
            }
            if ref_splitmer == query_splitmer {
                if query_kmer == ref_kmer {
                    matches += 1;
                } else {
                    mismatches += 1;
                    if mismatches > max_mismatches {
                        return (matches, mismatches);
                    }
                }
            }
            break;
        }
    }

    (matches, mismatches)
}

/// No consensus allele at this splitmer for this cluster.
const NO_ALLELE: u8 = 4;

/// The middle (variable) base of a SNPmer.
#[inline]
fn snpmer_allele(kmer: Kmer48, k: usize) -> u8 {
    ((kmer.to_u64() >> (k - 1)) & 3) as u8
}

/// Consensus SNPmer alleles of a group's clusters, laid out by splitmer: row `r` holds one
/// byte per cluster, the allele that cluster's consensus has at that splitmer (or
/// `NO_ALLELE`). Scoring a read against every cluster then streams one contiguous row per
/// read marker (vectorizable), instead of merge-scanning each cluster's consensus in turn.
struct AlleleMatrix {
    row_of: FxHashMap<u64, u32>,
    cells: Vec<u8>,
    num_clusters: usize,
}

impl AlleleMatrix {
    fn new<'a>(consensuses: impl Iterator<Item = Option<&'a [Marker]>>, num_clusters: usize, k: usize) -> Self {
        let mut row_of: FxHashMap<u64, u32> = FxHashMap::default();
        let mut cells: Vec<u8> = Vec::new();
        for (cluster_idx, consensus) in consensuses.enumerate() {
            for &(splitmer, kmer) in consensus.unwrap_or(&[]) {
                let next_row = row_of.len() as u32;
                let row = *row_of.entry(splitmer).or_insert(next_row);
                if row == next_row {
                    cells.resize(cells.len() + num_clusters, NO_ALLELE);
                }
                cells[row as usize * num_clusters + cluster_idx] = snpmer_allele(kmer, k);
            }
        }
        AlleleMatrix { row_of, cells, num_clusters }
    }

    #[inline]
    fn row(&self, splitmer: u64) -> Option<&[u8]> {
        self.row_of.get(&splitmer).map(|&row| {
            let start = row as usize * self.num_clusters;
            &self.cells[start..start + self.num_clusters]
        })
    }
}

/// Adds one marker with `allele` to the (matches, mismatches) of every cluster in `row`.
#[inline]
fn score_row(row: &[u8], allele: u8, matches: &mut [u16], mismatches: &mut [u16]) {
    for ((m, mm), &a) in matches.iter_mut().zip(mismatches.iter_mut()).zip(row) {
        *m += (a == allele) as u16;
        *mm += (a != allele && a != NO_ALLELE) as u16;
    }
}

/// Reassign reads to their best matching cluster based on polymorphic markers.
fn reassign_reads_to_best_cluster(
    mut current_clusters: Vec<MarkerCluster>,
    twin_reads: &[TwinRead],
    k: usize,
    l: usize,
    marker_type: PolyMarkerType,
    min_cluster_size: usize,
) -> (Vec<Vec<usize>>, usize) {
    if current_clusters.len() < 2 {
        let clusters = current_clusters
            .into_iter()
            .map(|cluster| cluster.reads)
            .filter(|reads| !reads.is_empty() && reads.len() >= min_cluster_size)
            .collect();
        return (clusters, 0);
    }

    current_clusters
        .par_iter_mut()
        .filter(|cluster| cluster.consensus.is_none())
        .for_each(|cluster| {
            cluster.consensus = Some(build_marker_consensus(
                &cluster.reads,
                twin_reads,
                k,
                l,
                marker_type,
            ));
        });

    // The matrix scorer counts in u16; a read with more markers than that (only possible with
    // a very large --max-read-length) takes the scan below instead.
    let counts_fit = || {
        current_clusters
            .iter()
            .flat_map(|cluster| &cluster.reads)
            .all(|&read_id| twin_reads[read_id].snpmer_count() <= u16::MAX as usize)
    };
    if marker_type == PolyMarkerType::Snpmer && counts_fit() {
        return reassign_snpmer_reads_by_matrix(current_clusters, twin_reads, k, min_cluster_size);
    }

    // Collect assignments per source cluster so we can distribute reads without a shared lock.
    let per_cluster: Vec<_> = current_clusters
        .par_iter()
        .enumerate()
        .map(|(cluster_idx, cluster)| {
            let mut assignments = Vec::with_capacity(cluster.reads.len());
            let mut read_markers = Vec::new();
            let mut reassigned = 0;
            for &read_id in &cluster.reads {
                read_markers.clear();
                for_each_marker(&twin_reads[read_id], k, l, marker_type, |marker| {
                    read_markers.push(marker);
                });
                read_markers.sort_unstable_by_key(|&(splitmer, _)| splitmer);

                let mut best_cluster = cluster_idx;
                let (mut best_mismatches, mut best_matches) = (usize::MAX, 0);
                for (candidate_idx, candidate) in current_clusters.iter().enumerate() {
                    let (matches, mismatches) = compare_markers(
                        &read_markers,
                        candidate.consensus.as_deref().unwrap(),
                        best_mismatches,
                    );
                    if mismatches < best_mismatches
                        || (mismatches == best_mismatches && matches > best_matches)
                    {
                        best_mismatches = mismatches;
                        best_matches = matches;
                        best_cluster = candidate_idx;
                    }
                }

                assignments.push((best_cluster, read_id));
                if best_cluster != cluster_idx {
                    reassigned += 1;
                    log::trace!(
                        "Reassigned read {} from cluster {} to cluster {} (mismatches: {}, matches: {})",
                        read_id,
                        cluster_idx,
                        best_cluster,
                        best_mismatches,
                        best_matches
                    );
                }
            }
            (assignments, reassigned)
        })
        .collect();

    let mut new_clusters: Vec<Vec<usize>> = vec![Vec::new(); current_clusters.len()];
    let mut num_reassignments = 0;
    for (assignments, reassigned) in per_cluster {
        num_reassignments += reassigned;
        for (best_cluster, read_id) in assignments {
            new_clusters[best_cluster].push(read_id);
        }
    }

    new_clusters.retain(|cluster| !cluster.is_empty() && cluster.len() >= min_cluster_size);
    for cluster in &mut new_clusters {
        cluster.sort_unstable();
    }
    log::trace!(
        "Reassignment complete: {} reads reassigned to better clusters",
        num_reassignments
    );

    (new_clusters, num_reassignments)
}

/// Same result as the merge-scan path in `reassign_reads_to_best_cluster`, for SNPmers:
/// each read is scored against all clusters at once through an `AlleleMatrix`, and the work
/// is split over reads rather than clusters so one big cluster is not left to one thread.
fn reassign_snpmer_reads_by_matrix(
    current_clusters: Vec<MarkerCluster>,
    twin_reads: &[TwinRead],
    k: usize,
    min_cluster_size: usize,
) -> (Vec<Vec<usize>>, usize) {
    let num_clusters = current_clusters.len();
    let mask = !(3u64 << (k - 1));
    let matrix = AlleleMatrix::new(
        current_clusters.iter().map(|cluster| cluster.consensus.as_deref()),
        num_clusters,
        k,
    );
    let jobs: Vec<(usize, usize)> = current_clusters
        .iter()
        .enumerate()
        .flat_map(|(cluster_idx, cluster)| cluster.reads.iter().map(move |&read_id| (cluster_idx, read_id)))
        .collect();

    let best_clusters: Vec<usize> = jobs
        .par_iter()
        .with_min_len(64)
        .map_init(
            || (vec![0u16; num_clusters], vec![0u16; num_clusters]),
            |(matches, mismatches), &(cluster_idx, read_id)| {
                let read = &twin_reads[read_id];
                matches.fill(0);
                mismatches.fill(0);
                for (_, kmer) in read.snpmers() {
                    if let Some(row) = matrix.row(kmer.to_u64() & mask) {
                        score_row(row, snpmer_allele(kmer, k), matches, mismatches);
                    }
                }
                // The scan's rule: fewest mismatches, then most matches, then lowest index.
                // (The scan's early exits only skip clusters that cannot win.)
                let mut best_cluster = cluster_idx;
                let (mut best_mismatches, mut best_matches) = (u32::MAX, 0);
                for candidate_idx in 0..num_clusters {
                    let (m, mm) = (matches[candidate_idx] as u32, mismatches[candidate_idx] as u32);
                    if mm < best_mismatches || (mm == best_mismatches && m > best_matches) {
                        best_mismatches = mm;
                        best_matches = m;
                        best_cluster = candidate_idx;
                    }
                }
                best_cluster
            },
        )
        .collect();

    let mut new_clusters: Vec<Vec<usize>> = vec![Vec::new(); num_clusters];
    let mut num_reassignments = 0;
    for (&(cluster_idx, read_id), &best_cluster) in jobs.iter().zip(&best_clusters) {
        if best_cluster != cluster_idx {
            num_reassignments += 1;
        }
        new_clusters[best_cluster].push(read_id);
    }
    new_clusters.retain(|cluster| !cluster.is_empty() && cluster.len() >= min_cluster_size);
    for cluster in &mut new_clusters {
        cluster.sort_unstable();
    }
    (new_clusters, num_reassignments)
}

/// Merge concordant clusters in order of descending depth.
fn recluster_one_round(
    current_clusters: Vec<Vec<usize>>,
    twin_reads: &[TwinRead],
    k: usize,
    l: usize,
    marker_type: PolyMarkerType,
) -> (Vec<MarkerCluster>, usize) {
    let mut clusters: Vec<_> = current_clusters
        .into_iter()
        .filter(|reads| !reads.is_empty())
        .map(|reads| MarkerCluster {
            reads,
            consensus: None,
        })
        .collect();
    if clusters.len() < 2 {
        for cluster in &mut clusters {
            cluster.reads.sort_unstable();
        }
        return (clusters, 0);
    }

    clusters.sort_by(|a, b| {
        b.reads
            .len()
            .cmp(&a.reads.len())
            .then_with(|| a.reads.first().cmp(&b.reads.first()))
    });
    clusters.par_iter_mut().for_each(|cluster| {
        cluster.consensus = Some(build_marker_consensus(
            &cluster.reads,
            twin_reads,
            k,
            l,
            marker_type,
        ));
    });

    // For SNPmers, score each cluster's consensus against all later ones at once.
    let num_clusters = clusters.len();
    // u16 pair counts: a consensus has one marker per splitmer, so this never binds in
    // practice; the scan below covers it if it ever does.
    let counts_fit = clusters
        .iter()
        .all(|c| c.consensus.as_ref().map_or(true, |m| m.len() <= u16::MAX as usize));
    let matrix = (marker_type == PolyMarkerType::Snpmer && counts_fit).then(|| {
        AlleleMatrix::new(clusters.iter().map(|c| c.consensus.as_deref()), num_clusters, k)
    });
    let mut pair_matches = vec![0u16; num_clusters];
    let mut pair_mismatches = vec![0u16; num_clusters];

    let mut num_merges = 0;
    for i in 0..clusters.len() {
        let (left, right) = clusters.split_at_mut(i + 1);
        let cluster = &mut left[i];
        if cluster.reads.is_empty() {
            continue;
        }

        let original_size = cluster.reads.len();
        let consensus = cluster.consensus.as_deref().unwrap();
        if let Some(matrix) = &matrix {
            // This cluster's consensus is fixed during its scan, so all its pair scores can
            // be computed up front.
            let (m, mm) = (&mut pair_matches[i + 1..], &mut pair_mismatches[i + 1..]);
            m.fill(0);
            mm.fill(0);
            for &(splitmer, kmer) in consensus {
                let row = matrix.row(splitmer).unwrap();
                score_row(&row[i + 1..], snpmer_allele(kmer, k), m, mm);
            }
        }
        for (j, candidate) in (i + 1..).zip(right) {
            if candidate.reads.is_empty() {
                continue;
            }

            let other_consensus = candidate.consensus.as_deref().unwrap();
            // With mismatches > 0 the pair is never concordant, so the scan's early exit
            // (partial match count) and the full count give the same decision.
            let (matches, mismatches) = if matrix.is_some() {
                (pair_matches[j] as usize, pair_mismatches[j] as usize)
            } else {
                compare_markers(consensus, other_consensus, 0)
            };
            let len = consensus.len();
            let other_len = other_consensus.len();
            let larger = cluster.reads.len().max(candidate.reads.len());
            let smaller = cluster.reads.len().min(candidate.reads.len());
            let concordant = mismatches == 0
                && ((matches >= len.min(other_len.max(2)) && matches >= other_len.min(len.max(2)))
                    || ((matches as f64) > len.min(other_len) as f64 * 0.975
                        && larger / smaller > 50)
                    || (larger / smaller > 500 && smaller <= 2));

            if concordant {
                let old_size = cluster.reads.len();
                let candidate_size = candidate.reads.len();
                cluster.reads.extend(std::mem::take(&mut candidate.reads));
                candidate.consensus = None;
                num_merges += 1;
                log::trace!(
                    "Merged cluster {} (size {}) into cluster {} (old size: {}, new size: {})",
                    j,
                    candidate_size,
                    i,
                    old_size,
                    cluster.reads.len()
                );
            } else {
                log::trace!("Clusters {} and {} not concordant, not merging", i, j);
            }
        }

        // Keep the original markers throughout the merge scan. Reassignment can rebuild
        // them if this cluster gained reads and there is still more than one destination.
        if cluster.reads.len() != original_size {
            cluster.consensus = None;
        }
        cluster.reads.sort_unstable();
    }

    clusters.retain(|cluster| !cluster.reads.is_empty());
    clusters.sort_by(|a, b| {
        b.reads
            .len()
            .cmp(&a.reads.len())
            .then_with(|| a.reads.first().cmp(&b.reads.first()))
    });

    (clusters, num_merges)
}

pub fn recluster_using_consensus_reps(
    clusters: FxHashMap<usize, Vec<Vec<usize>>>,
    twin_reads: &[TwinRead],
    args: &Cli,
) -> Vec<Vec<usize>> {
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
    log::info!(
        "Starting with {} k-mer groups containing {} total SNPmer clusters",
        total_groups,
        total_initial_clusters
    );

    // Step 2: Iteratively recluster until convergence
    let mut iteration = 0;
    // Groups whose last iteration had no merges and no reassignments. Such a group is at a
    // fixed point: its clusters, and therefore its consensuses, are unchanged, so another
    // round would reproduce it exactly. They are carried over instead of being recomputed.
    let mut settled: FxHashSet<usize> = FxHashSet::default();
    loop {
        if iteration >= args.max_iterations_recluster {
            log::info!(
                "Reached maximum reclustering iterations ({})",
                args.max_iterations_recluster
            );
            break;
        }
        iteration += 1;
        let total_merges = Mutex::new(0);
        let total_reassignments = Mutex::new(0);

        let new_clusters: Mutex<FxHashMap<usize, Vec<Vec<usize>>>> =
            Mutex::new(FxHashMap::default());
        let newly_settled: Mutex<Vec<usize>> = Mutex::new(Vec::new());
        let settled_ref = &settled;

        // Process each k-mer group independently
        current_clusters
            .into_par_iter()
            .for_each(|(kmer_cluster_id, snpmer_clusters)| {
                if settled_ref.contains(&kmer_cluster_id) {
                    new_clusters.lock().unwrap().insert(kmer_cluster_id, snpmer_clusters);
                    return;
                }
                log::trace!(
                    "Processing k-mer group {} with {} SNPmer clusters",
                    kmer_cluster_id,
                    snpmer_clusters.len()
                );

                // Step 2a: Merge clusters within this group
                let (merged_clusters, num_merges) = recluster_one_round(
                    snpmer_clusters,
                    twin_reads,
                    k,
                    args.blockmer_length,
                    marker_type,
                );

                *total_merges.lock().unwrap() += num_merges;
                log::trace!("Processing k-mer group {} --- merged", kmer_cluster_id);

                // Step 2b: Reassign reads to best matching clusters within this group
                let (reassigned_clusters, num_reassignments) = reassign_reads_to_best_cluster(
                    merged_clusters,
                    twin_reads,
                    k,
                    args.blockmer_length,
                    marker_type,
                    args.min_cluster_size,
                );
                *total_reassignments.lock().unwrap() += num_reassignments;
                if num_merges == 0 && num_reassignments == 0 {
                    newly_settled.lock().unwrap().push(kmer_cluster_id);
                }

                if !reassigned_clusters.is_empty() {
                    new_clusters
                        .lock()
                        .unwrap()
                        .insert(kmer_cluster_id, reassigned_clusters);
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
        settled.extend(newly_settled.into_inner().unwrap());

        // Check for convergence (no merges and no reassignments across all groups)
        if total_merges == 0 {
            log::info!("Convergence reached after {} iterations", iteration);
            break;
        }
    }

    // Step 3: Flatten the hierarchical structure for final output and debugging
    let final_group_count = current_clusters.len();
    let mut final_clusters: Vec<Vec<usize>> = Vec::new();
    let mut sorted_keys: Vec<usize> = current_clusters.keys().cloned().collect();
    sorted_keys.sort_unstable();
    for kmer_cluster_id in sorted_keys {
        let snpmer_clusters = current_clusters.remove(&kmer_cluster_id).unwrap();
        final_clusters.extend(
            snpmer_clusters
                .into_iter()
                .filter(|cluster| !cluster.is_empty()),
        );
    }

    // Sort by size descending
    final_clusters.sort_by(|a, b| {
        b.len()
            .cmp(&a.len())
            .then_with(|| a.first().cmp(&b.first()))
    });
    final_clusters.retain(|cluster| cluster.len() >= args.min_cluster_size);

    log::info!(
        "Final result: {} total clusters across {} k-mer groups",
        final_clusters.len(),
        final_group_count
    );

    // Step 4: Debugging - Compare all final clusters and output mismatch statistics
    if log::log_enabled!(log::Level::Trace) {
        log::info!("Performing detailed pairwise cluster comparison for debugging...");
        log::info!("=== Debugging: Pairwise cluster comparison ===");

        // Build consensus for each final cluster
        let final_consensus: Vec<Vec<Marker>> = final_clusters
            .iter()
            .map(|cluster| {
                build_marker_consensus(
                    cluster,
                    twin_reads,
                    k,
                    args.blockmer_length,
                    PolyMarkerType::Snpmer,
                )
            })
            .collect();

        // Compare all pairs of clusters
        for i in 0..final_clusters.len() {
            let rep_i = final_clusters[i][0]; // First member as representative

            for j in (i + 1)..final_clusters.len() {
                let rep_j = final_clusters[j][0];

                // Count matches and mismatches
                let (matches, mismatches) =
                    compare_markers(&final_consensus[i], &final_consensus[j], usize::MAX);

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn marker_scan_counts_repeated_read_markers_and_can_stop_on_mismatch() {
        let marker = |splitmer, allele| (splitmer, Kmer48::from_u64(allele));
        let query = [
            marker(1, 10),
            marker(2, 20),
            marker(2, 20),
            marker(2, 21),
            marker(4, 40),
        ];
        let consensus = [marker(2, 20), marker(3, 30), marker(4, 40)];

        assert_eq!(compare_markers(&query, &consensus, usize::MAX), (3, 1));
        assert_eq!(compare_markers(&query, &consensus, 1), (3, 1));
        assert_eq!(compare_markers(&query, &consensus, 0), (2, 1));
        assert_eq!(compare_markers(&query, &[marker(5, 50)], 0), (0, 0));
    }

    #[test]
    fn allele_matrix_matches_merge_scan() {
        let k = 17;
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut rand = move |n: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % n
        };
        let mid = 3u64 << (k - 1);
        // A small splitmer universe so clusters and reads overlap heavily.
        let universe: Vec<u64> = (0..60).map(|i| (i * 7_919_003) & !mid & ((1 << 34) - 1)).collect();
        let marker = |splitmer: u64, allele: u64| (splitmer, Kmer48::from_u64(splitmer | (allele << (k - 1))));
        for _trial in 0..300 {
            let num_clusters = 1 + rand(12) as usize;
            let mut consensuses: Vec<Vec<Marker>> = Vec::new();
            for _ in 0..num_clusters {
                let mut c: Vec<Marker> = Vec::new();
                for &s in &universe {
                    if rand(3) > 0 {
                        c.push(marker(s, rand(4)));
                    }
                }
                c.sort_unstable_by_key(|&(s, _)| s);
                consensuses.push(c);
            }
            let matrix = AlleleMatrix::new(consensuses.iter().map(|c| Some(c.as_slice())), num_clusters, k);
            for _read in 0..20 {
                let mut read: Vec<Marker> = Vec::new();
                for &s in &universe {
                    if rand(2) > 0 {
                        read.push(marker(s, rand(4)));
                    }
                }
                read.sort_unstable_by_key(|&(s, _)| s);
                let own = rand(num_clusters as u64) as usize;
                // Original scan.
                let mut scan_best = own;
                let (mut bm, mut bmatch) = (usize::MAX, 0);
                for (j, c) in consensuses.iter().enumerate() {
                    let (m, mm) = compare_markers(&read, c, bm);
                    if mm < bm || (mm == bm && m > bmatch) {
                        bm = mm;
                        bmatch = m;
                        scan_best = j;
                    }
                }
                // Matrix.
                let (mut m, mut mm) = (vec![0u16; num_clusters], vec![0u16; num_clusters]);
                for &(s, kmer) in &read {
                    if let Some(row) = matrix.row(s) {
                        score_row(row, snpmer_allele(kmer, k), &mut m, &mut mm);
                    }
                }
                let mut best = own;
                let (mut bm2, mut bmatch2) = (u32::MAX, 0);
                for j in 0..num_clusters {
                    if (mm[j] as u32) < bm2 || (mm[j] as u32 == bm2 && m[j] as u32 > bmatch2) {
                        bm2 = mm[j] as u32;
                        bmatch2 = m[j] as u32;
                        best = j;
                    }
                }
                assert_eq!(scan_best, best);
                // Full counts agree too.
                for j in 0..num_clusters {
                    assert_eq!(compare_markers(&read, &consensuses[j], usize::MAX), (m[j] as usize, mm[j] as usize));
                }
            }
            // Merge-round pair scores: consensus i against all clusters.
            for i in 0..num_clusters {
                let (mut m, mut mm) = (vec![0u16; num_clusters], vec![0u16; num_clusters]);
                for &(s, kmer) in &consensuses[i] {
                    score_row(matrix.row(s).unwrap(), snpmer_allele(kmer, k), &mut m, &mut mm);
                }
                for j in 0..num_clusters {
                    let (sm, smm) = compare_markers(&consensuses[i], &consensuses[j], 0);
                    assert_eq!(smm == 0, mm[j] == 0);
                    if smm == 0 {
                        assert_eq!(sm, m[j] as usize);
                    }
                }
            }
        }
    }

    #[test]
    fn zero_or_one_cluster_needs_no_read_markers() {
        let (clusters, merges) = recluster_one_round(
            vec![vec![3, 1, 2], vec![]],
            &[],
            17,
            3,
            PolyMarkerType::Snpmer,
        );
        assert_eq!(merges, 0);
        assert!(clusters[0].consensus.is_none());
        assert_eq!(
            reassign_reads_to_best_cluster(clusters, &[], 17, 3, PolyMarkerType::Snpmer, 3),
            (vec![vec![1, 2, 3]], 0)
        );

        let (clusters, _) = recluster_one_round(vec![vec![]], &[], 17, 3, PolyMarkerType::Snpmer);
        assert_eq!(
            reassign_reads_to_best_cluster(clusters, &[], 17, 3, PolyMarkerType::Snpmer, 3),
            (Vec::new(), 0)
        );
    }
}
