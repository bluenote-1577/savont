use crate::cli::Cli;
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

/// Compare two reads by SNPmers and count matches/mismatches
/// Returns (matches, mismatches)
fn compare_snpmers(
    read1_snpmers: &[(u32, Kmer48)],
    read2_snpmers: &[(u32, Kmer48)],
    k: usize,
) -> (usize, usize) {
    let mask = !(3 << (k - 1));

    // Build index of read2's splitmers (k-mer without middle base)
    let mut splitmer_index: FxHashMap<u64, Vec<(usize, Kmer48)>> = FxHashMap::default();
    for (idx, &(_pos, kmer)) in read2_snpmers.iter().enumerate() {
        let splitmer = kmer.to_u64() & mask;
        splitmer_index.entry(splitmer).or_insert_with(Vec::new).push((idx, kmer));
    }

    let mut matches = 0;
    let mut mismatches = 0;

    // For each SNPmer in read1, check if read2 has matching splitmer
    for &(_pos, kmer1) in read1_snpmers {
        let splitmer1 = kmer1.to_u64() & mask;

        if let Some(candidates) = splitmer_index.get(&splitmer1) {
            for &(_idx, kmer2) in candidates {
                // Splitmer matches, now check middle base
                if kmer1 == kmer2 {
                    matches += 1;
                } else {
                    mismatches += 1;
                }
            }
        }
    }

    (matches, mismatches)
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
/// Returns true if any mismatches found (not compatible)
fn all_snpmer_mismatch(
    index: &FxHashMap<u64, Vec<(usize, Kmer48)>>,
    query_snpmers: &[(u32, Kmer48)],
    k: usize,
) -> bool {
    let mask = !(3 << (k - 1));
    let mut clusters_with_mismatch = index.values().flat_map(|v| v.iter().map(|&(id, _)| (id, false))).collect::<FxHashMap<usize, bool>>();

    for &(_pos, query_kmer) in query_snpmers {
        let query_splitmer = query_kmer.to_u64() & mask;

        if let Some(candidates) = index.get(&query_splitmer) {
            for &(candidate_id, candidate_kmer) in candidates {
                // Splitmer matches - check if middle base differs
                if query_kmer != candidate_kmer {
                    clusters_with_mismatch.insert(candidate_id, true);
                }
            }
        }
    }

    // If any candidate has no mismatches, return false
    return clusters_with_mismatch.values().all(|&mismatch| mismatch);
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

    // Global cluster assignment: read_id -> snpmer_representative_read_id
    let mut snpmer_cluster_assignment: HashMap<usize, usize> = HashMap::new();
    let mut local_clusters_all: Vec<Vec<usize>> = Vec::new();

    let cluster_file = output_dir.join("snpmer_clusters.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&cluster_file).unwrap());

    writeln!(writer, "kmer_cluster_id\tsnpmer_cluster_id\tsize\trepresentative\tmembers").unwrap();

    let mut total_snpmer_clusters = 0;

    // Process each k-mer cluster independently
    for (kmer_cluster_id, kmer_cluster) in kmer_clusters.iter().enumerate() {
        if kmer_cluster.len() < 1 {
            continue;
        }

        // SNPmer inverted index for this k-mer cluster: splitmer -> Vec<(read_id, full_kmer)>
        let mut snpmer_index: FxHashMap<u64, Vec<(usize, Kmer48)>> = FxHashMap::default();

        // Local cluster assignments within this k-mer cluster
        let mut local_assignment: HashMap<usize, usize> = HashMap::new();

        // Process reads in this k-mer cluster sequentially
        for &read_id in kmer_cluster {
            let read_snpmers = twin_reads[read_id].snpmers_vec();

            // Query this read's SNPmers against the current index
            let all_mismatch = all_snpmer_mismatch(&snpmer_index, &read_snpmers, k);

            if all_mismatch {
                // Mismatch found - this read becomes a new SNPmer representative
                add_read_snpmers_to_index(&mut snpmer_index, read_id, &read_snpmers, k);
                local_assignment.insert(read_id, read_id);
            } else {
                // Some read exists without mismatches - assign to existing representative
                // Find the representative by checking which read is in the index
                let mut assigned = false;
                let mut representative_stats = vec![];
                for &other_read_id in kmer_cluster {
                    if local_assignment.get(&other_read_id) == Some(&other_read_id) {
                        // This is a representative - check if compatible
                        let other_snpmers = twin_reads[other_read_id].snpmers_vec();
                        let (matches, mismatches) = compare_snpmers(&read_snpmers, &other_snpmers, k);

                        log::trace!("Comparing read {} to representative {}: {} matches, {} mismatches", read_id, other_read_id, matches, mismatches);

                        if mismatches == 0 && matches > 0 {
                            representative_stats.push((mismatches, - (matches as i64), other_read_id));
                        }
                    }

                    
                }

                if !representative_stats.is_empty() {
                    // Choose the best representative (most matches, least mismatches)
                    representative_stats.sort();
                    let best_rep = representative_stats[0].2;
                    if representative_stats.len() > 1 && representative_stats[0].1 == representative_stats[1].1 {
                        continue;
                    }
                    local_assignment.insert(read_id, best_rep);
                    assigned = true;
                }
            }
        }

        // Build SNPmer subclusters within this k-mer cluster
        let mut local_clusters_map: HashMap<usize, Vec<usize>> = HashMap::new();
        for (read_id, rep_id) in &local_assignment {
            if read_id != rep_id {
                continue;
            }
            local_clusters_map.entry(*rep_id).or_insert_with(Vec::new).push(*rep_id);
            snpmer_cluster_assignment.insert(*read_id, *rep_id);
        }
        for (read_id, rep_id) in local_assignment.iter() {
            if read_id == rep_id {
                continue;
            }
            local_clusters_map.entry(*rep_id).or_insert_with(Vec::new).push(*read_id);
            snpmer_cluster_assignment.insert(*read_id, *rep_id);
        }

        let mut local_clusters: Vec<Vec<usize>> = local_clusters_map.into_values().collect();
        local_clusters.sort_by(|a, b| b.len().cmp(&a.len()));

        local_clusters.retain(|cluster| cluster.len() >= args.min_cluster_size);

        // Write SNPmer subclusters for this k-mer cluster
        for (local_snpmer_id, snpmer_cluster) in local_clusters.iter().enumerate() {
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
            total_snpmer_clusters += 1;
        }

        local_clusters_all.extend(local_clusters);

        if kmer_cluster_id % 100 == 0 && kmer_cluster_id > 0 {
            log::debug!(
                "Processed {} / {} k-mer clusters for SNPmer clustering",
                kmer_cluster_id,
                kmer_clusters.len()
            );
        }
    }

    log::info!(
        "SNPmer clustering complete. {} total SNPmer clusters from {} k-mer clusters",
        total_snpmer_clusters,
        kmer_clusters.len()
    );
    log::info!("Wrote SNPmer clusters to {}", cluster_file.display());

    let final_file = output_dir.join("final_clusters.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&final_file).unwrap());

    for (snpmer_rep_id, cluster) in local_clusters_all.iter().enumerate() {
        let representative = cluster[0];
        writeln!(
            writer,
            "final_cluster_{}\tsize_{}\trepresentative_{}\tmembers_{}",
            snpmer_rep_id,
            cluster.len(),
            representative,
            cluster.iter().map(|x| format!("{} {}", &twin_reads[*x].id, &twin_reads[*x].est_id.unwrap())).collect::<Vec<_>>().join("\n")
        ).unwrap();
    }

    return local_clusters_all;
}
