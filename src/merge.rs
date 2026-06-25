use crate::cli;
use crate::constants::ASV_FILE;
use crate::seeding::minimizer_seeds_positions;
use crate::taxonomy::load_fasta_with_needletail;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

// ── sequence hashing ─────────────────────────────────────────────────────────

fn djb2_hash(seq: &[u8]) -> u64 {
    let mut h: u64 = 5381;
    for &b in seq {
        h = h.wrapping_mul(33).wrapping_add(b.to_ascii_uppercase() as u64);
    }
    h
}

fn seq_hash(seq: &[u8]) -> String {
    let fwd = djb2_hash(seq);
    let rc_seq = crate::utils::reverse_complement(seq);
    let rev = djb2_hash(&rc_seq);
    format!("{:016x}", fwd.min(rev))
}

// ── helpers ───────────────────────────────────────────────────────────────────

fn sample_name_from_dir(dir: &Path) -> String {
    let ft = dir.join("feature-table.tsv");
    if let Ok(f) = std::fs::File::open(&ft) {
        for line in BufReader::new(f).lines().flatten() {
            if line.starts_with("#OTU ID") {
                if let Some(name) = line.split('\t').nth(1) {
                    return name.to_string();
                }
            }
        }
    }
    dir.file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("sample")
        .to_string()
}

/// Read feature-table.tsv and return (sample_names, otu_id → per-sample depths).
/// Handles both single-column and multi-column (pooled-samples) tables.
fn feature_table_from_dir(dir: &Path) -> Option<(Vec<String>, HashMap<String, Vec<u64>>)> {
    let ft = dir.join("feature-table.tsv");
    let f = std::fs::File::open(&ft).ok()?;
    let mut lines = BufReader::new(f).lines().flatten();

    // Skip comment lines until we hit the header
    let header_line = lines.find(|l| l.starts_with("#OTU ID"))?;
    let cols: Vec<&str> = header_line.split('\t').collect();
    // Columns 1.. are sample names
    let sample_names: Vec<String> = cols[1..].iter().map(|s| s.to_string()).collect();
    if sample_names.is_empty() {
        return None;
    }

    let n_samples = sample_names.len();
    let mut depths: HashMap<String, Vec<u64>> = HashMap::new();
    for line in lines {
        if line.is_empty() || line.starts_with('#') { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() { continue; }
        let otu_id = fields[0].to_string();
        let per_sample: Vec<u64> = (1..=n_samples)
            .map(|i| fields.get(i).and_then(|v| v.parse().ok()).unwrap_or(0))
            .collect();
        depths.insert(otu_id, per_sample);
    }

    Some((sample_names, depths))
}

fn depth_from_header_total(header: &str) -> u64 {
    let token = header.split_whitespace().next().unwrap_or("").split('_').last().unwrap_or("0");
    // Sum dash-separated per-sample depths if present (pooled format)
    token.split('-').filter_map(|s| s.parse::<u64>().ok()).sum::<u64>().max(0)
}

// ── taxonomy TSV parsing ──────────────────────────────────────────────────────

/// Read `asv_mappings.tsv` and return `(asv_header, qiime_lineage)` pairs.
/// The lineage is built from taxonomy columns in QIIME2 order
/// (superkingdom;phylum;class;order;family;genus;species).  Missing columns
/// are skipped so the function works with both classify and sintax mapping files.
fn read_asv_mapping_keys(path: &Path) -> std::io::Result<Vec<(String, String)>> {
    let f = std::fs::File::open(path)?;
    let mut lines = BufReader::new(f).lines();

    let header_line = lines
        .next()
        .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "empty file"))??;
    let cols: Vec<&str> = header_line.split('\t').collect();

    const QIIME_ORDER: &[&str] = &[
        "superkingdom", "phylum", "class", "order", "family", "genus", "species",
    ];
    let level_indices: Vec<Option<usize>> = QIIME_ORDER
        .iter()
        .map(|name| cols.iter().position(|&c| c == *name))
        .collect();

    let mut pairs = Vec::new();
    for line in lines.flatten() {
        if line.is_empty() { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() { continue; }
        let asv_header = fields[0].to_string();
        let lineage = level_indices.iter()
            .filter_map(|opt| opt.and_then(|i| fields.get(i)).map(|v| v.to_string()))
            .collect::<Vec<_>>()
            .join(";");
        pairs.push((asv_header, lineage));
    }
    Ok(pairs)
}

// ── writers ───────────────────────────────────────────────────────────────────

fn write_feature_table(
    table: &BTreeMap<String, (Vec<u8>, Vec<u64>)>,
    sample_names: &[String],
    path: &Path,
) -> std::io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    write!(f, "#OTU ID")?;
    for s in sample_names { write!(f, "\t{}", s)?; }
    writeln!(f)?;
    for (hash, (_, counts)) in table {
        write!(f, "{}", hash)?;
        for &c in counts { write!(f, "\t{}", c)?; }
        writeln!(f)?;
    }
    Ok(())
}

fn write_rep_seqs(
    table: &BTreeMap<String, (Vec<u8>, Vec<u64>)>,
    path: &Path,
) -> std::io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    for (hash, (seq, _)) in table {
        writeln!(f, ">{}", hash)?;
        writeln!(f, "{}", String::from_utf8_lossy(seq))?;
    }
    Ok(())
}

/// Write a transposed taxon count table: rows = samples, columns = unique lineage strings.
/// Counts are summed across all ASVs that share the same lineage.
fn write_taxon_table(
    asv_table: &BTreeMap<String, (Vec<u8>, Vec<u64>)>,
    hash_to_lineage: &HashMap<String, String>,
    sample_names: &[String],
    path: &Path,
) -> std::io::Result<()> {
    let mut lineage_counts: BTreeMap<String, Vec<u64>> = BTreeMap::new();
    for (hash, (_, counts)) in asv_table {
        let lineage = hash_to_lineage
            .get(hash)
            .map(|s| s.as_str())
            .unwrap_or("Unclassified")
            .to_string();
        let entry = lineage_counts
            .entry(lineage)
            .or_insert_with(|| vec![0u64; sample_names.len()]);
        for (a, &b) in entry.iter_mut().zip(counts) {
            *a += b;
        }
    }

    if lineage_counts.is_empty() {
        return Ok(());
    }

    let mut f = std::fs::File::create(path)?;

    write!(f, "taxon")?;
    for name in sample_names {
        write!(f, "\t{}", name)?;
    }
    writeln!(f)?;

    for (lineage, counts) in &lineage_counts {
        write!(f, "{}", lineage)?;
        for &c in counts {
            write!(f, "\t{}", c)?;
        }
        writeln!(f)?;
    }
    Ok(())
}

/// Write ASV-level taxonomy (hash IDs as Feature IDs) — for use with
/// the merged_feature_table in qiime taxa barplot.
/// Every hash in `asv_table` is written; hashes without a classification
/// get "Unclassified" so QIIME2 does not complain about missing IDs.
fn write_asv_taxonomy_file(
    asv_table: &BTreeMap<String, (Vec<u8>, Vec<u64>)>,
    hash_to_lineage: &HashMap<String, String>,
    path: &Path,
) -> std::io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "Feature ID\tTaxon")?;
    for hash in asv_table.keys() {
        let lineage = hash_to_lineage.get(hash).map(|s| s.as_str()).unwrap_or("Unclassified");
        writeln!(f, "{}\t{}", hash, lineage)?;
    }
    Ok(())
}

// ── fuzzy merge ───────────────────────────────────────────────────────────────

fn compute_minimizers(seq: &[u8]) -> Vec<u64> {
    let mut kmer_vec = Vec::new();
    let mut pos_vec = Vec::new();
    minimizer_seeds_positions(seq, &mut kmer_vec, &mut pos_vec, 28, 31);
    kmer_vec.sort_unstable();
    kmer_vec.dedup();
    kmer_vec
}

/// Merge near-identical ASVs (within MAX_LEN_DIFF bp) that share all minimizers
/// of the shorter sequence. Processes shortest-first so shorter sequences are
/// absorbed into longer representatives. Returns the number of ASVs absorbed.
fn fuzzy_merge_table(
    table: &mut BTreeMap<String, (Vec<u8>, Vec<u64>)>,
    hash_to_lineage: &mut HashMap<String, String>,
) -> usize {
    const MAX_LEN_DIFF: usize = 10;

    // Compute unique minimizers per sequence
    let minimizers: HashMap<String, Vec<u64>> = table.keys()
        .map(|h| (h.clone(), compute_minimizers(&table[h].0)))
        .collect();

    // Inverted index: minimizer kmer → set of hashes containing it
    let mut inverted: HashMap<u64, HashSet<String>> = HashMap::new();
    for (hash, kmers) in &minimizers {
        for &kmer in kmers {
            inverted.entry(kmer).or_default().insert(hash.clone());
        }
    }

    // Process shortest sequences first so they get absorbed into longer ones
    let mut sorted_hashes: Vec<String> = table.keys().cloned().collect();
    sorted_hashes.sort_by_key(|h| table[h].0.len());

    let n_before = sorted_hashes.len();
    let mut absorbed: HashSet<String> = HashSet::new();

    for hash in &sorted_hashes {
        if absorbed.contains(hash) { continue; }
        let kmers = &minimizers[hash];
        if kmers.is_empty() { continue; }
        let seq_len = table[hash].0.len();

        // Intersect candidate sets for all of the query's minimizers.
        // Start from the smallest list to minimise work.
        let seed_kmer = kmers.iter()
            .min_by_key(|&&k| inverted.get(&k).map_or(0, |s| s.len()))
            .unwrap();
        let mut candidates: HashSet<String> = inverted
            .get(seed_kmer)
            .cloned()
            .unwrap_or_default();

        for &kmer in kmers {
            if candidates.is_empty() { break; }
            match inverted.get(&kmer) {
                Some(set) => candidates.retain(|h| set.contains(h)),
                None => { candidates.clear(); break; }
            }
        }

        // Keep only unabsorbed candidates that are longer (or equal) and within 10 bp
        candidates.remove(hash);
        candidates.retain(|h| {
            if absorbed.contains(h) { return false; }
            let clen = table[h].0.len();
            clen >= seq_len && clen - seq_len <= MAX_LEN_DIFF
        });

        if candidates.is_empty() { continue; }
        else{
            log::debug!(
                "Merging {} candidates into {}: {} ({} bp) → {} candidates within {} bp",
                candidates.len(), hash, hash, seq_len, candidates.len(), MAX_LEN_DIFF
            );
        }

        // Representative = highest total depth among candidates
        let best = candidates.iter()
            .max_by_key(|h| table[h.as_str()].1.iter().sum::<u64>())
            .unwrap()
            .clone();

        // Merge counts into representative
        let my_counts = table[hash].1.clone();
        for (a, b) in table.get_mut(&best).unwrap().1.iter_mut().zip(&my_counts) {
            *a += b;
        }

        // Propagate lineage to representative if it lacks one
        if !hash_to_lineage.contains_key(&best) {
            if let Some(lineage) = hash_to_lineage.get(hash).cloned() {
                hash_to_lineage.insert(best.clone(), lineage);
            }
        }

        // Remove absorbed hash from inverted index so it can't be a future candidate
        for &kmer in kmers {
            if let Some(set) = inverted.get_mut(&kmer) {
                set.remove(hash);
            }
        }
        absorbed.insert(hash.clone());
    }

    for h in &absorbed {
        table.remove(h);
        hash_to_lineage.remove(h);
    }

    let n_absorbed = absorbed.len();
    if n_absorbed > 0 {
        log::info!(
            "Fuzzy merge: {} → {} unique ASVs ({} near-identical sequences absorbed)",
            n_before, table.len(), n_absorbed,
        );
    }
    n_absorbed
}

// ── public entry point ────────────────────────────────────────────────────────

pub fn export(args: &cli::ExportArgs) {
    let output_dir = Path::new(&args.output_dir);
    std::fs::create_dir_all(output_dir).expect("Failed to create output directory");

    let n_dirs = args.input_dirs.len();

    // ── first pass: determine column structure from feature-table.tsv ─────────
    // Each directory may contribute 1 (single-sample) or N (pooled) columns.
    let mut dir_col_offsets: Vec<usize> = Vec::with_capacity(n_dirs);
    let mut dir_col_counts: Vec<usize> = Vec::with_capacity(n_dirs);
    let mut sample_names: Vec<String> = Vec::new();

    for input_dir in &args.input_dirs {
        let dir = Path::new(input_dir);
        dir_col_offsets.push(sample_names.len());
        match feature_table_from_dir(dir) {
            Some((names, _)) => {
                dir_col_counts.push(names.len());
                sample_names.extend(names);
            }
            None => {
                dir_col_counts.push(1);
                sample_names.push(sample_name_from_dir(dir));
            }
        }
    }

    let total_cols = sample_names.len();
    let mut asv_table: BTreeMap<String, (Vec<u8>, Vec<u64>)> = BTreeMap::new();
    let mut total_reads: Vec<u64> = vec![0u64; total_cols];

    // hash → QIIME lineage string, built from asv_mappings
    let mut hash_to_lineage: HashMap<String, String> = HashMap::new();

    // ── second pass: fill in depths from feature-table.tsv ───────────────────
    for (dir_idx, input_dir) in args.input_dirs.iter().enumerate() {
        let dir = Path::new(input_dir);
        let col_start = dir_col_offsets[dir_idx];
        let n_cols = dir_col_counts[dir_idx];

        // Read per-ASV depths from feature-table.tsv (preferred) or fall back to FASTA header
        let ft_depths: HashMap<String, Vec<u64>> = feature_table_from_dir(dir)
            .map(|(_, d)| d)
            .unwrap_or_default();

        // ASVs — build a per-sample token→hash map for joining with asv_mappings
        let mut token_to_hash: HashMap<String, String> = HashMap::new();
        let fasta_path = dir.join(ASV_FILE);
        match load_fasta_with_needletail(&fasta_path) {
            Ok(seqs) => {
                for (header, seq) in &seqs {
                    let token = header.trim_start_matches('>')
                        .split_whitespace()
                        .next()
                        .unwrap_or("")
                        .to_string();
                    let hash = seq_hash(seq);
                    token_to_hash.insert(token.clone(), hash.clone());

                    // Use feature-table.tsv depths; fall back to FASTA header total depth
                    let per_col_depths: Vec<u64> = ft_depths.get(&token)
                        .cloned()
                        .unwrap_or_else(|| vec![depth_from_header_total(header)]);

                    let entry = asv_table
                        .entry(hash)
                        .or_insert_with(|| (seq.clone(), vec![0u64; total_cols]));
                    for (col_i, &depth) in per_col_depths.iter().enumerate().take(n_cols) {
                        total_reads[col_start + col_i] += depth;
                        entry.1[col_start + col_i] += depth;
                    }
                }
            }
            Err(e) => {
                log::error!("Could not read {}: {}", fasta_path.display(), e);
                continue;
            }
        }

        // asv_mappings.tsv — join header tokens → hashes → taxonomy lineage
        let mappings_path = dir.join("asv_mappings.tsv");
        if mappings_path.exists() {
            match read_asv_mapping_keys(&mappings_path) {
                Ok(pairs) => {
                    for (header_token, lineage) in pairs {
                        if let Some(hash) = token_to_hash.get(&header_token) {
                            hash_to_lineage.entry(hash.clone())
                                .or_insert_with(|| lineage.clone());
                        }
                    }
                }
                Err(e) => log::warn!("Could not read {}: {}", mappings_path.display(), e),
            }
        }
    }

    log::info!("Loaded {} input directories ({} total sample columns), {} unique ASVs",
        n_dirs, total_cols, asv_table.len());

    // ── apply --relabel ───────────────────────────────────────────────────────
    if let Some(ref labels) = args.relabel {
        if labels.len() != total_cols {
            log::error!(
                "--relabel: {} label(s) provided for {} total sample column(s); counts must match",
                labels.len(), total_cols,
            );
            std::process::exit(1);
        }
        sample_names = labels.clone();
        log::info!("Sample names overridden via --relabel: {:?}", sample_names);
    }

    // ── warn on duplicate sample names ───────────────────────────────────────
    {
        let mut seen: HashSet<&str> = HashSet::new();
        let mut dups: Vec<&str> = Vec::new();
        for name in &sample_names {
            if !seen.insert(name.as_str()) {
                dups.push(name.as_str());
            }
        }
        if !dups.is_empty() {
            dups.sort_unstable();
            dups.dedup();
            log::warn!(
                "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            log::warn!(
                 "WARNING: DUPLICATE SAMPLE NAMES DETECTED: [{}]", dups.join(", "));
            log::warn!("Duplicate column names in merged outputs will produce incorrect results in downstream tools.");
            log::warn!("Use --relabel to assign unique names.");
            log::warn!("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        }
    }

    // ── fuzzy merge ───────────────────────────────────────────────────────────
    if !args.no_fuzzy {
        fuzzy_merge_table(&mut asv_table, &mut hash_to_lineage);
    }

    // ── write ASV outputs ─────────────────────────────────────────────────────
    let out = output_dir.display();

    let ft_path = output_dir.join("merged_feature_table.tsv");
    write_feature_table(&asv_table, &sample_names, &ft_path)
        .expect("Failed to write merged_feature_table.tsv");
    log::info!("Wrote {}", ft_path.display());

    let rs_path = output_dir.join("merged_rep_seqs.fasta");
    write_rep_seqs(&asv_table, &rs_path)
        .expect("Failed to write merged_rep_seqs.fasta");
    log::info!("Wrote {}", rs_path.display());

    // ── write taxonomy outputs ────────────────────────────────────────────────
    let asv_tax_path = output_dir.join("merged_asv_taxonomy.tsv");
    write_asv_taxonomy_file(&asv_table, &hash_to_lineage, &asv_tax_path)
        .expect("Failed to write merged_asv_taxonomy.tsv");
    log::info!("Wrote {} ({} ASVs classified)", asv_tax_path.display(), hash_to_lineage.len());

    let taxon_counts_path = output_dir.join("merged_taxon_counts.tsv");
    write_taxon_table(&asv_table, &hash_to_lineage, &sample_names, &taxon_counts_path)
        .expect("Failed to write merged_taxon_counts.tsv");
    log::info!("Wrote {}", taxon_counts_path.display());

    log::info!(
        "To import into QIIME2:\n\
         \n\
         # Feature table\n\
         biom convert -i {out}/merged_feature_table.tsv -o feature-table.biom --table-type='OTU table' --to-hdf5\n\
         qiime tools import --type 'FeatureTable[Frequency]' --input-path feature-table.biom --output-path feature-table.qza\n\
         \n\
         # Representative sequences\n\
         qiime tools import --type 'FeatureData[Sequence]' \\\n\
           --input-path {out}/merged_rep_seqs.fasta --output-path rep-seqs.qza\n\
         \n\
         # If `savont classify / sintax` was run: ASV-level taxonomy (use with feature-table.qza for taxa barplot)\n\
         qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat \\\n\
           --input-path {out}/merged_asv_taxonomy.tsv --output-path taxonomy.qza\n\
         \n\
         # If `savont classify / sintax` was run: Taxonomy bar plot\n\
         qiime taxa barplot --i-table feature-table.qza --i-taxonomy taxonomy.qza \\\n\
           --o-visualization taxa-bar-plots.qzv\\\n",
        out = out,
    );
    log::info!("Export complete.");
}

// ── unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    // ── seq_hash ──────────────────────────────────────────────────────────────

    #[test]
    fn test_seq_hash_deterministic() {
        let seq = b"ACGTACGTACGT";
        assert_eq!(seq_hash(seq), seq_hash(seq));
    }

    #[test]
    fn test_seq_hash_different_sequences() {
        let h1 = seq_hash(b"ACGTACGT");
        let h2 = seq_hash(b"TGCATGCA");
        assert_ne!(h1, h2);
    }

    #[test]
    fn test_seq_hash_case_insensitive() {
        assert_eq!(seq_hash(b"ACGT"), seq_hash(b"acgt"));
    }

    #[test]
    fn test_seq_hash_rc_canonical() {
        let fwd = b"ACGTTGCAACGT";
        let rc = crate::utils::reverse_complement(fwd);
        assert_eq!(seq_hash(fwd), seq_hash(&rc),
            "seq_hash must return the same value for a sequence and its reverse complement");
    }

    // ── depth_from_header ─────────────────────────────────────────────────────

    #[test]
    fn test_depth_from_header_normal() {
        assert_eq!(depth_from_header_total("final_consensus_0_depth_42"), 42);
        assert_eq!(depth_from_header_total(">final_consensus_0_depth_42"), 42);
    }

    #[test]
    fn test_depth_from_header_missing() {
        assert_eq!(depth_from_header_total("some_header_no_number"), 0);
    }

    // ── compute_minimizers ────────────────────────────────────────────────────

    #[test]
    fn test_compute_minimizers_nonempty() {
        let seq: Vec<u8> = b"ACGT".iter().cycle().take(200).cloned().collect();
        let mins = compute_minimizers(&seq);
        assert!(!mins.is_empty(), "200 bp sequence should yield non-empty minimizers");
        // Verify deduplication: no duplicates
        let mut sorted = mins.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(mins.len(), sorted.len(), "minimizers should be deduplicated");
    }

    // ── fuzzy_merge_table ─────────────────────────────────────────────────────

    #[test]
    fn test_fuzzy_merge_absorbs_shorter_prefix() {
        let s1: Vec<u8> = b"ACGT".iter().cycle().take(100).cloned().collect();
        let s2: Vec<u8> = {
            let mut v: Vec<u8> = b"ACGT".iter().cycle().take(100).cloned().collect();
            v.extend_from_slice(b"ACGTACG");
            v
        };
        let h1 = seq_hash(&s1);
        let h2 = seq_hash(&s2);

        let mut table: BTreeMap<String, (Vec<u8>, Vec<u64>)> = BTreeMap::new();
        table.insert(h1.clone(), (s1.clone(), vec![3, 0]));
        table.insert(h2.clone(), (s2.clone(), vec![0, 5]));

        let mut hash_to_lineage: HashMap<String, String> = HashMap::new();
        hash_to_lineage.insert(h1.clone(), "Bacteria;Firmicutes".to_string());

        let n_absorbed = fuzzy_merge_table(&mut table, &mut hash_to_lineage);

        assert_eq!(n_absorbed, 1, "s1 should be absorbed into s2");
        assert!(!table.contains_key(&h1), "h1 should be removed after absorption");
        assert!(table.contains_key(&h2), "h2 should remain as representative");

        let (_, counts) = &table[&h2];
        assert_eq!(counts[0], 3, "counts[0] should be summed");
        assert_eq!(counts[1], 5, "counts[1] should be unchanged");

        // lineage propagates from h1 to h2 when h2 had no lineage
        assert_eq!(
            hash_to_lineage.get(&h2).map(|s| s.as_str()),
            Some("Bacteria;Firmicutes"),
            "lineage should propagate from absorbed h1 to h2"
        );
    }

    #[test]
    fn test_fuzzy_merge_no_merge_different_sequences() {
        // s1: ACGT repeating (offset 0), s2: TGCA repeating (offset 2)
        let s1: Vec<u8> = b"ACGT".iter().cycle().take(100).cloned().collect();
        let s2: Vec<u8> = b"TGCA".iter().cycle().take(100).cloned().collect();
        let h1 = seq_hash(&s1);
        let h2 = seq_hash(&s2);

        let mut table: BTreeMap<String, (Vec<u8>, Vec<u64>)> = BTreeMap::new();
        table.insert(h1.clone(), (s1, vec![10]));
        table.insert(h2.clone(), (s2, vec![10]));

        let mut hash_to_lineage: HashMap<String, String> = HashMap::new();
        let n_absorbed = fuzzy_merge_table(&mut table, &mut hash_to_lineage);

        assert_eq!(n_absorbed, 0, "distinct sequences should not be merged");
        assert_eq!(table.len(), 2, "both sequences should remain");
    }

    #[test]
    fn test_fuzzy_merge_no_merge_beyond_length_diff() {
        // s2 is s1 + 15 extra bytes, exceeding MAX_LEN_DIFF=10
        let s1: Vec<u8> = b"ACGT".iter().cycle().take(100).cloned().collect();
        let s2: Vec<u8> = {
            let mut v: Vec<u8> = b"ACGT".iter().cycle().take(100).cloned().collect();
            v.extend_from_slice(b"ACGTACGTACGTACG"); // 15 extra bytes
            v
        };
        let h1 = seq_hash(&s1);
        let h2 = seq_hash(&s2);

        let mut table: BTreeMap<String, (Vec<u8>, Vec<u64>)> = BTreeMap::new();
        table.insert(h1.clone(), (s1, vec![10]));
        table.insert(h2.clone(), (s2, vec![10]));

        let mut hash_to_lineage: HashMap<String, String> = HashMap::new();
        let n_absorbed = fuzzy_merge_table(&mut table, &mut hash_to_lineage);

        assert_eq!(n_absorbed, 0, "length diff > MAX_LEN_DIFF should prevent merge");
        assert_eq!(table.len(), 2, "both sequences should remain");
    }

    #[test]
    fn test_fuzzy_merge_lineage_not_overwritten() {
        let s1: Vec<u8> = b"ACGT".iter().cycle().take(100).cloned().collect();
        let s2: Vec<u8> = {
            let mut v: Vec<u8> = b"ACGT".iter().cycle().take(100).cloned().collect();
            v.extend_from_slice(b"ACGT");
            v
        };
        let h1 = seq_hash(&s1);
        let h2 = seq_hash(&s2);

        let mut table: BTreeMap<String, (Vec<u8>, Vec<u64>)> = BTreeMap::new();
        table.insert(h1.clone(), (s1, vec![3]));
        table.insert(h2.clone(), (s2, vec![5]));

        let mut hash_to_lineage: HashMap<String, String> = HashMap::new();
        hash_to_lineage.insert(h1.clone(), "Bacteria;Firmicutes".to_string());
        hash_to_lineage.insert(h2.clone(), "Bacteria;Proteobacteria".to_string());

        let n_absorbed = fuzzy_merge_table(&mut table, &mut hash_to_lineage);

        assert_eq!(n_absorbed, 1, "s1 should be absorbed into s2");
        assert_eq!(
            hash_to_lineage.get(&h2).map(|s| s.as_str()),
            Some("Bacteria;Proteobacteria"),
            "h2's original lineage must not be overwritten"
        );
    }

    // ── read_asv_mapping_keys ─────────────────────────────────────────────────

    #[test]
    fn test_read_asv_mapping_keys() {
        let mut tmp = NamedTempFile::new().unwrap();
        writeln!(
            tmp,
            "asv_header\tdepth\talignment_identity\tnumber_mismatches\ttax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom\treference"
        ).unwrap();
        writeln!(
            tmp,
            "final_consensus_0_depth_42\t42\t99.50\t2\t12345\tEscherichia coli\tEscherichia\tEnterobacteriaceae\tEnterobacterales\tGammaproteobacteria\tPseudomonadota\t\tBacteria\tref_seq"
        ).unwrap();
        writeln!(
            tmp,
            "final_consensus_1_depth_20\t20\tNA\tNA\tNA\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED"
        ).unwrap();
        tmp.flush().unwrap();

        let pairs = read_asv_mapping_keys(tmp.path()).unwrap();
        assert_eq!(pairs.len(), 2);

        assert_eq!(pairs[0].0, "final_consensus_0_depth_42");
        assert_eq!(
            pairs[0].1,
            "Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli"
        );

        assert_eq!(pairs[1].0, "final_consensus_1_depth_20");
        assert_eq!(
            pairs[1].1,
            "UNCLASSIFIED;UNCLASSIFIED;UNCLASSIFIED;UNCLASSIFIED;UNCLASSIFIED;UNCLASSIFIED;UNCLASSIFIED"
        );
    }
}
