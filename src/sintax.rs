use crate::cli;
use crate::constants::ASV_FILE;
use crate::taxonomy;
use crate::types::BYTE_TO_SEQ;
use fxhash::FxHashMap;
use rayon::iter::IntoParallelRefIterator;
use std::collections::HashMap;
use std::collections::HashSet;
use std::path::Path;
use std::sync::Mutex;
use rayon::prelude::*;

const K: usize = 12;
const SUBSAMPLE: usize = 32;

// ── tiny xorshift RNG ────────────────────────────────────────────────────────

struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self { Rng(seed.max(1)) }

    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn next_usize(&mut self, n: usize) -> usize {
        (self.next() as usize) % n
    }
}

// ── k-mer encoding ───────────────────────────────────────────────────────────

fn extract_kmers(seq: &[u8]) -> Vec<u32> {
    let mask = ((1u32 << (2 * K)) - 1) as u32; // 0xFFFF for K=8
    let rc_shift = (2 * (K - 1)) as u32;         // 14 for K=8

    let mut kmers = Vec::with_capacity(seq.len().saturating_sub(K - 1));
    let mut kmer_f: u32 = 0;
    let mut kmer_r: u32 = 0;

    for (i, &b) in seq.iter().enumerate() {
        let enc = BYTE_TO_SEQ[b as usize] as u32;
        let comp = 3 - enc;
        kmer_f = ((kmer_f << 2) | enc) & mask;
        kmer_r = (kmer_r >> 2) | (comp << rc_shift);
        if i + 1 >= K {
            kmers.push(kmer_f.min(kmer_r));
        }
    }
    kmers
}

// ── SINTAX hit ───────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct SintaxHit {
    pub asv_header: String,
    pub depth: usize,
    pub abundance: f64,
    pub species: String,
    pub species_boot: f64,
    pub genus: String,
    pub genus_boot: f64,
    pub family: String,
    pub family_boot: f64,
    pub order: String,
    pub order_boot: f64,
    pub class: String,
    pub class_boot: f64,
    pub phylum: String,
    pub phylum_boot: f64,
    pub superkingdom: String,
    pub superkingdom_boot: f64,
}

// ── output helpers ───────────────────────────────────────────────────────────

fn hit_to_classification(
    hit: &SintaxHit,
    min_bootstrap: f64,
    detailed_unclassified: bool,
) -> taxonomy::AsvClassification {
    let unclassified = if detailed_unclassified {
        format!("UNCLASSIFIED-({})", hit.asv_header)
    } else {
        "UNCLASSIFIED".to_string()
    };

    let apply = |boot: f64, name: &str| -> String {
        if boot >= min_bootstrap { name.to_string() } else { unclassified.clone() }
    };

    let taxonomy = Some(taxonomy::TaxonomyAssignment {
        tax_id: String::new(),
        species:      unclassified.clone(), // sintax is genus-level only
        genus:        apply(hit.genus_boot,        &hit.genus),
        family:       apply(hit.family_boot,       &hit.family),
        order:        apply(hit.order_boot,        &hit.order),
        class:        apply(hit.class_boot,        &hit.class),
        phylum:       apply(hit.phylum_boot,       &hit.phylum),
        clade: String::new(),
        superkingdom: apply(hit.superkingdom_boot, &hit.superkingdom),
        subspecies: String::new(),
        species_subgroup: String::new(),
        species_group: String::new(),
    });

    taxonomy::AsvClassification {
        asv_id: hit.asv_header.clone(),
        asv_header: hit.asv_header.clone(),
        hit_reference_id: String::new(),
        abundance: hit.abundance,
        best_hit_tax_id: None,
        identity: None,
        nm: None,
        taxonomy,
    }
}

fn write_sintax_asv_mappings(
    hits: &[Option<SintaxHit>],
    min_bootstrap: f64,
    path: &Path,
) -> std::io::Result<()> {
    use std::io::Write;
    let mut file = std::fs::File::create(path)?;
    writeln!(
        file,
        "asv_header\tdepth\t\
         species_bootstrap\tgenus_bootstrap\tfamily_bootstrap\t\
         order_bootstrap\tclass_bootstrap\tphylum_bootstrap\tsuperkingdom_bootstrap\t\
         species\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom"
    )?;
    let unc = "UNCLASSIFIED";
    let apply = |boot: f64, name: &str| -> String {
        if boot >= min_bootstrap { name.to_string() } else { unc.to_string() }
    };
    for hit in hits.iter().flatten() {
        writeln!(
            file,
            "{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t\
             {}\t{}\t{}\t{}\t{}\t{}\t{}",
            hit.asv_header,
            hit.depth,
            hit.species_boot,
            hit.genus_boot,
            hit.family_boot,
            hit.order_boot,
            hit.class_boot,
            hit.phylum_boot,
            hit.superkingdom_boot,
            unc,                                      // species: sintax is genus-level max
            apply(hit.genus_boot,        &hit.genus),
            apply(hit.family_boot,       &hit.family),
            apply(hit.order_boot,        &hit.order),
            apply(hit.class_boot,        &hit.class),
            apply(hit.phylum_boot,       &hit.phylum),
            apply(hit.superkingdom_boot, &hit.superkingdom),
        )?;
    }
    Ok(())
}

// ── public entry point ───────────────────────────────────────────────────────

pub fn sintax(args: &cli::SintaxArgs, db: &taxonomy::Database) {
    let input_fasta = Path::new(&args.input_dir).join(ASV_FILE);

    if !input_fasta.exists() {
        log::error!("Input FASTA not found: {}", input_fasta.display());
        std::process::exit(1);
    }

    let sequences = match taxonomy::load_fasta_with_needletail(&input_fasta) {
        Ok(s) if !s.is_empty() => s,
        Ok(_) => { log::warn!("No sequences in {}", input_fasta.display()); return; }
        Err(e) => { log::error!("Failed to load ASV sequences: {}", e); std::process::exit(1); }
    };

    let n_asvs = sequences.len();
    let n_iter = args.n_iter;
    let n_pairs = n_asvs * n_iter;

    let asv_depths = taxonomy::extract_depths_from_headers(&sequences);
    let total_reads: usize = asv_depths.iter().sum();

    // ── Phase 1: Build query map from ASV subsamples ──────────────────────────
    // kmer_query[kmer] = list of (asv_idx, iter_idx) pairs that sampled this kmer
    log::info!("Building SINTAX query map from {} ASVs ({} iterations × {} k-mers)",
        n_asvs, n_iter, SUBSAMPLE);

    let mut kmer_query: FxHashMap<u32, Vec<(u32, u32)>> = FxHashMap::default();

    // Store per-ASV kmer lists temporarily so each iteration can resample
    let asv_kmers: Vec<Vec<u32>> = sequences.iter()
        .map(|(_, seq)| extract_kmers(seq))
        .collect();

    for (asv_i, kmers) in asv_kmers.iter().enumerate() {
        if kmers.is_empty() { continue; }
        for iter_j in 0..n_iter {
            let seed = (asv_i as u64) * (n_iter as u64) + (iter_j as u64) + 1;
            let mut rng = Rng::new(seed);
            for _ in 0..SUBSAMPLE {
                let ki = rng.next_usize(kmers.len());
                let kmer = kmers[ki];
                kmer_query.entry(kmer).or_default().push((asv_i as u32, iter_j as u32));
            }
        }
    }
    drop(asv_kmers);

    log::info!("Query map has {} distinct k-mers", kmer_query.len());

    // ── Phase 2: Stream database, update best scores per (asv, iter) ──────────
    log::info!("Streaming database {}", db.fasta_path.display());

    let best_scores:   Vec<Mutex<u16>> = (0..n_pairs).map(|_| Mutex::new(0u16)).collect();
    let best_taxonomy: Vec<Mutex<Option<taxonomy::TaxonomyEntry>>> = (0..n_pairs).map(|_| Mutex::new(None)).collect();
    let counter = Mutex::new(0usize);

    let db_seqs = taxonomy::load_fasta_with_needletail(&db.fasta_path)
        .expect("Failed to load database FASTA for SINTAX");

    // Parallel
    db_seqs.par_iter()
        .for_each(|(header, seq)| {
        let mut ref_hit_counts: Vec<u16> = vec![0u16; n_pairs];
        let mut touched: Vec<usize> = Vec::new();

        let header = header.trim_start_matches('>');
        let key = match (db.extract_key)(header) {
            Some(k) => k,
            None => return,
        };
        let entry = match db.taxonomy.get(&key) {
            Some(e) => e,
            None => return,
        };

        // Deduplicate reference k-mers so each is counted at most once per reference
        let ref_kmer_set: HashSet<u32> = extract_kmers(seq).into_iter().collect();

        for kmer in ref_kmer_set {
            if let Some(pairs) = kmer_query.get(&kmer) {
                for &(asv_i, iter_j) in pairs {
                    let idx = asv_i as usize * n_iter + iter_j as usize;
                    if ref_hit_counts[idx] == 0 { touched.push(idx); }
                    ref_hit_counts[idx] += 1;
                }
            }
        }

        // Update best scores
        for &idx in &touched {
            let mut best_score = best_scores[idx].lock().unwrap();
            let mut best_entry = best_taxonomy[idx].lock().unwrap();
            if ref_hit_counts[idx] > *best_score { 
                *best_score = ref_hit_counts[idx];
                *best_entry = Some(entry.clone());
            }
        }

        let mut count = counter.lock().unwrap();
        *count += 1;
        if *count % 1000 == 0 {
            log::info!("Processed {} / {} reference sequences...", *count, db_seqs.len());
        }
    });

    log::info!("Finished streaming database...");
    // ── Phase 3: Aggregate votes per ASV ─────────────────────────────────────
    let mut all_hits: Vec<Option<SintaxHit>> = Vec::with_capacity(n_asvs);
    let best_scores = best_scores.into_iter().map(|m| m.into_inner().unwrap()).collect::<Vec<u16>>();
    let best_taxonomy = best_taxonomy.into_iter().map(|m| m.into_inner().unwrap()).collect::<Vec<Option<taxonomy::TaxonomyEntry>>>();

    for asv_i in 0..n_asvs {
        let base = asv_i * n_iter;

        let mut species_votes:     HashMap<String, usize> = HashMap::new();
        let mut genus_votes:       HashMap<String, usize> = HashMap::new();
        let mut family_votes:      HashMap<String, usize> = HashMap::new();
        let mut order_votes:       HashMap<String, usize> = HashMap::new();
        let mut class_votes:       HashMap<String, usize> = HashMap::new();
        let mut phylum_votes:      HashMap<String, usize> = HashMap::new();
        let mut superkingdom_votes: HashMap<String, usize> = HashMap::new();
        let mut classified = 0usize;

        for iter_j in 0..n_iter {
            if let Some(ref e) = best_taxonomy[base + iter_j] {
                if best_scores[base + iter_j] > 0 {
                    classified += 1;
                    *species_votes.entry(e.species.clone()).or_insert(0) += 1;
                    *genus_votes.entry(e.genus.clone()).or_insert(0) += 1;
                    *family_votes.entry(e.family.clone()).or_insert(0) += 1;
                    *order_votes.entry(e.order.clone()).or_insert(0) += 1;
                    *class_votes.entry(e.class.clone()).or_insert(0) += 1;
                    *phylum_votes.entry(e.phylum.clone()).or_insert(0) += 1;
                    *superkingdom_votes.entry(e.superkingdom.clone()).or_insert(0) += 1;
                }
            }
        }

        if classified == 0 {
            all_hits.push(None);
            continue;
        }

        // Debug output for first 10 ASVs

            if asv_i < 10 {
                log::debug!("ASV {}: classified in {}/{} iterations", asv_i, classified, n_iter);
                log::debug!("  Species votes: {:?}", species_votes);
                log::debug!("  Genus votes:   {:?}", genus_votes);
                log::debug!("  Family votes:  {:?}", family_votes);
                log::debug!("  Order votes:   {:?}", order_votes);
                log::debug!("  Class votes:   {:?}", class_votes);
                log::debug!("  Phylum votes:  {:?}", phylum_votes);
                log::debug!("  Superkingdom votes: {:?}", superkingdom_votes);
            }

        let top = |votes: &HashMap<String, usize>| -> (String, f64) {
            votes.iter()
                .max_by_key(|(_, &c)| c)
                .map(|(name, &count)| (name.clone(), count as f64 / n_iter as f64))
                .unwrap_or_default()
        };

        let asv_header = sequences[asv_i].0
            .trim_start_matches('>')
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();

        let (species, species_boot)           = top(&species_votes);
        let (genus, genus_boot)               = top(&genus_votes);
        let (family, family_boot)             = top(&family_votes);
        let (order, order_boot)               = top(&order_votes);
        let (class, class_boot)               = top(&class_votes);
        let (phylum, phylum_boot)             = top(&phylum_votes);
        let (superkingdom, superkingdom_boot) = top(&superkingdom_votes);

        let depth = asv_depths[asv_i];
        let abundance = if total_reads > 0 { depth as f64 / total_reads as f64 } else { 0.0 };

        all_hits.push(Some(SintaxHit {
            asv_header,
            depth,
            abundance,
            species, species_boot,
            genus, genus_boot,
            family, family_boot,
            order, order_boot,
            class, class_boot,
            phylum, phylum_boot,
            superkingdom, superkingdom_boot,
        }));
    }

    // Sort by abundance descending
    all_hits.sort_by(|a, b| {
        let ab = a.as_ref().map_or(0.0, |h| h.abundance);
        let bb = b.as_ref().map_or(0.0, |h| h.abundance);
        bb.partial_cmp(&ab).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Build AsvClassification for reuse of existing abundance writers
    let classifications: Vec<taxonomy::AsvClassification> = all_hits.iter().enumerate()
        .map(|(i, hit)| {
            if let Some(h) = hit {
                hit_to_classification(h, args.min_bootstrap, args.detailed_unclassified)
            } else {
                let header = sequences[i].0
                    .trim_start_matches('>')
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
                taxonomy::AsvClassification {
                    asv_id: header.clone(),
                    asv_header: header,
                    hit_reference_id: String::new(),
                    abundance: asv_depths[i] as f64 / total_reads.max(1) as f64,
                    best_hit_tax_id: None,
                    identity: None,
                    nm: None,
                    taxonomy: None,
                }
            }
        })
        .collect();

    let output_dir = args.output_dir.as_deref()
        .map(Path::new)
        .unwrap_or_else(|| Path::new(&args.input_dir));

    taxonomy::write_genus_abundance(&classifications, &output_dir.join("genus_abundance.tsv"))
        .expect("Failed to write genus_abundance.tsv");
    log::info!("Wrote genus_abundance.tsv");

    write_sintax_asv_mappings(&all_hits, args.min_bootstrap, &output_dir.join("asv_mappings.tsv"))
        .expect("Failed to write asv_mappings.tsv");
    log::info!("Wrote asv_mappings.tsv");

    let classified = all_hits.iter().filter(|h| h.is_some()).count();
    log::info!("SINTAX complete: {}/{} ASVs classified", classified, n_asvs);
}
