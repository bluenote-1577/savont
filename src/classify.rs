use crate::cli;
use crate::constants::ASV_FILE;
use std::path::Path;
use crate::taxonomy;
use std::sync::Mutex;
use rayon::prelude::*;
use std::collections::HashMap;

/// Represents a mapping from an ASV to a database sequence
#[derive(Debug, Clone)]
struct AsvMapping {
    asv_idx: usize,
    _tax_idx: String,  // Index into the unique tax_id list
    hit_reference_id: String,
    index: usize, 
    identity: f64,
    nm: u32,
    depth: usize,
    species: String,
}

/// Run EM algorithm to estimate taxonomic abundances given multi-mapped ASVs
fn run_em_algorithm(
    mappings: &[AsvMapping],
    num_taxa: usize,
    total_reads: usize,
    convergence_threshold: f64,
) -> Vec<f64> {
    // Initialize abundance estimates uniformly
    let mut tax_abundances = vec![1.0 / num_taxa as f64; num_taxa];
    let mut iteration = 0;
    const MAX_ITERATIONS: usize = 1000;

    let mut id_to_species_map = HashMap::new();
    for mapping in mappings {
        id_to_species_map.insert(mapping.index, mapping.species.clone());
    }

    log::info!("Starting EM algorithm with {} taxa and {} mappings", num_taxa, mappings.len());

    loop {
        iteration += 1;
        let mut new_tax_abundances = vec![0.0; num_taxa];

        // E-step + M-step combined:
        // For each ASV, distribute its reads proportionally across its mapped taxa
        // based on current abundance estimates
        let mut asv_mappings_grouped: HashMap<usize, Vec<&AsvMapping>> = HashMap::new();
        for mapping in mappings {
            asv_mappings_grouped
                .entry(mapping.asv_idx)
                .or_insert_with(Vec::new)
                .push(mapping);
        }

        for (_asv_idx, asv_maps) in &asv_mappings_grouped {
            // Calculate denominator: sum of (depth * current_abundance) for all taxa this ASV maps to
            let denominator: f64 = asv_maps.iter()
                .map(|m| tax_abundances[m.index])
                .sum();

            if denominator > 0.0 {
                // Distribute ASV's depth to each taxon proportionally
                for mapping in asv_maps.iter() {
                    let contribution = (mapping.depth as f64) * tax_abundances[mapping.index] / denominator;
                    new_tax_abundances[mapping.index] += contribution;
                }
            }
        }

        // Normalize to get new abundance estimates
        let total_assigned: f64 = new_tax_abundances.iter().sum();
        if total_assigned > 0.0 {
            for abundance in new_tax_abundances.iter_mut() {
                *abundance /= total_reads as f64;
            }
        }

        // Check convergence
        let max_change = tax_abundances.iter()
            .zip(new_tax_abundances.iter())
            .map(|(old, new)| (old - new).abs())
            .fold(0.0, f64::max);

        log::debug!("EM iteration {}: max abundance change = {:.6e}", iteration, max_change);

        // Write the 5 top most changed taxa for debugging
        let mut changes: Vec<(usize, f64)> = tax_abundances.iter()
            .zip(new_tax_abundances.iter())
            .enumerate()
            .map(|(idx, (old, new))| (idx, (old - new).abs()))
            .collect();
        changes.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        log::debug!("Top 5 abundance changes:");
        for (idx, change) in changes.iter().take(5) {
            log::debug!("  Taxon {}, Species {} : change = {:.6e}", idx, id_to_species_map.get(&idx).unwrap(), change);
        }

        tax_abundances = new_tax_abundances;

        if max_change < convergence_threshold || iteration >= MAX_ITERATIONS {
            log::info!("EM converged after {} iterations (max change: {:.6e})", iteration, max_change);
            break;
        }
    }

    // Filter out very low abundance taxa
    let min_abundance = convergence_threshold;
    for abundance in tax_abundances.iter_mut() {
        if *abundance < min_abundance {
            *abundance = 0.0;
        }
    }

    tax_abundances
}

/// Collect all best mappings (those with minimum NM) for each ASV
fn collect_best_mappings(
    consensus_sequences: &[(String, Vec<u8>)],
    asv_depths: &[usize],
    db: &taxonomy::Database,
    args: &cli::ClassifyArgs,
) -> Vec<(usize, String, f64, u32, usize, String, String)> {

    let aligner = minimap2::Aligner::builder()
        .map_ont()
        .with_index_threads(args.threads)
        .with_cigar()
        .with_index(db.fasta_path.to_str().unwrap(), None)
        .expect("Failed to create aligner with database index");

    log::info!("Aligning {} consensus sequences to database", consensus_sequences.len());


    let all_mappings = Mutex::new(Vec::new());

    consensus_sequences.par_iter().enumerate().for_each(|(asv_idx, (header, sequence))| {
        let asv_header = header.trim_start_matches('>').to_string();

        // Align to database
        let alignment_result = aligner.map(sequence, true, false, None, None, None);

        if let Ok(mappings) = alignment_result {
            if !mappings.is_empty() {
                // Find the minimum NM value (best alignment quality)
                let min_nm = mappings.first().and_then(|m| m.alignment.as_ref().map(|a| a.nm));

                if let Some(min_nm) = min_nm {
                    // Collect ALL mappings with the minimum NM value
                    for mapping in mappings.iter() {
                        if let Some(alignment) = &mapping.alignment {
                            if alignment.nm == min_nm {
                                let alignment_length = mapping.query_end - mapping.query_start;
                                let identity = 100.0 * (1.0 - (alignment.nm as f64 / alignment_length as f64));

                                if let Some(db_header) = &mapping.target_name {
                                    // Extract ID based on database type
                                    let db_key = if args.db_type.silva_db.is_some() {
                                        taxonomy::extract_silva_accession_from_header(db_header)
                                    } else {
                                        taxonomy::extract_tax_id_from_header(db_header)
                                    };

                                    if let Some(key) = db_key {
                                        if db.taxonomy.contains_key(&key) {
                                            let mapping_record = (asv_idx, key, identity, alignment.nm as u32, asv_depths[asv_idx], 
                                                asv_header.clone(), (*(*mapping.target_name.as_ref().unwrap())).clone());
                                            all_mappings.lock().unwrap().push(mapping_record);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    });

    all_mappings.into_inner().unwrap()
}

pub fn classify(args: &cli::ClassifyArgs, db: &taxonomy::Database) {
    // Step 2: Load consensus sequences from clustering output
    let input_fasta = Path::new(&args.input_dir).join(ASV_FILE);
    if !input_fasta.exists() {
        eprintln!("ERROR [savont] Input FASTA not found: {}", input_fasta.display());
        std::process::exit(1);
    }

    log::info!("Loading consensus sequences from {}", input_fasta.display());
    let consensus_sequences = taxonomy::load_fasta_with_needletail(&input_fasta)
        .expect("Failed to load consensus sequences");
    log::info!("Loaded {} consensus sequences", consensus_sequences.len());

    // Step 3: Build aligner using database FASTA file
    log::info!("Building minimap2 index from database FASTA: {}", db.fasta_path.display());
    // Step 5: Collect all best mappings for each ASV
    let asv_depths = taxonomy::extract_depths_from_headers(&consensus_sequences);
    let total_reads: usize = asv_depths.iter().sum();

    let all_mappings = collect_best_mappings(&consensus_sequences, &asv_depths, &db, &args);
    log::info!("Collected {} total mappings from {} ASVs", all_mappings.len(), consensus_sequences.len());

    // Step 6: Build tax_id index and mapping matrix
    let mut tax_id_to_idx: HashMap<String, usize> = HashMap::new();
    let mut idx_to_tax_id: Vec<String> = Vec::new();

    for (_, tax_id, _, _, _, _, _) in &all_mappings {
        if !tax_id_to_idx.contains_key(tax_id) {
            let idx = idx_to_tax_id.len();
            tax_id_to_idx.insert(tax_id.clone(), idx);
            idx_to_tax_id.push(tax_id.clone());
        }
    }

    log::info!("Found {} unique taxonomic IDs", idx_to_tax_id.len());

    // Convert mappings to indexed format
    let mappings: Vec<AsvMapping> = all_mappings.iter()
        .map(|(asv_idx, tax_id, identity, nm, depth, _, hit_reference_id)| AsvMapping {
            asv_idx: *asv_idx,
            _tax_idx: tax_id.clone(),
            index: *tax_id_to_idx.get(tax_id).unwrap(),
            hit_reference_id: hit_reference_id.clone(),
            identity: *identity,
            nm: *nm,
            depth: *depth,
            species: db.taxonomy.get(tax_id).unwrap().species.clone(),
        })
        .collect();

    // Step 7: Run EM algorithm to distribute abundances
    log::info!("Running EM algorithm to distribute abundances");
    let convergence_threshold = 0.1 / total_reads as f64;
    let tax_abundances = run_em_algorithm(&mappings, idx_to_tax_id.len(), total_reads, convergence_threshold);

    // Step 8: Build classifications from EM results
    let mut classifications: Vec<taxonomy::AsvClassification> = Vec::new();
    let mut secondary_classifications: Vec<taxonomy::AsvClassification> = Vec::new();

    for asv_idx in 0..consensus_sequences.len() {
        let (header, _) = &consensus_sequences[asv_idx];
        let asv_id = format!("ASV_{}", asv_idx);
        let asv_header = header.trim_start_matches('>').split_whitespace().next().unwrap_or("").to_string();

        // Find all mappings for this ASV
        let asv_mappings: Vec<&AsvMapping> = mappings.iter()
            .filter(|m| m.asv_idx == asv_idx)
            .collect();

        if !asv_mappings.is_empty() {
            // Log all top hits for this ASV
            log::debug!("ASV {} ({} depth, {} total hits):", asv_id, asv_depths[asv_idx], asv_mappings.len());

            // Sort mappings by EM abundance (descending) to show best hits first
            let mut sorted_mappings = asv_mappings.clone();
            sorted_mappings.sort_by(|a, b| {
                tax_abundances[b.index]
                    .partial_cmp(&tax_abundances[a.index])
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            for (rank, mapping) in sorted_mappings.iter().enumerate() {
                let tax_id = &idx_to_tax_id[mapping.index];
                log::debug!(
                    "  Hit #{}: tax_id={}, species={}, identity={:.2}%, nm={}, EM_abundance={:.6}",
                    rank + 1, tax_id, mapping.species, mapping.identity, mapping.nm, tax_abundances[mapping.index]
                );

                let taxonomy_entry = db.taxonomy.get(tax_id).unwrap();
                let taxonomy_assignment = taxonomy::TaxonomyAssignment::from_taxonomy_entry(
                    taxonomy_entry,
                    mapping.identity,
                    args.species_threshold,
                    args.genus_threshold,
                    &asv_header,
                );
                secondary_classifications.push(taxonomy::AsvClassification {
                    asv_id: asv_id.clone(),
                    asv_header: asv_header.clone(),
                    abundance: asv_depths[asv_idx] as f64 / total_reads as f64,
                    best_hit_tax_id: Some(tax_id.clone()),
                    identity: Some(mapping.identity),
                    taxonomy: Some(taxonomy_assignment),
                    nm: Some(mapping.nm as usize),
                    hit_reference_id: mapping.hit_reference_id.clone(),
                });
            }

            // Take the mapping with highest EM-estimated abundance (or first if tied)
            let best_mapping = asv_mappings.iter()
                .max_by(|a, b| {
                    tax_abundances[a.index]
                        .partial_cmp(&tax_abundances[b.index])
                        .unwrap()
                })
                .unwrap();

            let tax_id = &idx_to_tax_id[best_mapping.index];
            let taxonomy_entry = db.taxonomy.get(tax_id).unwrap();
            let taxonomy_assignment = taxonomy::TaxonomyAssignment::from_taxonomy_entry(
                taxonomy_entry,
                best_mapping.identity,
                args.species_threshold,
                args.genus_threshold,
                &asv_header,
            );

            let asv_depth = asv_depths[asv_idx];
            let abundance = asv_depth as f64 / total_reads as f64;

            classifications.push(taxonomy::AsvClassification {
                asv_id: asv_id.clone(),
                asv_header,
                abundance,
                best_hit_tax_id: Some(tax_id.clone()),
                identity: Some(best_mapping.identity),
                taxonomy: Some(taxonomy_assignment),
                nm: Some(best_mapping.nm as usize),
                hit_reference_id: best_mapping.hit_reference_id.clone(),
            });

            log::debug!("Classified ASV {}: tax_id_acc={}, tax_id_tax={}, species={}, genus ={}, abundance={:.6}",
                asv_id, tax_id, taxonomy_entry.tax_id, taxonomy_entry.species, taxonomy_entry.genus, abundance);
        } else {
            // No alignment found
            classifications.push(taxonomy::AsvClassification {
                asv_id,
                asv_header,
                abundance: asv_depths[asv_idx] as f64 / total_reads as f64,
                hit_reference_id: String::new(),
                nm: None,
                best_hit_tax_id: None,
                identity: None,
                taxonomy: None,
            });
        }
    }

    classifications.sort_by(|a, b| b.abundance.partial_cmp(&a.abundance).unwrap());

    // Step 7: Write output files
    //let output_dir_path = Path::new(args.output_dir.as_ref().unwrap_or(args.input_dir.clone()));
    let output_dir_path = if let Some(ref out_dir) = args.output_dir {
        Path::new(out_dir)
    } else {
        Path::new(&args.input_dir)
    };

    let species_file = output_dir_path.join("species_abundance.tsv");
    taxonomy::write_species_abundance(&classifications, &species_file)
        .expect("Failed to write species abundance file");
    log::info!("Wrote species abundance table to {}", species_file.display());

    let genus_file = output_dir_path.join("genus_abundance.tsv");
    taxonomy::write_genus_abundance(&classifications, &genus_file)
        .expect("Failed to write genus abundance file");
    log::info!("Wrote genus abundance table to {}", genus_file.display());

    let mappings_file = output_dir_path.join("asv_mappings.tsv");
    taxonomy::write_asv_mappings(&secondary_classifications, &mappings_file)
        .expect("Failed to write ASV mappings file");
    log::info!("Wrote ASV mappings table to {}", mappings_file.display());

    log::info!("Classification complete! Classified {}/{} ASVs",
        classifications.iter().filter(|c| c.taxonomy.is_some()).count(),
        classifications.len());

    // Classified X species at species and genus level
    log::info!("Classified {}/{} ASVs at species level",
        classifications.iter()
            .filter(|c| {
                if let Some(tax) = &c.taxonomy {
                    !tax.species.is_empty() && !tax.species.contains("UNCLASSIFIED")
                } else {
                    false
                }

            })
            .count(),
        classifications.iter().filter(|c| c.taxonomy.is_some()).count(),
        );

    log::info!("Classified {}/{} ASVs at genus level",
        classifications.iter()
            .filter(|c| {
                if let Some(tax) = &c.taxonomy {
                    !tax.genus.is_empty() && !tax.genus.contains("UNCLASSIFIED")
                } else {
                    false
                }
            })
            .count(), 
            classifications.iter().filter(|c| c.taxonomy.is_some()).count(),
        );


}