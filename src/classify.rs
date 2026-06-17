use crate::cli;
use crate::constants::ASV_FILE;
use std::path::Path;
use crate::taxonomy;
use std::sync::Mutex;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};

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

    let fasta_str = db.fasta_path.to_str().unwrap();
    let mmi_path = format!("{}.mmi", fasta_str);
    let aligner = if Path::new(&mmi_path).exists() {
        log::info!("Loading pre-built minimap2 index: {}", mmi_path);
        minimap2::Aligner::builder()
            .map_ont()
            .with_index_threads(args.threads)
            .with_cigar()
            .with_index(&mmi_path, None)
            .expect("Failed to load minimap2 index")
    } else {
        log::info!("Building minimap2 index from {} (saving to {})", fasta_str, mmi_path);
        minimap2::Aligner::builder()
            .map_ont()
            .with_index_threads(args.threads)
            .with_cigar()
            .with_index(fasta_str, Some(&mmi_path))
            .expect("Failed to build minimap2 index")
    };

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
                                    let db_key = (db.extract_key)(db_header);

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

/// Read feature-table.tsv and return (sample_names, depths_per_asv[asv_idx][sample_k]).
/// ASV ordering matched by OTU ID token to FASTA header token.
fn read_feature_table_for_classify(
    ft_path: &Path,
    consensus_sequences: &[(String, Vec<u8>)],
) -> Option<(Vec<String>, Vec<Vec<usize>>)> {
    let f = std::fs::File::open(ft_path).ok()?;
    let mut lines = BufReader::new(f).lines().flatten();

    let header_line = lines.find(|l| l.starts_with("#OTU ID"))?;
    let cols: Vec<&str> = header_line.split('\t').collect();
    let sample_names: Vec<String> = cols[1..].iter().map(|s| s.to_string()).collect();
    let n_samples = sample_names.len();
    if n_samples == 0 { return None; }

    let mut otu_depths: HashMap<String, Vec<usize>> = HashMap::new();
    for line in lines {
        if line.is_empty() || line.starts_with('#') { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() { continue; }
        let otu_id = fields[0].to_string();
        let depths: Vec<usize> = (1..=n_samples)
            .map(|i| fields.get(i).and_then(|v| v.parse().ok()).unwrap_or(0))
            .collect();
        otu_depths.insert(otu_id, depths);
    }

    let per_asv: Vec<Vec<usize>> = consensus_sequences.iter().map(|(header, _)| {
        let token = header.trim_start_matches('>').split_whitespace().next().unwrap_or("");
        otu_depths.get(token).cloned().unwrap_or_else(|| vec![0; n_samples])
    }).collect();

    Some((sample_names, per_asv))
}

/// Write wide-format species abundance table with one column per sample.
fn write_species_abundance_pooled(
    classifications: &[taxonomy::AsvClassification],
    per_asv_per_sample: &[Vec<usize>],
    sample_names: &[String],
    output_path: &Path,
) -> std::io::Result<()> {
    use std::collections::BTreeMap;
    let mut file = std::fs::File::create(output_path)?;

    write!(file, "species\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom")?;
    for name in sample_names { write!(file, "\t{}", name)?; }
    writeln!(file)?;

    let n_samples = sample_names.len();
    let sample_totals: Vec<usize> = (0..n_samples)
        .map(|k| per_asv_per_sample.iter().map(|s| s.get(k).copied().unwrap_or(0)).sum())
        .collect();

    let mut taxon_per_sample: BTreeMap<String, (taxonomy::TaxonomyAssignment, Vec<f64>)> = BTreeMap::new();
    for cls in classifications {
        if let Some(ref tax) = cls.taxonomy {
            let key = format!("{}|{}|{}|{}|{}|{}|{}|{}", tax.species, tax.genus, tax.family, tax.order, tax.class, tax.phylum, tax.clade, tax.superkingdom);
            let asv_idx = cls.asv_id.trim_start_matches("ASV_").parse::<usize>().unwrap_or(0);
            let entry = taxon_per_sample.entry(key).or_insert_with(|| (tax.clone(), vec![0.0; n_samples]));
            for k in 0..n_samples {
                let depth = per_asv_per_sample.get(asv_idx).and_then(|s| s.get(k)).copied().unwrap_or(0);
                if sample_totals[k] > 0 {
                    entry.1[k] += depth as f64 / sample_totals[k] as f64;
                }
            }
        }
    }

    let mut sorted: Vec<_> = taxon_per_sample.into_iter().collect();
    sorted.sort_by(|a, b| {
        let sum_a: f64 = a.1.1.iter().sum();
        let sum_b: f64 = b.1.1.iter().sum();
        sum_b.partial_cmp(&sum_a).unwrap()
    });

    for (_, (tax, abundances)) in sorted {
        write!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", tax.species, tax.genus, tax.family, tax.order, tax.class, tax.phylum, tax.clade, tax.superkingdom)?;
        for &a in &abundances { write!(file, "\t{:.6}", a)?; }
        writeln!(file)?;
    }
    Ok(())
}

/// Write wide-format genus abundance table with one column per sample.
fn write_genus_abundance_pooled(
    classifications: &[taxonomy::AsvClassification],
    per_asv_per_sample: &[Vec<usize>],
    sample_names: &[String],
    output_path: &Path,
) -> std::io::Result<()> {
    use std::collections::BTreeMap;
    let mut file = std::fs::File::create(output_path)?;

    write!(file, "genus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom")?;
    for name in sample_names { write!(file, "\t{}", name)?; }
    writeln!(file)?;

    let n_samples = sample_names.len();
    let sample_totals: Vec<usize> = (0..n_samples)
        .map(|k| per_asv_per_sample.iter().map(|s| s.get(k).copied().unwrap_or(0)).sum())
        .collect();

    let mut taxon_per_sample: BTreeMap<String, (taxonomy::TaxonomyAssignment, Vec<f64>)> = BTreeMap::new();
    for cls in classifications {
        if let Some(ref tax) = cls.taxonomy {
            let key = format!("{}|{}|{}|{}|{}|{}", tax.genus, tax.family, tax.order, tax.class, tax.phylum, tax.clade);
            let asv_idx = cls.asv_id.trim_start_matches("ASV_").parse::<usize>().unwrap_or(0);
            let entry = taxon_per_sample.entry(key).or_insert_with(|| (tax.clone(), vec![0.0; n_samples]));
            for k in 0..n_samples {
                let depth = per_asv_per_sample.get(asv_idx).and_then(|s| s.get(k)).copied().unwrap_or(0);
                if sample_totals[k] > 0 {
                    entry.1[k] += depth as f64 / sample_totals[k] as f64;
                }
            }
        }
    }

    let mut sorted: Vec<_> = taxon_per_sample.into_iter().collect();
    sorted.sort_by(|a, b| {
        let sum_a: f64 = a.1.1.iter().sum();
        let sum_b: f64 = b.1.1.iter().sum();
        sum_b.partial_cmp(&sum_a).unwrap()
    });

    for (_, (tax, abundances)) in sorted {
        write!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}", tax.genus, tax.family, tax.order, tax.class, tax.phylum, tax.clade, tax.superkingdom)?;
        for &a in &abundances { write!(file, "\t{:.6}", a)?; }
        writeln!(file)?;
    }
    Ok(())
}

pub fn classify(args: &cli::ClassifyArgs, db: &taxonomy::Database) {
    // Step 2: Load consensus sequences from clustering output
    let input_fasta = Path::new(&args.input_dir).join(ASV_FILE);
    if !input_fasta.exists() {
        eprintln!("ERROR [savont] Input FASTA not found: {}", input_fasta.display());
        std::process::exit(1);
    }

    log::info!("Loading consensus sequences from {}", input_fasta.display());
    let consensus_sequences_result = taxonomy::load_fasta_with_needletail(&input_fasta);
    let consensus_sequences = match consensus_sequences_result {
        Ok(seqs) => seqs,
        Err(e) => {
            log::warn!("WARN [savont] Failed to load consensus sequences: {}. Either the input FASTA is empty or there was an error during parsing.", e);
            Vec::new()
        }
    };
    
    log::info!("Loaded {} consensus sequences", consensus_sequences.len());

    // Step 3: Build aligner using database FASTA file
    log::info!("Building minimap2 index from database FASTA: {}", db.fasta_path.display());

    // Step 5: Collect all best mappings for each ASV
    // Read depths from feature-table.tsv if available (authoritative for pooled samples);
    // fall back to parsing FASTA headers.
    let ft_path = Path::new(&args.input_dir).join("feature-table.tsv");
    let (ft_sample_names, per_asv_per_sample) = read_feature_table_for_classify(&ft_path, &consensus_sequences)
        .unwrap_or_else(|| {
            let depths = taxonomy::extract_depths_from_headers(&consensus_sequences);
            (vec!["sample".to_string()], depths.iter().map(|&d| vec![d]).collect())
        });

    let asv_depths: Vec<usize> = per_asv_per_sample.iter().map(|s| s.iter().sum()).collect();
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
                    args.detailed_unclassified,
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
                args.detailed_unclassified,
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
    if ft_sample_names.len() > 1 {
        write_species_abundance_pooled(&classifications, &per_asv_per_sample, &ft_sample_names, &species_file)
            .expect("Failed to write species abundance file");
    } else {
        taxonomy::write_species_abundance(&classifications, &species_file)
            .expect("Failed to write species abundance file");
    }
    log::info!("Wrote species abundance table to {}", species_file.display());

    let genus_file = output_dir_path.join("genus_abundance.tsv");
    if ft_sample_names.len() > 1 {
        write_genus_abundance_pooled(&classifications, &per_asv_per_sample, &ft_sample_names, &genus_file)
            .expect("Failed to write genus abundance file");
    } else {
        taxonomy::write_genus_abundance(&classifications, &genus_file)
            .expect("Failed to write genus abundance file");
    }
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