use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

/// Represents a taxonomic entry from the database
#[derive(Debug, Clone)]
pub struct TaxonomyEntry {
    pub tax_id: String,
    pub species: String,
    pub genus: String,
    pub family: String,
    pub order: String,
    pub class: String,
    pub phylum: String,
    pub clade: String,
    pub superkingdom: String,
    pub subspecies: String,
    pub species_subgroup: String,
    pub species_group: String,
}

/// Represents an EMU database with sequences and taxonomy
pub struct Database {
    pub fasta_path: PathBuf,
    pub taxonomy: HashMap<String, TaxonomyEntry>,
}

impl Database {
    /// Load a Database from a directory
    pub fn load_emu(db_dir: &Path) -> Result<Self, std::io::Error> {
        let fasta_path = db_dir.join("species_taxid.fasta");
        let taxonomy_path = db_dir.join("taxonomy.tsv");

        if !fasta_path.exists() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("FASTA file not found: {}", fasta_path.display()),
            ));
        }

        if !taxonomy_path.exists() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("Taxonomy file not found: {}", taxonomy_path.display()),
            ));
        }

        log::info!("Loading taxonomy from {}", taxonomy_path.display());
        let taxonomy = Self::load_taxonomy(&taxonomy_path)?;
        log::info!("Loaded {} taxonomy entries", taxonomy.len());

        Ok(Database {
            fasta_path,
            taxonomy,
        })
    }

    /// Load EMU taxonomy from a TSV file
    fn load_taxonomy(path: &Path) -> Result<HashMap<String, TaxonomyEntry>, std::io::Error> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut taxonomy = HashMap::new();

        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;

            // Skip header line
            if line_num == 0 {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 12 {
                log::warn!("Skipping malformed line {}: insufficient fields", line_num + 1);
                continue;
            }

            let entry = TaxonomyEntry {
                tax_id: fields[0].to_string(),
                species: fields[1].to_string(),
                genus: fields[2].to_string(),
                family: fields[3].to_string(),
                order: fields[4].to_string(),
                class: fields[5].to_string(),
                phylum: fields[6].to_string(),
                clade: fields[7].to_string(),
                superkingdom: fields[8].to_string(),
                subspecies: fields[9].to_string(),
                species_subgroup: fields[10].to_string(),
                species_group: fields[11].to_string(),
            };

            taxonomy.insert(entry.tax_id.clone(), entry);
        }

        Ok(taxonomy)
    }

    /// Load Silva database from a directory
    pub fn load_silva(db_dir: &Path) -> Result<Self, std::io::Error> {
        // Find FASTA file (could be .fasta or .fasta.gz)
        let fasta_path = std::fs::read_dir(db_dir)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .find(|p| {
                p.extension().and_then(|e| e.to_str()) == Some("fasta") ||
                p.file_name().and_then(|n| n.to_str()).map_or(false, |n| n.ends_with(".fasta.gz")) || 
                p.file_name().and_then(|n| n.to_str()).map_or(false, |n| n.ends_with(".fa.gz"))
            })
            .ok_or_else(|| std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("No FASTA file found in {}", db_dir.display()),
            ))?;

        // Find taxonomy TSV file (taxmap_*.txt or taxmap_*.txt.gz)
        let taxonomy_path = std::fs::read_dir(db_dir)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .find(|p| {
                p.file_name()
                    .and_then(|n| n.to_str())
                    .map_or(false, |n| n.starts_with("taxmap_") && (n.ends_with(".txt") || n.ends_with(".txt.gz")))
            })
            .ok_or_else(|| std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("No taxmap file found in {}", db_dir.display()),
            ))?;

        log::info!("Loading Silva taxonomy from {}", taxonomy_path.display());
        let taxonomy = Self::load_silva_taxonomy(&taxonomy_path)?;
        log::info!("Loaded {} taxonomy entries", taxonomy.len());

        Ok(Database {
            fasta_path,
            taxonomy,
        })
    }

    /// Load Silva taxonomy from a TSV file
    /// Format: primaryAccession  start  stop  path  organism_name  taxid
    fn load_silva_taxonomy(path: &Path) -> Result<HashMap<String, TaxonomyEntry>, std::io::Error> {
        use std::io::BufRead;

        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut taxonomy = HashMap::new();

        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;

            // Skip header line
            if line_num == 0 {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 6 {
                log::warn!("Skipping malformed Silva line {}: insufficient fields", line_num + 1);
                continue;
            }

            let accession = fields[0].to_string();
            let path_str = fields[3];
            let organism_name = fields[4].to_string();
            let taxid = fields[5].to_string();

            // Parse taxonomy path (semicolon-separated)
            let tax_levels: Vec<&str> = path_str.split(';').collect();

            // Map to standard taxonomy fields
            // Silva: Kingdom → ... → Genus (variable depth)
            // Standard: Superkingdom, Phylum, Class, Order, Family, Genus, Species
            let superkingdom = tax_levels.get(0).unwrap_or(&"UNKNOWN").trim().to_string();
            let phylum = tax_levels.get(1).unwrap_or(&"UNKNOWN").trim().to_string();
            let class = tax_levels.get(2).unwrap_or(&"UNKNOWN").trim().to_string();
            let order = tax_levels.get(3).unwrap_or(&"UNKNOWN").trim().to_string();
            let family = tax_levels.get(4).unwrap_or(&"UNKNOWN").trim().to_string();
            let genus = tax_levels.get(5).unwrap_or(&"UNKNOWN").trim().to_string();

            let entry = TaxonomyEntry {
                tax_id: taxid,
                species: organism_name,
                genus,
                family,
                order,
                class,
                phylum,
                clade: String::new(), // Not in Silva
                superkingdom,
                subspecies: String::new(),
                species_subgroup: String::new(),
                species_group: String::new(),
            };

            taxonomy.insert(accession, entry);
        }

        Ok(taxonomy)
    }
}

/// Represents the classification result for a single ASV
#[derive(Debug, Clone)]
pub struct AsvClassification {
    pub asv_id: String,
    pub asv_header: String,
    pub hit_reference_id: String,
    pub abundance: f64,
    pub best_hit_tax_id: Option<String>,
    pub identity: Option<f64>,
    pub nm: Option<usize>,
    pub taxonomy: Option<TaxonomyAssignment>,
}

/// Represents a taxonomic assignment with UNCLASSIFIED markers
#[derive(Debug, Clone)]
pub struct TaxonomyAssignment {
    pub tax_id: String,
    pub species: String,
    pub genus: String,
    pub family: String,
    pub order: String,
    pub class: String,
    pub phylum: String,
    pub clade: String,
    pub superkingdom: String,
    pub subspecies: String,
    pub species_subgroup: String,
    pub species_group: String,
}

impl TaxonomyAssignment {
    /// Create a taxonomy assignment based on identity thresholds
    pub fn from_taxonomy_entry(
        entry: &TaxonomyEntry,
        identity: f64,
        species_threshold: f64,
        genus_threshold: f64,
        asv_header: &str,
    ) -> Self {
        let unclassified_marker = format!("UNCLASSIFIED-({})", asv_header);

        if identity >= species_threshold {
            // Species-level classification
            Self {
                tax_id: entry.tax_id.clone(),
                species: entry.species.clone(),
                genus: entry.genus.clone(),
                family: entry.family.clone(),
                order: entry.order.clone(),
                class: entry.class.clone(),
                phylum: entry.phylum.clone(),
                clade: entry.clade.clone(),
                superkingdom: entry.superkingdom.clone(),
                subspecies: entry.subspecies.clone(),
                species_subgroup: entry.species_subgroup.clone(),
                species_group: entry.species_group.clone(),
            }
        } else if identity >= genus_threshold {
            // Genus-level classification
            Self {
                tax_id: entry.tax_id.clone(),
                species: unclassified_marker.clone(),
                genus: entry.genus.clone(),
                family: entry.family.clone(),
                order: entry.order.clone(),
                class: entry.class.clone(),
                phylum: entry.phylum.clone(),
                clade: entry.clade.clone(),
                superkingdom: entry.superkingdom.clone(),
                subspecies: String::new(),
                species_subgroup: String::new(),
                species_group: String::new(),
            }
        } else if identity >= 75.0 {
            // Family-level classification
            Self {
                tax_id: entry.tax_id.clone(),
                species: unclassified_marker.clone(),
                genus: unclassified_marker.clone(),
                family: entry.family.clone(),
                order: entry.order.clone(),
                class: entry.class.clone(),
                phylum: entry.phylum.clone(),
                clade: entry.clade.clone(),
                superkingdom: entry.superkingdom.clone(),
                subspecies: String::new(),
                species_subgroup: String::new(),
                species_group: String::new(),
            }
        }
        else {
            // Below genus threshold - unclassified
            Self {
                tax_id: entry.tax_id.clone(),
                species: unclassified_marker.clone(),
                genus: unclassified_marker.clone(),
                family: unclassified_marker.clone(),
                order: unclassified_marker.clone(),
                class: unclassified_marker.clone(),
                phylum: unclassified_marker.clone(),
                clade: unclassified_marker.clone(),
                superkingdom: unclassified_marker,
                subspecies: String::new(),
                species_subgroup: String::new(),
                species_group: String::new(),
            }
        }
    }
}

/// Extract tax_id from EMU database FASTA header
/// Format: >2420510:emu_db:1 [...]
pub fn extract_tax_id_from_header(header: &str) -> Option<String> {
    let header = header.trim_start_matches('>');
    header.split(':').next().map(|s| s.to_string())
}

/// Extract accession from Silva database FASTA header
/// Format: >AY846372.1.1779 Eukaryota;...
/// Returns: AY846372
pub fn extract_silva_accession_from_header(header: &str) -> Option<String> {
    let header = header.trim_start_matches('>');
    // Split by space first to get the accession part
    let accession_part = header.split_whitespace().next()?;
    // Split by '.' and take the first token
    accession_part.split('.').next().map(|s| s.to_string())
}

/// Write species-level taxonomy abundance table to TSV file
pub fn write_species_abundance(
    classifications: &[AsvClassification],
    output_path: &Path,
) -> std::io::Result<()> {
    let mut file = File::create(output_path)?;

    // Write header
    writeln!(
        file,
        "tax_id\tabundance\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom"
    )?;

    // Aggregate abundances by taxonomy
    let mut taxonomy_abundances: HashMap<String, (TaxonomyAssignment, f64)> = HashMap::new();

    for classification in classifications {
        if let Some(ref taxonomy) = classification.taxonomy {
            // Create a unique key for this taxonomic assignment (species-level)
            let key = format!(
                "{}|{}|{}|{}|{}|{}|{}|{}|{}",
                taxonomy.species,
                taxonomy.genus,
                taxonomy.family,
                taxonomy.order,
                taxonomy.class,
                taxonomy.phylum,
                taxonomy.clade,
                taxonomy.superkingdom,
                taxonomy.tax_id
            );

            taxonomy_abundances
                .entry(key)
                .and_modify(|(_, abundance)| *abundance += classification.abundance)
                .or_insert((taxonomy.clone(), classification.abundance));
        }
    }

    // Sort by abundance (descending)
    let mut sorted_taxa: Vec<_> = taxonomy_abundances.into_iter().collect();
    sorted_taxa.sort_by(|a, b| b.1 .1.partial_cmp(&a.1 .1).unwrap());

    // Write sorted taxonomy entries
    for (_, (taxonomy, abundance)) in sorted_taxa {
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            taxonomy.tax_id,
            abundance,
            taxonomy.species,
            taxonomy.genus,
            taxonomy.family,
            taxonomy.order,
            taxonomy.class,
            taxonomy.phylum,
            taxonomy.clade,
            taxonomy.superkingdom,
        )?;
    }

    Ok(())
}

/// Write genus-level taxonomy abundance table to TSV file
pub fn write_genus_abundance(
    classifications: &[AsvClassification],
    output_path: &Path,
) -> std::io::Result<()> {
    let mut file = File::create(output_path)?;

    // Write header
    writeln!(
        file,
        "abundance\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom"
    )?;

    // Aggregate abundances by genus
    let mut genus_abundances: HashMap<String, (String, String, String, String, String, String, String, f64)> = HashMap::new();

    for classification in classifications {
        if let Some(ref taxonomy) = classification.taxonomy {
            // Create a unique key for genus-level (genus + higher levels)
            let key = format!(
                "{}|{}|{}|{}|{}|{}|{}",
                taxonomy.genus,
                taxonomy.family,
                taxonomy.order,
                taxonomy.class,
                taxonomy.phylum,
                taxonomy.clade,
                taxonomy.superkingdom
            );

            genus_abundances
                .entry(key)
                .and_modify(|(_, _, _, _, _, _, _, abundance)| *abundance += classification.abundance)
                .or_insert((
                    taxonomy.genus.clone(),
                    taxonomy.family.clone(),
                    taxonomy.order.clone(),
                    taxonomy.class.clone(),
                    taxonomy.phylum.clone(),
                    taxonomy.clade.clone(),
                    taxonomy.superkingdom.clone(),
                    classification.abundance
                ));
        }
    }

    // Sort by abundance (descending)
    let mut sorted_genera: Vec<_> = genus_abundances.into_iter().collect();
    sorted_genera.sort_by(|a, b| b.1 .7.partial_cmp(&a.1 .7).unwrap());

    // Write sorted genus entries
    for (_, (genus, family, order, class, phylum, clade, superkingdom, abundance)) in sorted_genera {
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            abundance,
            genus,
            family,
            order,
            class,
            phylum,
            clade,
            superkingdom,
        )?;
    }

    Ok(())
}

/// Write ASV mapping details to TSV file
pub fn write_asv_mappings(
    classifications: &[AsvClassification],
    output_path: &Path,
) -> std::io::Result<()> {
    let mut file = File::create(output_path)?;

    // Write header
    writeln!(
        file,
        "asv_header\tdepth\talignment_identity\tnumber_mismatches\ttax_id\tspecies\tgenus\treference"
    )?;

    for classification in classifications {
        // Extract depth from header (format: final_consensus_0_depth_42)
        let depth_str = extract_depth_string(&classification.asv_header);

        if let Some(ref taxonomy) = classification.taxonomy {
            if let Some(identity) = classification.identity {
                // Calculate approximate values
                // We don't have the exact mapping length and NM stored, so we'll need to add them
                // For now, write what we have
                writeln!(
                    file,
                    "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}",
                    classification.asv_header,
                    depth_str,
                    identity,
                    classification.nm.unwrap_or(0),
                    classification.best_hit_tax_id.as_ref().unwrap_or(&String::from("NA")),
                    taxonomy.species,
                    taxonomy.genus,
                    classification.hit_reference_id,
                )?;
            }
        } else {
            // Unclassified ASV
            writeln!(
                file,
                "{}\t{}\tNA\tNA\tNA\tNA\tUNCLASSIFIED\tUNCLASSIFIED",
                classification.asv_header,
                depth_str,
            )?;
        }
    }

    Ok(())
}


/// Load FASTA sequences using needletail
pub fn load_fasta_with_needletail(path: &Path) -> std::io::Result<Vec<(String, Vec<u8>)>> {
    let mut reader = needletail::parse_fastx_file(path)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

    let mut sequences = Vec::new();

    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let header = String::from_utf8_lossy(rec.id()).to_string();
        let seq = rec.seq().to_vec();
        sequences.push((format!(">{}", header), seq));
    }

    Ok(sequences)
}

/// Extract depth values from FASTA headers (format: >prefix_depth_N)
pub fn extract_depths_from_headers(sequences: &[(String, Vec<u8>)]) -> Vec<usize> {
    sequences.iter().map(|(header, _)| {
        // Parse depth from header like ">final_consensus_0_depth_42"
        let first_non_whitespace = header.split_whitespace().next().unwrap_or(header);
        first_non_whitespace.split('_')
            .last()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(1)
    }).collect()
}

pub fn extract_depth_string(sequence: &str) -> String {
    // Parse depth from header like ">final_consensus_0_depth_42"
    let first_non_whitespace = sequence.split_whitespace().next().unwrap_or(sequence);
    first_non_whitespace.split('_')
        .last()
        .unwrap_or("1")
        .to_string()
}