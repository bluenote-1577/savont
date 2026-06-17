use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use flate2::read::GzDecoder;

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
    /// Extracts the taxonomy-map lookup key from a minimap2 target-name string.
    pub extract_key: fn(&str) -> Option<String>,
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
            extract_key: extract_tax_id_from_header,
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
            extract_key: extract_silva_accession_from_header,
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

    /// Load GTDB r232 SSU database from a directory containing ssu_all_r232.fna.gz
    pub fn load_gtdb(db_dir: &Path) -> Result<Self, std::io::Error> {
        let fasta_path = std::fs::read_dir(db_dir)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .find(|p| {
                p.file_name()
                    .and_then(|n| n.to_str())
                    .map_or(false, |n| n.ends_with(".fna.gz") || n.ends_with(".fna") || n.ends_with(".fa.gz") || n.ends_with(".fasta.gz"))
            })
            .ok_or_else(|| std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("No .fna.gz FASTA file found in {}", db_dir.display()),
            ))?;

        log::info!("Loading GTDB taxonomy from FASTA headers: {}", fasta_path.display());
        let taxonomy = Self::load_gtdb_taxonomy_from_fasta(&fasta_path)?;
        log::info!("Loaded {} GTDB taxonomy entries", taxonomy.len());

        Ok(Database { fasta_path, taxonomy, extract_key: extract_gtdb_key_from_header })
    }

    /// Parse taxonomy from GTDB FASTA headers.
    /// Header format: >REF_NAME d__Domain;p__Phylum;...;s__Genus species [location=...] ...
    fn load_gtdb_taxonomy_from_fasta(path: &Path) -> Result<HashMap<String, TaxonomyEntry>, std::io::Error> {
        let file = File::open(path)?;
        let reader: Box<dyn BufRead> = if path.to_str().map_or(false, |s| s.ends_with(".gz")) {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut taxonomy = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            if !line.starts_with('>') {
                continue;
            }
            let header = &line[1..];
            let mut tokens = header.splitn(2, ' ');
            let ref_name = match tokens.next() {
                Some(n) if !n.is_empty() => n.to_string(),
                _ => continue,
            };
            let rest = tokens.next().unwrap_or("");

            // Taxonomy string ends at first " [" annotation
            let tax_str = if let Some(idx) = rest.find(" [") {
                &rest[..idx]
            } else {
                rest.trim()
            };

            let mut superkingdom = String::new();
            let mut phylum = String::new();
            let mut class = String::new();
            let mut order = String::new();
            let mut family = String::new();
            let mut genus = String::new();
            let mut species = String::new();

            for level in tax_str.split(';') {
                let level = level.trim();
                if let Some(val) = level.strip_prefix("d__") {
                    superkingdom = val.to_string();
                } else if let Some(val) = level.strip_prefix("p__") {
                    phylum = val.to_string();
                } else if let Some(val) = level.strip_prefix("c__") {
                    class = val.to_string();
                } else if let Some(val) = level.strip_prefix("o__") {
                    order = val.to_string();
                } else if let Some(val) = level.strip_prefix("f__") {
                    family = val.to_string();
                } else if let Some(val) = level.strip_prefix("g__") {
                    genus = val.to_string();
                } else if let Some(val) = level.strip_prefix("s__") {
                    species = val.to_string();
                }
            }

            let entry = TaxonomyEntry {
                tax_id: ref_name.clone(),
                species,
                genus,
                family,
                order,
                class,
                phylum,
                clade: String::new(),
                superkingdom,
                subspecies: String::new(),
                species_subgroup: String::new(),
                species_group: String::new(),
            };
            taxonomy.insert(ref_name, entry);
        }

        Ok(taxonomy)
    }

    /// Load GreenGenes2 database from a directory containing gg2_*.fa.gz
    /// Header format: >d__Domain;p__Phylum;...;s__epithet;  (taxonomy IS the header)
    pub fn load_gg2(db_dir: &Path) -> Result<Self, std::io::Error> {
        let fasta_path = std::fs::read_dir(db_dir)?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .find(|p| {
                p.file_name()
                    .and_then(|n| n.to_str())
                    .map_or(false, |n| {
                        n.ends_with(".fa.gz") || n.ends_with(".fasta.gz") || n.ends_with(".fa")
                    })
            })
            .ok_or_else(|| std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("No .fa.gz file found in {}", db_dir.display()),
            ))?;

        log::info!("Loading GreenGenes2 taxonomy from FASTA headers: {}", fasta_path.display());
        let taxonomy = Self::load_gg2_taxonomy_from_fasta(&fasta_path)?;
        log::info!("Loaded {} GreenGenes2 taxonomy entries", taxonomy.len());

        Ok(Database { fasta_path, taxonomy, extract_key: extract_gg2_key_from_header })
    }

    /// Parse GreenGenes2 taxonomy from FASTA headers.
    /// Header format: >d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__epithet;
    /// The full header string (without '>') is used as the lookup key.
    fn load_gg2_taxonomy_from_fasta(path: &Path) -> Result<HashMap<String, TaxonomyEntry>, std::io::Error> {
        let file = File::open(path)?;
        let reader: Box<dyn BufRead> = if path.to_str().map_or(false, |s| s.ends_with(".gz")) {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut taxonomy = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            if !line.starts_with('>') { continue; }

            // Key = full header without '>'
            let key = line[1..].trim().to_string();
            if key.is_empty() { continue; }

            let mut superkingdom = String::new();
            let mut phylum = String::new();
            let mut class = String::new();
            let mut order = String::new();
            let mut family = String::new();
            let mut genus = String::new();
            let mut species_epithet = String::new();

            for level in key.split(';') {
                let level = level.trim();
                if let Some(val) = level.strip_prefix("d__") {
                    superkingdom = val.to_string();
                } else if let Some(val) = level.strip_prefix("p__") {
                    phylum = val.to_string();
                } else if let Some(val) = level.strip_prefix("c__") {
                    class = val.to_string();
                } else if let Some(val) = level.strip_prefix("o__") {
                    order = val.to_string();
                } else if let Some(val) = level.strip_prefix("f__") {
                    family = val.to_string();
                } else if let Some(val) = level.strip_prefix("g__") {
                    genus = val.to_string();
                } else if let Some(val) = level.strip_prefix("s__") {
                    species_epithet = val.to_string();
                }
            }

            // Build full species name from genus + epithet when both present
            let species = if !genus.is_empty() && !species_epithet.is_empty() {
                format!("{} {}", genus, species_epithet)
            } else {
                species_epithet
            };

            const UNANNOTATED: &str = "Greengenes_unannotated";
            let fill = |s: String| if s.is_empty() { UNANNOTATED.to_string() } else { s };

            let entry = TaxonomyEntry {
                tax_id: key.clone(),
                species:      fill(species),
                genus:        fill(genus),
                family:       fill(family),
                order:        fill(order),
                class:        fill(class),
                phylum:       fill(phylum),
                clade: String::new(),
                superkingdom: fill(superkingdom),
                subspecies: String::new(),
                species_subgroup: String::new(),
                species_group: String::new(),
            };
            taxonomy.insert(key, entry);
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
        detailed_unclassified: bool,
    ) -> Self {
        let unclassified_marker = if detailed_unclassified {
            format!("UNCLASSIFIED-({})", asv_header)
        } else {
            "UNCLASSIFIED".to_string()
        };

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
        } else if identity >= 86.5 {
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
        } else if identity >= 82.0 {
            // Order-level classification
            Self {
                tax_id: entry.tax_id.clone(),
                species: unclassified_marker.clone(),
                genus: unclassified_marker.clone(),
                family: unclassified_marker.clone(),
                order: entry.order.clone(),
                class: entry.class.clone(),
                phylum: entry.phylum.clone(),
                clade: entry.clade.clone(),
                superkingdom: entry.superkingdom.clone(),
                subspecies: String::new(),
                species_subgroup: String::new(),
                species_group: String::new(),
            }
        } else if identity >= 78.5 {
            // Class-level classification
            Self {
                tax_id: entry.tax_id.clone(),
                species: unclassified_marker.clone(),
                genus: unclassified_marker.clone(),
                family: unclassified_marker.clone(),
                order: unclassified_marker.clone(),
                class: entry.class.clone(),
                phylum: entry.phylum.clone(),
                clade: entry.clade.clone(),
                superkingdom: entry.superkingdom.clone(),
                subspecies: String::new(),
                species_subgroup: String::new(),
                species_group: String::new(),
            }
        } else if identity >= 75.0 {
            // Phylum-level classification
            Self {
                tax_id: entry.tax_id.clone(),
                species: unclassified_marker.clone(),
                genus: unclassified_marker.clone(),
                family: unclassified_marker.clone(),
                order: unclassified_marker.clone(),
                class: unclassified_marker.clone(),
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

/// Extract reference name from GTDB FASTA header (the full first token)
/// Format: >RS_GCF_002517985.1~NZ_NOCN01000152.1 d__Bacteria;...
/// Returns: RS_GCF_002517985.1~NZ_NOCN01000152.1
pub fn extract_gtdb_key_from_header(header: &str) -> Option<String> {
    let header = header.trim_start_matches('>');
    header.split_whitespace().next().map(|s| s.to_string())
}

/// Extract key from GreenGenes2 FASTA header.
/// Format: >d__Bacteria;p__...;s__epithet;  (no separate ref name — taxonomy IS the header)
/// Returns the full header string (trimmed), which is the taxonomy map key.
pub fn extract_gg2_key_from_header(header: &str) -> Option<String> {
    let header = header.trim_start_matches('>').trim();
    if header.is_empty() { None } else { Some(header.to_string()) }
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
        "abundance\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom"
    )?;

    // Aggregate abundances by taxonomy
    let mut taxonomy_abundances: HashMap<String, (TaxonomyAssignment, f64)> = HashMap::new();

    for classification in classifications {
        if let Some(ref taxonomy) = classification.taxonomy {
            // Create a unique key for this taxonomic assignment (species-level)
            let key = format!(
                "{}|{}|{}|{}|{}|{}|{}|{}",
                taxonomy.species,
                taxonomy.genus,
                taxonomy.family,
                taxonomy.order,
                taxonomy.class,
                taxonomy.phylum,
                taxonomy.clade,
                taxonomy.superkingdom,
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
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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

    writeln!(
        file,
        "asv_header\tdepth\talignment_identity\tnumber_mismatches\ttax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tclade\tsuperkingdom\treference"
    )?;

    for classification in classifications {
        let depth_str = extract_depth_string(&classification.asv_header);

        if let Some(ref taxonomy) = classification.taxonomy {
            if let Some(identity) = classification.identity {
                writeln!(
                    file,
                    "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    classification.asv_header,
                    depth_str,
                    identity,
                    classification.nm.unwrap_or(0),
                    classification.best_hit_tax_id.as_ref().unwrap_or(&String::from("NA")),
                    taxonomy.species,
                    taxonomy.genus,
                    taxonomy.family,
                    taxonomy.order,
                    taxonomy.class,
                    taxonomy.phylum,
                    taxonomy.clade,
                    taxonomy.superkingdom,
                    classification.hit_reference_id,
                )?;
            }
        } else {
            writeln!(
                file,
                "{}\t{}\tNA\tNA\tNA\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED\tUNCLASSIFIED",
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

/// Extract total depth from a FASTA header token (supports both "42" and "42-15-27" formats).
fn parse_depth_token(token: &str) -> usize {
    let depth: usize = token.split('-')
        .filter_map(|s| s.parse::<usize>().ok())
        .sum();
    depth.max(1)
}

/// Extract depth values from FASTA headers (format: >prefix_depth_N or >prefix_depth_N-M-K).
pub fn extract_depths_from_headers(sequences: &[(String, Vec<u8>)]) -> Vec<usize> {
    sequences.iter().map(|(header, _)| {
        let first_token = header.split_whitespace().next().unwrap_or(header);
        let depth_token = first_token.split('_').last().unwrap_or("1");
        parse_depth_token(depth_token)
    }).collect()
}

pub fn extract_depth_string(sequence: &str) -> String {
    let first_token = sequence.split_whitespace().next().unwrap_or(sequence);
    first_token.split('_').last().unwrap_or("1").to_string()
}