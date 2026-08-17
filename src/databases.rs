use std::path::Path;
use std::process::Command;
use crate::taxonomy::{self, Database};

const MARKER_FILE: &str = ".savont_db";
/// Flat keyword list used by the CLI `possible_values` validator.
pub const KEYWORDS: &[&str] = &["emu-1", "silva-138.2", "greengenes2-2024.09", "unite-10.0"];


/// Describes one versioned database: how to download it, load it, and extract
/// a lookup key from a minimap2 target-name string.
pub struct DatabaseDef {
    pub keyword:     &'static str,
    pub description: &'static str,
    pub download:    fn(dest: &Path) -> Result<(), String>,
    pub load:        fn(dir:  &Path) -> Result<Database, std::io::Error>,
    pub extract_key: fn(header: &str) -> Option<String>,
}

/// Every supported database version.  Adding a new one = one new row here.
pub const ALL: &[DatabaseDef] = &[
    DatabaseDef {
        keyword:     KEYWORDS[0],
        description: "EMU default 16S rRNA database",
        download:    download_emu,
        load:        Database::load_emu,
        extract_key: taxonomy::extract_tax_id_from_header,
    },
    DatabaseDef {
        keyword:     KEYWORDS[1],
        description: "SILVA SSU Ref NR99 v138.2",
        download:    download_silva,
        load:        Database::load_silva,
        extract_key: taxonomy::extract_silva_accession_from_header,
    },
    DatabaseDef {
        keyword:     KEYWORDS[2],
        description: "GreenGenes2 2024.09 species-level trainset from DADA2",
        download:    download_gg2,
        load:        Database::load_gg2,
        extract_key: taxonomy::extract_gg2_key_from_header,
    },
    DatabaseDef {
        keyword:     KEYWORDS[3],
        description: "UNITE ITS 10.0 — dynamic divergence thresholding, all eukaryotes",
        download:    download_unite,
        load:        Database::load_unite,
        extract_key: taxonomy::extract_unite_key_from_header,
    },
];


/// Look up a database definition by keyword.
pub fn find(keyword: &str) -> Option<&'static DatabaseDef> {
    ALL.iter().find(|d| d.keyword == keyword)
}

/// Human-readable list of all available keywords.
pub fn keyword_list() -> String {
    KEYWORDS.join(", ")
}

// ── marker file helpers ──────────────────────────────────────────────────────

/// Write a `.savont_db` marker so the directory is self-identifying.
pub fn write_marker(dir: &Path, keyword: &str) -> std::io::Result<()> {
    std::fs::write(dir.join(MARKER_FILE), keyword)
}

/// Read the `.savont_db` marker, if present.
pub fn read_marker(dir: &Path) -> Option<String> {
    std::fs::read_to_string(dir.join(MARKER_FILE))
        .ok()
        .map(|s| s.trim().to_string())
}

// ── database loader ──────────────────────────────────────────────────────────

/// Auto-detect the database type from the directory (marker file, then
/// basename), look it up in the registry, and load it.
pub fn load_database(dir: &Path) -> Result<Database, std::io::Error> {
    let keyword = read_marker(dir)
        .or_else(|| {
            dir.file_name()
               .and_then(|n| n.to_str())
               .map(|s| s.to_string())
        })
        .ok_or_else(|| std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("Cannot determine database type for '{}'", dir.display()),
        ))?;

    let def = find(&keyword).ok_or_else(|| std::io::Error::new(
        std::io::ErrorKind::InvalidInput,
        format!(
            "Unknown database keyword '{}'. Available: {}",
            keyword,
            keyword_list()
        ),
    ))?;

    log::info!("Detected database type '{}' for {}", keyword, dir.display());
    (def.load)(dir)
}

// ── download functions ───────────────────────────────────────────────────────

fn download_emu(dest: &Path) -> Result<(), String> {
    log::info!("Downloading EMU database...");
    let tar = dest.join("emu_default.tar.gz");

    let s = Command::new("wget")
        .arg("--content-disposition")
        .arg("https://osf.io/8qcwd/download")
        .arg("-O").arg(&tar)
        .status()
        .map_err(|e| format!("wget failed: {}", e))?;
    if !s.success() { return Err("wget returned non-zero for EMU download".into()); }

    // Extract into dest/ — tar creates dest/emu_default/
    let s = Command::new("tar")
        .arg("-xzf").arg(&tar)
        .arg("-C").arg(dest)
        .status()
        .map_err(|e| format!("tar failed: {}", e))?;
    if !s.success() { return Err("tar returned non-zero for EMU extraction".into()); }

    std::fs::remove_file(&tar).ok();

    // Move contents of dest/emu_default/ up into dest/ and remove the subdir.
    let subdir = dest.join("emu_default");
    for entry in std::fs::read_dir(&subdir)
        .map_err(|e| format!("Failed to read emu_default/: {}", e))?
    {
        let entry = entry.map_err(|e| format!("Directory entry error: {}", e))?;
        let target = dest.join(entry.file_name());
        std::fs::rename(entry.path(), &target)
            .map_err(|e| format!("mv {:?} -> {:?} failed: {}", entry.path(), target, e))?;
    }
    std::fs::remove_dir(&subdir).ok();

    Ok(())
}

fn download_silva(dest: &Path) -> Result<(), String> {
    log::info!("Downloading SILVA database...");
    let fasta_url = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz";
    let tax_url   = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.2.txt.gz";

    let s = Command::new("wget").arg(fasta_url).arg("-P").arg(dest)
        .status().map_err(|e| format!("wget failed: {}", e))?;
    if !s.success() { return Err("wget returned non-zero for SILVA FASTA".into()); }

    let s = Command::new("wget").arg(tax_url).arg("-P").arg(dest)
        .status().map_err(|e| format!("wget failed: {}", e))?;
    if !s.success() { return Err("wget returned non-zero for SILVA taxonomy".into()); }

    // Decompress the taxonomy file so load_silva can read it as plain text.
    Command::new("gzip")
        .arg("-d")
        .arg(dest.join("taxmap_slv_ssu_ref_nr_138.2.txt.gz"))
        .status()
        .map_err(|e| format!("gzip failed: {}", e))?;

    Ok(())
}

fn download_gg2(dest: &Path) -> Result<(), String> {
    log::info!("Downloading GreenGenes2 2024.09 species-level database...");
    let url = "https://zenodo.org/records/14169078/files/gg2_2024_09_toSpecies_trainset.fa.gz";

    let s = Command::new("wget").arg(url).arg("-P").arg(dest)
        .status().map_err(|e| format!("wget failed: {}", e))?;
    if !s.success() { return Err("wget returned non-zero for GreenGenes2 download".into()); }

    Ok(())
}

fn download_unite(dest: &Path) -> Result<(), String> {
    log::info!("Downloading UNITE ITS 10.0 database (dynamic divergence thresholding, all eukaryotes)...");
    let url = "https://s3.hpc.ut.ee/plutof-public/original/e861a3d6-54f4-42dc-882a-5f129beac39a.tgz";
    let tar = dest.join("unite-10.0.tgz");

    let s = Command::new("wget")
        .arg(url)
        .arg("-O").arg(&tar)
        .status()
        .map_err(|e| format!("wget failed: {}", e))?;
    if !s.success() { return Err("wget returned non-zero for UNITE download".into()); }

    let s = Command::new("tar")
        .arg("-xzf").arg(&tar)
        .arg("-C").arg(dest)
        .status()
        .map_err(|e| format!("tar failed: {}", e))?;
    if !s.success() { return Err("tar returned non-zero for UNITE extraction".into()); }

    std::fs::remove_file(&tar).ok();

    // Find the extracted .fasta file and preprocess it: minimap2 rejects reference names
    // longer than 255 chars, but many UNITE headers exceed that. We write a short-header
    // FASTA (SH accession only) and a companion taxonomy TSV so the full lineage is kept.
    let orig_fasta = std::fs::read_dir(dest)
        .map_err(|e| format!("Cannot read dest dir: {}", e))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .find(|p| p.extension().and_then(|x| x.to_str()) == Some("fasta"))
        .ok_or_else(|| "No .fasta file found after UNITE extraction".to_string())?;

    log::info!("Preprocessing UNITE FASTA (shortening headers, extracting taxonomy)...");
    preprocess_unite(&orig_fasta, dest)?;
    std::fs::remove_file(&orig_fasta).ok();

    log::info!("UNITE database ready in {}", dest.display());
    Ok(())
}

/// Rewrites a full UNITE FASTA into `unite_sequences.fasta` (SH-number headers only)
/// and writes `unite_taxonomy.tsv` (sh_id → taxonomy string), so minimap2 stays under
/// the 255-character reference name limit.
fn preprocess_unite(src: &Path, dest: &Path) -> Result<(), String> {
    use std::io::{BufRead, BufReader, Write, BufWriter};

    let reader = BufReader::new(
        std::fs::File::open(src).map_err(|e| format!("Cannot open {}: {}", src.display(), e))?
    );
    let mut fa = BufWriter::new(
        std::fs::File::create(dest.join("unite_sequences.fasta"))
            .map_err(|e| format!("Cannot create unite_sequences.fasta: {}", e))?
    );
    let mut tax = BufWriter::new(
        std::fs::File::create(dest.join("unite_taxonomy.tsv"))
            .map_err(|e| format!("Cannot create unite_taxonomy.tsv: {}", e))?
    );

    writeln!(tax, "sh_id\ttaxonomy").map_err(|e| e.to_string())?;

    for line in reader.lines() {
        let line = line.map_err(|e| format!("IO error reading FASTA: {}", e))?;
        if line.starts_with('>') {
            let header = &line[1..];
            // Header: Species|Accession|SH_id|type|k__X;p__X;...;s__X
            let parts: Vec<&str> = header.splitn(5, '|').collect();
            if parts.len() == 5 {
                let sh_id  = parts[2];
                let tax_str = parts[4];
                writeln!(fa,  ">{}", sh_id).map_err(|e| e.to_string())?;
                writeln!(tax, "{}\t{}", sh_id, tax_str).map_err(|e| e.to_string())?;
            } else {
                // Malformed header — truncate to 255 chars as a fallback
                let id = &header[..header.len().min(255)];
                writeln!(fa, ">{}", id).map_err(|e| e.to_string())?;
            }
        } else {
            writeln!(fa, "{}", line).map_err(|e| e.to_string())?;
        }
    }

    Ok(())
}
