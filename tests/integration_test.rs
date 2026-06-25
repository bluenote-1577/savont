use assert_cmd::Command;
use minimap2;
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

const REF_FA: &str = "tests/data/zymo_ref_asvs.fa.gz";
const READS_FQ: &str = "tests/data/ont_zymo_1000.trimmed.fq.gz";
const READS_FQ_2: &str = "tests/data/ont_zymo_1000_2.trimmed.fq.gz";
const EMU_DB_PATH: &str = "tests/data/emu-1";
const SILVA_DB_PATH: &str = "tests/data/silva-138.2";
const GG2_DB_PATH: &str = "tests/data/greengenes2-2024.09";

// ── shared helpers ────────────────────────────────────────────────────────────

/// Run `savont asv` on a reads file and return the output TempDir.
fn run_asv(reads: &str) -> TempDir {
    let tmp = TempDir::new().unwrap();
    Command::cargo_bin("savont")
        .unwrap()
        .args([
            "asv", reads,
            "-o", tmp.path().to_str().unwrap(),
            "-t", "4",
            "--min-cluster-size", "5",
        ])
        .assert()
        .success();
    tmp
}

/// Return the EMU database path if it is already present at `tests/data/emu-1/`.
/// Returns None so the calling test soft-skips. To enable these tests, run:
///   savont download --location tests/data --dbs emu-1
fn emu_db() -> Option<PathBuf> {
    let db = Path::new(EMU_DB_PATH);
    if db.join("species_taxid.fasta").exists() {
        Some(db.to_path_buf())
    } else {
        eprintln!(
            "Skipping: EMU database not found at {}. \
             Run: savont download --location tests/data --dbs emu-1",
            EMU_DB_PATH
        );
        None
    }
}

/// Return the SILVA database path if it is already present at `tests/data/silva-138.2/`.
/// To enable these tests, run: savont download --location tests/data --dbs silva-138.2
fn silva_db() -> Option<PathBuf> {
    let db = Path::new(SILVA_DB_PATH);
    let has_fasta = fs::read_dir(db).ok()?.filter_map(|e| e.ok())
        .any(|e| e.file_name().to_str()
            .map_or(false, |n| n.ends_with(".fasta.gz") || n.ends_with(".fasta")));
    if has_fasta {
        Some(db.to_path_buf())
    } else {
        eprintln!(
            "Skipping: SILVA database not found at {}. \
             Run: savont download --location tests/data --dbs silva-138.2",
            SILVA_DB_PATH
        );
        None
    }
}

/// Return the GreenGenes2 database path if it is already present at `tests/data/greengenes2-2024.09/`.
/// To enable these tests, run: savont download --location tests/data --dbs greengenes2-2024.09
fn gg2_db() -> Option<PathBuf> {
    let db = Path::new(GG2_DB_PATH);
    let has_fasta = fs::read_dir(db).ok()?.filter_map(|e| e.ok())
        .any(|e| e.file_name().to_str()
            .map_or(false, |n| n.ends_with(".fa.gz") || n.ends_with(".fa")));
    if has_fasta {
        Some(db.to_path_buf())
    } else {
        eprintln!(
            "Skipping: GreenGenes2 database not found at {}. \
             Run: savont download --location tests/data --dbs greengenes2-2024.09",
            GG2_DB_PATH
        );
        None
    }
}

/// Run `savont asv` on the bundled 1000-read zymo dataset, then verify that
/// every output ASV aligns to the zymo reference ASVs with zero mismatches,
/// using the minimap2 Rust bindings (no external binary required).
#[test]
fn test_asv_generation_and_perfect_alignment() {
    let tmp = TempDir::new().unwrap();
    let out_dir = tmp.path().to_str().unwrap();

    // --- 1. run savont asv ---
    Command::cargo_bin("savont")
        .unwrap()
        .args([
            "asv", READS_FQ,
            "-o", out_dir,
            "-t", "4",
        ])
        .assert()
        .success();

    let asv_fasta = tmp.path().join("final_asvs.fasta");
    assert!(asv_fasta.exists(), "final_asvs.fasta was not created");

    // --- 2. load ASV sequences ---
    let asvs = savont::taxonomy::load_fasta_with_needletail(&asv_fasta)
        .expect("failed to read final_asvs.fasta");
    assert!(!asvs.is_empty(), "savont produced zero ASVs");

    // --- 3. build minimap2 index from reference, align each ASV ---
    let aligner = minimap2::Aligner::builder()
        .map_ont()
        .with_cigar()
        .with_index(REF_FA, None)
        .expect("failed to build minimap2 index from zymo reference");

    let mut mapped = 0usize;
    let mut imperfect: Vec<(String, i32)> = Vec::new();

    for (header, seq) in &asvs {
        let hits = aligner
            .map(seq, true, false, None, None, None)
            .expect("minimap2 mapping failed");

        // Keep only the primary hit (lowest NM / highest score comes first)
        let primary = hits
            .iter()
            .find(|h| h.is_primary);

        match primary {
            None => { /* unmapped – will be caught by the mapped == asvs.len() assert */ }
            Some(hit) => {
                mapped += 1;
                let nm = hit.alignment.as_ref().map(|a| a.nm).unwrap_or(i32::MAX);
                if nm > 0 {
                    imperfect.push((header.clone(), nm));
                }
            }
        }
    }

    assert!(mapped > 0, "no ASVs aligned to the zymo reference");
    assert_eq!(
        mapped,
        asvs.len(),
        "only {}/{} ASVs mapped to the zymo reference",
        mapped,
        asvs.len()
    );
    assert!(
        imperfect.is_empty(),
        "ASVs with NM > 0 (not perfect): {:?}",
        imperfect
    );
}

// ── database download + load tests ──────────────────────────────────────────
// These are skipped by default (`cargo test`) because they download large
// files from the internet.  Run them explicitly with:
//   cargo test -- --ignored

/// Download the EMU database, check the expected files are present, and load
/// the taxonomy via the registry auto-detection path.
#[test]
#[ignore]
fn test_download_and_load_emu_database() {
    let tmp = TempDir::new().unwrap();
    let loc = tmp.path().to_str().unwrap();

    Command::cargo_bin("savont")
        .unwrap()
        .args(["download", "--location", loc, "--dbs", "emu"])
        .assert()
        .success();

    let db_dir = tmp.path().join("emu");
    assert!(db_dir.join("species_taxid.fasta").exists(), "species_taxid.fasta missing");
    assert!(db_dir.join("taxonomy.tsv").exists(), "taxonomy.tsv missing");
    assert!(db_dir.join(".savont_db").exists(), "marker file missing");

    let db = savont::databases::load_database(&db_dir)
        .expect("Failed to load EMU database via registry");
    assert!(!db.taxonomy.is_empty(), "EMU taxonomy map is empty");
    assert!(db.fasta_path.exists(), "FASTA path does not exist on disk");
}

/// Download the SILVA database, check the expected files are present, and load
/// the taxonomy via the registry auto-detection path.
#[test]
#[ignore]
fn test_download_and_load_silva_database() {
    let tmp = TempDir::new().unwrap();
    let loc = tmp.path().to_str().unwrap();

    Command::cargo_bin("savont")
        .unwrap()
        .args(["download", "--location", loc, "--dbs", "silva-138.2"])
        .assert()
        .success();

    let db_dir = tmp.path().join("silva-138.2");
    let has_fasta = fs::read_dir(&db_dir).unwrap().filter_map(|e| e.ok())
        .any(|e| e.file_name().to_str()
            .map_or(false, |n| n.ends_with(".fasta.gz") || n.ends_with(".fasta")));
    assert!(has_fasta, "no FASTA file found after download");

    let has_taxmap = fs::read_dir(&db_dir).unwrap().filter_map(|e| e.ok())
        .any(|e| e.file_name().to_str()
            .map_or(false, |n| n.starts_with("taxmap_") && n.ends_with(".txt")));
    assert!(has_taxmap, "no taxmap file found after download");

    let db = savont::databases::load_database(&db_dir)
        .expect("Failed to load SILVA database via registry");
    assert!(!db.taxonomy.is_empty(), "SILVA taxonomy map is empty");
}

/// Download the GTDB r232 SSU database, check the expected file is present,
/// and load the taxonomy via the registry auto-detection path.
#[test]
#[ignore]
fn test_download_and_load_gtdb_database() {
    let tmp = TempDir::new().unwrap();
    let loc = tmp.path().to_str().unwrap();

    Command::cargo_bin("savont")
        .unwrap()
        .args(["download", "--location", loc, "--dbs", "gtdb-r232"])
        .assert()
        .success();

    let db_dir = tmp.path().join("gtdb-r232");
    let has_fna = fs::read_dir(&db_dir).unwrap().filter_map(|e| e.ok())
        .any(|e| e.file_name().to_str().map_or(false, |n| n.ends_with(".fna.gz")));
    assert!(has_fna, "no .fna.gz file found after download");

    let db = savont::databases::load_database(&db_dir)
        .expect("Failed to load GTDB database via registry");
    assert!(!db.taxonomy.is_empty(), "GTDB taxonomy map is empty");

    let bad: Vec<_> = db.taxonomy.values()
        .filter(|e| e.superkingdom.is_empty())
        .take(5).map(|e| e.tax_id.clone()).collect();
    assert!(bad.is_empty(), "GTDB entries missing superkingdom: {:?}", bad);
}

/// Parse a small hand-crafted GTDB FASTA header without hitting the network,
/// to confirm the taxonomy parser is correct.
#[test]
fn test_gtdb_taxonomy_parser_unit() {
    use std::io::Write;

    let tmp = TempDir::new().unwrap();

    // Write a two-entry mock GTDB FASTA into a temp directory
    let fna_path = tmp.path().join("mock_gtdb.fna");
    {
        let mut f = fs::File::create(&fna_path).unwrap();
        writeln!(
            f,
            ">RS_GCF_000001405.40~NC_000001.11 d__Bacteria;p__Pseudomonadota;\
c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;\
g__Escherichia;s__Escherichia coli [location=1..1500] [ssu_len=1500]"
        )
        .unwrap();
        writeln!(f, "ACGT").unwrap();
        writeln!(
            f,
            ">GB_GCA_000007185.1~AE017221.1 d__Archaea;p__Thermoproteota;\
c__Thermoprotei;o__Thermoproteales;f__Thermoproteaceae;\
g__Thermoproteus;s__Thermoproteus tenax [location=1..1200] [ssu_len=1200]"
        )
        .unwrap();
        writeln!(f, "TTTT").unwrap();
    }

    let db = savont::taxonomy::Database::load_gtdb(tmp.path())
        .expect("load_gtdb failed on mock file");

    assert_eq!(db.taxonomy.len(), 2);

    let ecoli = db
        .taxonomy
        .get("RS_GCF_000001405.40~NC_000001.11")
        .expect("E. coli entry missing");
    assert_eq!(ecoli.superkingdom, "Bacteria");
    assert_eq!(ecoli.phylum, "Pseudomonadota");
    assert_eq!(ecoli.class, "Gammaproteobacteria");
    assert_eq!(ecoli.order, "Enterobacterales");
    assert_eq!(ecoli.family, "Enterobacteriaceae");
    assert_eq!(ecoli.genus, "Escherichia");
    assert_eq!(ecoli.species, "Escherichia coli");

    let archaea = db
        .taxonomy
        .get("GB_GCA_000007185.1~AE017221.1")
        .expect("Archaea entry missing");
    assert_eq!(archaea.superkingdom, "Archaea");
    assert_eq!(archaea.genus, "Thermoproteus");
    assert_eq!(archaea.species, "Thermoproteus tenax");
}

// ── merge tests ───────────────────────────────────────────────────────────────

/// Run `savont asv` on both zymo replicates, merge them, and verify the ASV
/// feature table and representative sequences are structurally correct.
/// This test does not require a reference database.
#[test]
fn test_merge_feature_table() {
    let tmp1 = run_asv(READS_FQ);
    let tmp2 = run_asv(READS_FQ_2);
    let merge_out = TempDir::new().unwrap();

    Command::cargo_bin("savont")
        .unwrap()
        .args([
            "export",
            "-i", tmp1.path().to_str().unwrap(), tmp2.path().to_str().unwrap(),
            "-o", merge_out.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    // --- feature table structure ---
    let ft_path = merge_out.path().join("merged_feature_table.tsv");
    assert!(ft_path.exists(), "merged_feature_table.tsv not created");

    let content = fs::read_to_string(&ft_path).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert!(lines.len() >= 2, "feature table has fewer than 2 lines");
    assert!(lines[0].starts_with("#OTU ID\t"), "header row should start with #OTU ID");

    let n_cols = lines[0].split('\t').count();
    assert_eq!(n_cols, 3, "expected #OTU ID + 2 sample columns, got {}", n_cols);

    // All data rows: 3 tab-separated fields, last two are non-negative integers.
    for line in &lines[1..] {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 3, "data row has wrong column count: {:?}", line);
        fields[1].parse::<u64>().expect("sample-1 count is not an integer");
        fields[2].parse::<u64>().expect("sample-2 count is not an integer");
    }

    // --- rep seqs exist and IDs match feature table ---
    let rs_path = merge_out.path().join("merged_rep_seqs.fasta");
    assert!(rs_path.exists(), "merged_rep_seqs.fasta not created");

    let ft_ids: HashSet<String> = lines[1..].iter()
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    let rs_ids: HashSet<String> = fs::read_to_string(&rs_path)
        .unwrap()
        .lines()
        .filter(|l| l.starts_with('>'))
        .map(|l| l[1..].split_whitespace().next().unwrap().to_string())
        .collect();
    assert_eq!(ft_ids, rs_ids, "feature table ASV IDs must match rep_seqs IDs");

    // --- at least some ASVs are present in both samples ---
    let shared = lines[1..].iter().filter(|l| {
        let f: Vec<&str> = l.split('\t').collect();
        f[1].parse::<u64>().unwrap_or(0) > 0 && f[2].parse::<u64>().unwrap_or(0) > 0
    }).count();
    assert!(shared > 0, "no ASVs found in both samples — merge may be broken");
}

/// Run `savont classify` (EMU) on both replicates, merge, and verify the
/// species-level count and relab tables are present with the correct structure.
/// Downloads EMU to tests/data/emu-1/ on first run; soft-skips if unavailable.
#[test]
fn test_merge_with_classify() {
    let Some(db_dir) = emu_db() else {
        eprintln!("Skipping test_merge_with_classify: EMU database unavailable");
        return;
    };
    let db_str = db_dir.to_str().unwrap();

    let tmp1 = run_asv(READS_FQ);
    let tmp2 = run_asv(READS_FQ_2);

    for tmp in [&tmp1, &tmp2] {
        Command::cargo_bin("savont")
            .unwrap()
            .args(["classify", "-i", tmp.path().to_str().unwrap(), "-d", db_str, "-t", "4"])
            .assert()
            .success();
    }

    let merge_out = TempDir::new().unwrap();
    Command::cargo_bin("savont")
        .unwrap()
        .args([
            "export",
            "-i", tmp1.path().to_str().unwrap(), tmp2.path().to_str().unwrap(),
            "-o", merge_out.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    // Feature table must exist with correct structure
    let ft_path = merge_out.path().join("merged_feature_table.tsv");
    assert!(ft_path.exists(), "merged_feature_table.tsv not created");
    let ft_content = fs::read_to_string(&ft_path).unwrap();
    let ft_lines: Vec<&str> = ft_content.lines().collect();
    assert!(ft_lines[0].starts_with("#OTU ID\t"), "feature table header malformed");
    assert_eq!(ft_lines[0].split('\t').count(), 3, "expected 2 sample columns");

    // ASV taxonomy file must exist with correct header
    let asv_tax_path = merge_out.path().join("merged_asv_taxonomy.tsv");
    assert!(asv_tax_path.exists(), "merged_asv_taxonomy.tsv not created");
    let asv_tax = fs::read_to_string(&asv_tax_path).unwrap();
    let asv_tax_lines: Vec<&str> = asv_tax.lines().collect();
    assert_eq!(asv_tax_lines[0], "Feature ID\tTaxon", "asv taxonomy header wrong");
    assert!(asv_tax_lines.len() > 1, "merged_asv_taxonomy.tsv has no data rows");
}

/// Run `savont sintax` (genus-level) on both replicates, merge, and verify
/// that the merge falls back to genus-level tables (no species_abundance.tsv).
/// Downloads EMU to tests/data/emu-1/ on first run; soft-skips if unavailable.
#[test]
fn test_merge_with_sintax() {
    let Some(db_dir) = emu_db() else {
        eprintln!("Skipping test_merge_with_sintax: EMU database unavailable");
        return;
    };
    let db_str = db_dir.to_str().unwrap();

    let tmp1 = run_asv(READS_FQ);
    let tmp2 = run_asv(READS_FQ_2);

    for tmp in [&tmp1, &tmp2] {
        Command::cargo_bin("savont")
            .unwrap()
            .args(["sintax", "-i", tmp.path().to_str().unwrap(), "-d", db_str, "-t", "4"])
            .assert()
            .success();
    }

    let merge_out = TempDir::new().unwrap();
    Command::cargo_bin("savont")
        .unwrap()
        .args([
            "export",
            "-i", tmp1.path().to_str().unwrap(), tmp2.path().to_str().unwrap(),
            "-o", merge_out.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    // sintax does not produce species_abundance.tsv — merged_species_counts.tsv should not exist
    assert!(
        !merge_out.path().join("merged_species_counts.tsv").exists(),
        "merged_species_counts.tsv should not exist when only sintax was run"
    );

    // ASV taxonomy file must exist
    assert!(
        merge_out.path().join("merged_asv_taxonomy.tsv").exists(),
        "merged_asv_taxonomy.tsv not created"
    );
}

// ── SILVA classify tests ──────────────────────────────────────────────────────

fn check_asv_mappings_columns(dir: &Path) {
    let path = dir.join("asv_mappings.tsv");
    assert!(path.exists(), "asv_mappings.tsv not created");
    let first_line = fs::read_to_string(&path).unwrap();
    let header = first_line.lines().next().unwrap();
    let cols: Vec<&str> = header.split('\t').collect();
    for expected in &["asv_header", "species", "genus", "family", "order", "class", "phylum", "superkingdom"] {
        assert!(cols.contains(expected), "asv_mappings.tsv missing column '{}'", expected);
    }
}

/// Run `savont classify` with the SILVA database on the Zymo dataset and verify
/// output structure. Soft-skips if the database is not present at tests/data/silva-138.2/.
#[test]
fn test_classify_with_silva() {
    let Some(db_dir) = silva_db() else { return; };
    let db_str = db_dir.to_str().unwrap();

    let tmp = run_asv(READS_FQ);

    Command::cargo_bin("savont")
        .unwrap()
        .args(["classify", "-i", tmp.path().to_str().unwrap(), "-d", db_str, "-t", "4"])
        .assert()
        .success();

    // Both abundance tables should exist
    assert!(tmp.path().join("species_abundance.tsv").exists(), "species_abundance.tsv missing");
    assert!(tmp.path().join("genus_abundance.tsv").exists(), "genus_abundance.tsv missing");

    // asv_mappings.tsv should have all lineage columns
    check_asv_mappings_columns(tmp.path());

    // At least one Zymo genus should appear in the species table
    let species = fs::read_to_string(tmp.path().join("species_abundance.tsv")).unwrap();
    let known = ["Listeria", "Pseudomonas", "Escherichia", "Salmonella", "Staphylococcus"];
    assert!(
        known.iter().any(|g| species.contains(g)),
        "none of the expected Zymo genera found in SILVA species table"
    );
}

/// Run `savont classify` (SILVA) on both replicates, merge, and verify the
/// species-level tables and ASV taxonomy file are structurally correct.
/// Soft-skips if the database is not present at tests/data/silva-138.2/.
#[test]
fn test_merge_with_silva() {
    let Some(db_dir) = silva_db() else { return; };
    let db_str = db_dir.to_str().unwrap();

    let tmp1 = run_asv(READS_FQ);
    let tmp2 = run_asv(READS_FQ_2);

    for tmp in [&tmp1, &tmp2] {
        Command::cargo_bin("savont")
            .unwrap()
            .args(["classify", "-i", tmp.path().to_str().unwrap(), "-d", db_str, "-t", "4"])
            .assert()
            .success();
    }

    let merge_out = TempDir::new().unwrap();
    Command::cargo_bin("savont")
        .unwrap()
        .args([
            "export",
            "-i", tmp1.path().to_str().unwrap(), tmp2.path().to_str().unwrap(),
            "-o", merge_out.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    // ASV taxonomy file: all hashes should have a lineage containing known ranks
    let asv_tax_path = merge_out.path().join("merged_asv_taxonomy.tsv");
    assert!(asv_tax_path.exists(), "merged_asv_taxonomy.tsv not created");
    let asv_tax = fs::read_to_string(&asv_tax_path).unwrap();
    let asv_tax_lines: Vec<&str> = asv_tax.lines().collect();
    assert_eq!(asv_tax_lines[0], "Feature ID\tTaxon", "asv taxonomy header wrong");
    assert!(asv_tax_lines.len() > 1, "merged_asv_taxonomy.tsv has no data rows");

    // Every classified row should have a semicolon-separated lineage (not bare "Unclassified")
    let unclassified_count = asv_tax_lines[1..].iter()
        .filter(|l| l.split('\t').nth(1).map_or(false, |t| !t.contains(';')))
        .count();
    let total = asv_tax_lines.len() - 1;
    assert!(
        unclassified_count < total,
        "all ASVs are bare Unclassified in merged_asv_taxonomy — lineage join may be broken"
    );

    // Feature table and rep seqs must have matching IDs
    let ft = fs::read_to_string(merge_out.path().join("merged_feature_table.tsv")).unwrap();
    let ft_ids: HashSet<String> = ft.lines().skip(1)
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    let tax_ids: HashSet<String> = asv_tax_lines[1..].iter()
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(ft_ids, tax_ids, "feature table IDs must match ASV taxonomy IDs");
}

// ── GreenGenes2 classify tests ────────────────────────────────────────────────

/// Run `savont classify` with the GreenGenes2 database on the Zymo dataset and
/// verify output structure. Soft-skips if not present at tests/data/greengenes2-2024.09/.
#[test]
fn test_classify_with_greengenes2() {
    let Some(db_dir) = gg2_db() else { return; };
    let db_str = db_dir.to_str().unwrap();

    let tmp = run_asv(READS_FQ);

    Command::cargo_bin("savont")
        .unwrap()
        .args(["classify", "-i", tmp.path().to_str().unwrap(), "-d", db_str, "-t", "4"])
        .assert()
        .success();

    assert!(tmp.path().join("species_abundance.tsv").exists(), "species_abundance.tsv missing");
    assert!(tmp.path().join("genus_abundance.tsv").exists(), "genus_abundance.tsv missing");

    check_asv_mappings_columns(tmp.path());

    // GreenGenes2 uses d__/p__ prefixes internally but stores plain names in TaxonomyEntry
    let species = fs::read_to_string(tmp.path().join("species_abundance.tsv")).unwrap();
    let known = ["Listeria", "Pseudomonas", "Escherichia", "Salmonella", "Staphylococcus"];
    assert!(
        known.iter().any(|g| species.contains(g)),
        "none of the expected Zymo genera found in GreenGenes2 species table"
    );
}

/// Run `savont classify` (GreenGenes2) on both replicates, merge, and verify
/// that classified lineages in merged_asv_taxonomy.tsv are non-trivial.
/// Soft-skips if not present at tests/data/greengenes2-2024.09/.
#[test]
fn test_merge_with_greengenes2() {
    let Some(db_dir) = gg2_db() else { return; };
    let db_str = db_dir.to_str().unwrap();

    let tmp1 = run_asv(READS_FQ);
    let tmp2 = run_asv(READS_FQ_2);

    for tmp in [&tmp1, &tmp2] {
        Command::cargo_bin("savont")
            .unwrap()
            .args(["classify", "-i", tmp.path().to_str().unwrap(), "-d", db_str, "-t", "4"])
            .assert()
            .success();
    }

    let merge_out = TempDir::new().unwrap();
    Command::cargo_bin("savont")
        .unwrap()
        .args([
            "export",
            "-i", tmp1.path().to_str().unwrap(), tmp2.path().to_str().unwrap(),
            "-o", merge_out.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    assert!(
        merge_out.path().join("merged_asv_taxonomy.tsv").exists(),
        "merged_asv_taxonomy.tsv not created"
    );

    // ASV taxonomy: IDs must match feature table, classified lineages must have semicolons
    let ft = fs::read_to_string(merge_out.path().join("merged_feature_table.tsv")).unwrap();
    let ft_ids: HashSet<String> = ft.lines().skip(1)
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();

    let asv_tax = fs::read_to_string(merge_out.path().join("merged_asv_taxonomy.tsv")).unwrap();
    let asv_tax_lines: Vec<&str> = asv_tax.lines().collect();
    let tax_ids: HashSet<String> = asv_tax_lines[1..].iter()
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(ft_ids, tax_ids, "feature table IDs must match ASV taxonomy IDs");

    let classified = asv_tax_lines[1..].iter()
        .filter(|l| l.split('\t').nth(1).map_or(false, |t| t.contains(';')))
        .count();
    assert!(classified > 0, "no classified lineages (with semicolons) in merged_asv_taxonomy");
}

// ── pooled-samples tests ──────────────────────────────────────────────────────

/// Run `savont asv --pooled-samples` on both zymo replicates and verify:
/// - FASTA headers carry dash-separated per-sample depths
/// - feature-table.tsv has exactly 2 sample columns with integer depths
#[test]
fn test_pooled_samples_asv() {
    let tmp = TempDir::new().unwrap();
    Command::cargo_bin("savont")
        .unwrap()
        .args([
            "asv", "--pooled-samples", READS_FQ, READS_FQ_2,
            "-o", tmp.path().to_str().unwrap(),
            "-t", "4",
            "--min-cluster-size", "5",
        ])
        .assert()
        .success();

    // FASTA headers must use dash-separated depths in pooled mode
    let fasta_path = tmp.path().join("final_asvs.fasta");
    assert!(fasta_path.exists(), "final_asvs.fasta not created");
    let fasta = fs::read_to_string(&fasta_path).unwrap();
    let headers: Vec<&str> = fasta.lines().filter(|l| l.starts_with('>')).collect();
    assert!(!headers.is_empty(), "no ASVs produced by --pooled-samples");
    for h in &headers {
        let depth_part = h.split("_depth_").nth(1)
            .unwrap_or_else(|| panic!("header missing _depth_: {}", h));
        assert!(
            depth_part.contains('-'),
            "pooled header should have dash-separated depths: {}", h
        );
    }

    // feature-table.tsv: header row + 2 sample columns + integer depths in each data row
    let ft_path = tmp.path().join("feature-table.tsv");
    assert!(ft_path.exists(), "feature-table.tsv not created");
    let ft = fs::read_to_string(&ft_path).unwrap();
    let ft_lines: Vec<&str> = ft.lines().collect();
    assert!(ft_lines[0].starts_with("#OTU ID\t"), "feature table header malformed");
    assert_eq!(
        ft_lines[0].split('\t').count(), 3,
        "expected #OTU ID + 2 sample columns in pooled feature table"
    );
    for line in &ft_lines[1..] {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 3, "data row has wrong column count: {:?}", line);
        fields[1].parse::<u64>().expect("sample-1 depth is not an integer");
        fields[2].parse::<u64>().expect("sample-2 depth is not an integer");
    }
}

/// Run `savont classify` (EMU) on a `--pooled-samples` run and verify that
/// species/genus abundance tables are produced and contain known Zymo genera.
/// Soft-skips if EMU database is unavailable.
#[test]
fn test_pooled_samples_classify() {
    let Some(db_dir) = emu_db() else {
        eprintln!("Skipping test_pooled_samples_classify: EMU database unavailable");
        return;
    };
    let db_str = db_dir.to_str().unwrap();

    let tmp = TempDir::new().unwrap();
    Command::cargo_bin("savont")
        .unwrap()
        .args([
            "asv", "--pooled-samples", READS_FQ, READS_FQ_2,
            "-o", tmp.path().to_str().unwrap(),
            "-t", "4",
            "--min-cluster-size", "5",
        ])
        .assert()
        .success();

    Command::cargo_bin("savont")
        .unwrap()
        .args(["classify", "-i", tmp.path().to_str().unwrap(), "-d", db_str, "-t", "4"])
        .assert()
        .success();

    let species_path = tmp.path().join("species_abundance.tsv");
    assert!(species_path.exists(), "species_abundance.tsv not created");
    let genus_path = tmp.path().join("genus_abundance.tsv");
    assert!(genus_path.exists(), "genus_abundance.tsv not created");

    // Both tables should have 2 sample columns (one per pooled input file)
    for path in [&species_path, &genus_path] {
        let content = fs::read_to_string(path).unwrap();
        let header = content.lines().next().unwrap_or("");
        // abundance + 2 sample names + taxonomy columns — at least 4 fields total
        assert!(
            header.split('\t').count() >= 4,
            "{} header should have ≥4 columns for pooled classify: {}",
            path.display(), header
        );
    }

    // At least one known Zymo genus must appear
    let species = fs::read_to_string(&species_path).unwrap();
    let known = ["Listeria", "Pseudomonas", "Escherichia", "Salmonella", "Staphylococcus"];
    assert!(
        known.iter().any(|g| species.contains(g)),
        "none of the expected Zymo genera found in pooled species_abundance.tsv"
    );
}

/// Run `savont export` on a `--pooled-samples` output directory and verify that
/// the merged feature table expands to 2 columns and IDs match rep seqs.
#[test]
fn test_pooled_samples_export() {
    let tmp = TempDir::new().unwrap();
    Command::cargo_bin("savont")
        .unwrap()
        .args([
            "asv", "--pooled-samples", READS_FQ, READS_FQ_2,
            "-o", tmp.path().to_str().unwrap(),
            "-t", "4",
            "--min-cluster-size", "5",
        ])
        .assert()
        .success();

    let export_out = TempDir::new().unwrap();
    Command::cargo_bin("savont")
        .unwrap()
        .args([
            "export",
            "-i", tmp.path().to_str().unwrap(),
            "-o", export_out.path().to_str().unwrap(),
        ])
        .assert()
        .success();

    // merged_feature_table.tsv must carry 2 sample columns (expanded from pooled input)
    let ft_path = export_out.path().join("merged_feature_table.tsv");
    assert!(ft_path.exists(), "merged_feature_table.tsv not created");
    let ft = fs::read_to_string(&ft_path).unwrap();
    let ft_lines: Vec<&str> = ft.lines().collect();
    assert!(ft_lines[0].starts_with("#OTU ID\t"), "merged feature table header malformed");
    assert_eq!(
        ft_lines[0].split('\t').count(), 3,
        "expected 2 sample columns when exporting a pooled-samples directory"
    );

    // ASV IDs in feature table and rep seqs must agree
    let ft_ids: HashSet<String> = ft_lines[1..].iter()
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    let rs_path = export_out.path().join("merged_rep_seqs.fasta");
    assert!(rs_path.exists(), "merged_rep_seqs.fasta not created");
    let rs_ids: HashSet<String> = fs::read_to_string(&rs_path).unwrap()
        .lines()
        .filter(|l| l.starts_with('>'))
        .map(|l| l[1..].split_whitespace().next().unwrap().to_string())
        .collect();
    assert_eq!(ft_ids, rs_ids, "feature table IDs must match rep_seqs IDs in pooled export");
}
