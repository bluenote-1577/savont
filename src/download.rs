use std::process::Command;
use crate::cli;

/// Download EMU database
fn download_emu_database(location: &str) -> Result<(), String> {
    log::info!("Downloading EMU database to {}", location);

    // Create location directory if it doesn't exist
    std::fs::create_dir_all(location)
        .map_err(|e| format!("Failed to create directory {}: {}", location, e))?;

    log::info!("Downloading emu_default.tar.gz...");

    // Download the database
    let download_status = Command::new("wget")
        .arg("--content-disposition")
        .arg("https://osf.io/8qcwd/download")
        .arg("-P")
        .arg(location)
        .status()
        .map_err(|e| format!("Failed to execute wget: {}. Is wget installed?", e))?;

    if !download_status.success() {
        return Err("Failed to download EMU database".to_string());
    }

    log::info!("Extracting emu_default.tar.gz...");

    // Extract the tar.gz file
    let tar_path = format!("{}/emu_default.tar.gz", location);
    let extract_status = Command::new("tar")
        .arg("-xzf")
        .arg(&tar_path)
        .arg("-C")
        .arg(location)
        .status()
        .map_err(|e| format!("Failed to execute tar: {}. Is tar installed?", e))?;

    if !extract_status.success() {
        return Err("Failed to extract EMU database".to_string());
    }

    log::info!("EMU database successfully downloaded and extracted to {}", location);

    Ok(())
}

/// Download SILVA database
fn download_silva_database(location: &str) -> Result<(), String> {
    log::info!("Downloading SILVA database to {}", location);

    // Create silva_db directory
    let silva_dir = format!("{}/silva_db", location);
    std::fs::create_dir_all(&silva_dir)
        .map_err(|e| format!("Failed to create directory {}: {}", silva_dir, e))?;

    log::info!("Downloading SILVA FASTA file...");

    // Download FASTA file
    let fasta_url = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz";
    let fasta_status = Command::new("wget")
        .arg(fasta_url)
        .arg("-P")
        .arg(&silva_dir)
        .status()
        .map_err(|e| format!("Failed to execute wget: {}. Is wget installed?", e))?;

    if !fasta_status.success() {
        return Err("Failed to download SILVA FASTA file".to_string());
    }

    log::info!("Downloading SILVA taxonomy file...");

    // Download taxonomy file
    let tax_url = "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.2.txt.gz";
    let tax_status = Command::new("wget")
        .arg(tax_url)
        .arg("-P")
        .arg(&silva_dir)
        .status()
        .map_err(|e| format!("Failed to execute wget: {}", e))?;

    if !tax_status.success() {
        return Err("Failed to download SILVA taxonomy file".to_string());
    }

    Command::new("gzip")
        .arg("-d")
        .arg(format!("{}/taxmap_slv_ssu_ref_nr_138.2.txt.gz", silva_dir))
        .status()
        .map_err(|e| format!("Failed to execute gzip: {}", e))?;

    log::info!("SILVA database successfully downloaded to {}", silva_dir);

    Ok(())
}

/// Main download function called from main.rs
pub fn download(args: &cli::DownloadArgs) {
    if args.db_type.emu_db {
        match download_emu_database(&args.location) {
            Ok(_) => {
                log::info!("EMU database download complete!");
                log::info!("You can now use it with: savont classify --emu-db {}/emu_default", args.location);
            }
            Err(e) => {
                log::error!("Failed to download EMU database: {}", e);
                std::process::exit(1);
            }
        }
    } 
    if args.db_type.silva_db {
        match download_silva_database(&args.location) {
            Ok(_) => {
                log::info!("SILVA database download complete!");
                log::info!("You can now use it with: savont classify --silva-db {}/silva_db", args.location);
            }
            Err(e) => {
                log::error!("Failed to download SILVA database: {}", e);
                std::process::exit(1);
            }
        }
    }
}
