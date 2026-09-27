use crate::cli;
use crate::databases;
use std::path::Path;

pub fn download(args: &cli::DownloadArgs) {
    for kw in &args.dbs {
        // The keyword was already validated by clap's PossibleValuesParser.
        let def = databases::find(kw)
            .unwrap_or_else(|| panic!("BUG: unknown keyword '{}' passed CLI validation", kw));

        let dest = Path::new(&args.location).join(kw);
        if let Err(e) = std::fs::create_dir_all(&dest) {
            log::error!("Failed to create directory {}: {}", dest.display(), e);
            std::process::exit(1);
        }

        log::info!(
            "Downloading '{}' ({}) to {} ...",
            kw,
            def.description,
            dest.display()
        );

        match (def.download)(&dest) {
            Ok(()) => {
                databases::write_marker(&dest, kw).ok();
                log::info!("'{}' downloaded successfully.", kw);
                log::info!("Use with: savont classify -d {}", dest.display());
            }
            Err(e) => {
                log::error!("Failed to download '{}': {}", kw, e);
                std::process::exit(1);
            }
        }
    }
}
