//! Checks that every database download URL is still reachable, without downloading.
//! Needs network access and `wget` (the same tool `savont download` uses).
//! Run on its own with `cargo test --test database_urls`.

use savont::databases;
use std::process::Command;

#[test]
fn database_download_urls_are_reachable() {
    let mut broken = Vec::new();
    for def in databases::ALL {
        assert!(!def.urls.is_empty(), "{} lists no download URLs", def.keyword);
        for url in def.urls {
            let status = Command::new("wget")
                .args(["--spider", "-q", "--tries=2", "--timeout=30", url])
                .status()
                .expect("failed to run wget");
            if !status.success() {
                broken.push(format!("{} -> {} (wget exit {:?})", def.keyword, url, status.code()));
            }
        }
    }
    assert!(broken.is_empty(), "unreachable database URLs:\n{}", broken.join("\n"));
}
