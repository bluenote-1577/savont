# Changelog

All notable changes to savont will be documented in this file.

## [0.5.1] - 2026-6-14

- Added a way to autodetect samples with low polymorphism. In this case, we don't proceed with SNPmer clustering and just proceed with k-mer clusters + polishing.
- Added `-m` and `-M` for min read length and max read length CLI arguments. 

## [0.5.0] - 2026-05-11

### Changed

- **`savont merge`** command added. Allows merging of multiple savont runs with outputs that are easily useable by QIIME. Merges ASVs across multiple samples with a simple fuzzy clustering algorithm. 
- **Fixed ASV sequence boundary issues**: ASV ends were getting cut off during consensus generation. This is now fixed -- ASV boundaries should be more consistent across samples.
- **Fixed `savont classify` thresholds**: Using taxonomic thresholds for sequence identity for `savont classify` from Yarza et al. by default. 

## [0.4.0] - 2026-04-28

### Changed

- **Classify / Download command** usage changed. Databases are refactored, download logic are changed, and greengenes2 is added. 
- Tax IDs are removed from `species_abundance.tsv`. 

### Added

- **`savont sintax` command**: new k-mer bootstrap classification subcommand implementing the SINTAX algorithm for genus-level taxonomic classification. Uses 12-mers with 100 bootstrap iterations (32 k-mers sampled per iteration) to produce per-rank confidence scores. Key options: `--min-bootstrap` (default 0.8), `--n-iter` (default 100).
- **GreenGenes2 2024.09 database**: support for the GreenGenes2 species-level trainset via `savont download --dbs greengenes2-2024.09`. Ranks missing annotation are reported as `Greengenes_unannotated`.
- **Multiple database downloads**: `savont download --dbs` now accepts multiple database keywords in one invocation.
- **Auto-detection of database type**: `savont classify` and `savont sintax` auto-detect the database type from a `.savont_db` marker file or the directory name, removing the need to specify the database format manually.
- **`--detailed-unclassified` flag** for `savont classify` and `savont sintax`: outputs `UNCLASSIFIED-(asv_header)` instead of `UNCLASSIFIED` for ranks below the confidence threshold. Previously, `...-(asv_header)` was the standard behavior. 
- **`feature-table.tsv` output** from `savont asv`: QIIME2-compatible feature table (TSV format) written alongside `final_asvs.fasta`. Import with `biom convert` and `qiime tools import --type 'FeatureTable[Frequency]'`.


## [0.3.2] - 2025-1-10

### Added

- Detects cutadapt "rc" flag in fastq file and appropriately handles reverse complements.
- Added `--use-blockmers` flag to enable experimental blockmer-based polymorphic marker clustering
- Slightly tweaked consensus generation algorithm. Not too much of a change
- Fixed bugs with detecting chimeras. Should yield much better results for super high depth stuff. 
- Much better logging.  
- Added presets for operon processing. 

## [0.2.0] - 2025-12-29

### Added
- Added `--min-read-length` and `--max-read-length` parameters for flexible read length filtering
- Added `--posterior-threshold-ln` parameter for consensus quality control (default: 30.0)
- Added `--max-iterations-recluster` parameter to control reclustering iterations (default: 10)

### Changed
- **Significantly improved chimera detection algorithm**: Enhanced depth-aware, alignment-based detection
- **Major performance improvements**: Optimized runtime across all clustering steps

### Removed
- Removed `--not-full-16s` flag 


## [0.1.0] - 2024-12-XX

### Added
- Initial release of savont
