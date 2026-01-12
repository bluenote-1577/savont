# Changelog

All notable changes to savont will be documented in this file.

## [0.3.0] - 2025-1-10

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
