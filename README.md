# savont - Amplicon Sequencing Variants (ASVs) and taxonomic profiling for long-read amplicon sequencing 

**Savont** generates **Amplicon Sequence Variants (ASVs)** from Oxford Nanopore (ONT) R10.4 or PacBio HiFi amplicon sequencing data. 

Unlike mapping-based approaches (e.g. Emu or ONT's epi2me workflow), savont follows the Reads -> ASV -> Classification paradigm (just like DADA2). This can lead to more confident species classifications and exact ASV mapping information that is lost from read-level analysis. 

## What can savont do?

Savont is designed for high-accuracy long-read amplicon sequencing (>98% accuracy). It can:

* Generate high-resolution ASVs from ONT R10.4 or PacBio HiFi amplicon reads (e.g., full-length 16S or rRNA operon)
* Classify ASVs against EMU or SILVA reference databases
* Output species-level and genus-level taxonomic abundance tables
* Provide detailed ASV mapping information for quality control

> [!NOTE]
> Savont is optimized for long reads with >98% accuracy. R10.4 SUP ONT reads or HiFi are preferred.
> For lower quality reads (e.g. R9.4 ONT data or HAC/FAST base-called data) savont may **not** be useful.

## Install via conda or build from source

### Option 1: Build from source

Requirements:
1. [rust](https://www.rust-lang.org/tools/install) (tested for > v1.88) programming language
2. Standard linux toolchain (tar, gzip, wget, C++, gcc)
3. cmake

```sh
git clone https://github.com/bluenote-1577/savont
cd savont

# Build and install
cargo install --path .
savont --help
```

### Option 2: Conda

```sh
mamba install -c bioconda savont
# or use conda instead of mamba
```

## Quick start

### Step 1: Download a reference database

```sh
# Download EMU database
savont download --location databases --emu-db

# Or download SILVA database
savont download --location databases --silva-db

```

### Step 2: Generate ASVs from reads

```sh
# Cluster full-length 16S rRNA reads into ASVs
savont asv reads.fastq.gz -o savont-out -t 20 --min-read-length 1100 --max-read-length 2000 

# Cluster full bacterial rRNA operon amplicons into ASVs
savont asv reads.fastq.gz -o savont-out -t 20 --min-read-length 4000 --max-read-length 6500 

# For single-stranded protocols
savont asv --single-strand -o savont-out -t 20

ls savont-out/final_asvs.fasta
```

### (Optional) Step 3: Classify prokaryotic rRNA ASVs

```sh
# Classify using EMU database
savont classify -i savont-out -o classification-out --emu-db databases/emu_default -t 20

# Classify using SILVA database
savont classify -i savont-out -o classification-out --silva-db databases/silva_db -t 20

# Adjust identity thresholds
savont classify -i savont-out --emu-db databases/emu_default \
    --species-threshold 99.9 --genus-threshold 90.0
```

## Output

### ASV Clustering Output

The `savont asv` command produces:

1. **final_asvs.fasta** - Final ASV sequences (high-quality, chimera-filtered)
2. **final_clusters.tsv** - Cluster assignments mapping reads to ASVs
3. **temp/** - Directory containing intermediate files:

### Classification Output

The `savont classify` command produces three output files similar to Emu:

#### 1. species_abundance.tsv

Species-level taxonomic abundance table:

```
abundance       tax_id  species         genus   family  order   class   phylum  clade   superkingdom
0.45123         562     Escherichia_coli        Escherichia     Enterobacteriaceae      ...
0.23456         1280    Staphylococcus_aureus   Staphylococcus  Staphylococcaceae       ...
```

- `abundance` - Relative abundance estimated by EM algorithm
- `tax_id` - Taxonomic identifier
- Full taxonomic lineage from species to superkingdom

#### 2. genus_abundance.tsv

Genus-level collapsed abundance table (species aggregated to genus):

```
abundance       genus   family  order   class   phylum  clade   superkingdom
0.50123         Escherichia     Enterobacteriaceae      Enterobacterales        ...
0.30456         Staphylococcus  Staphylococcaceae       Bacillales              ...
```

- Aggregates all species within each genus
- Useful for genus-level community analysis

#### 3. asv_mappings.tsv

Individual ASV mapping details:

```
asv_header      depth     alignment_identity      number_mismatches       tax_id  species genus   reference
final_consensus_0_depth_5936    5936    99.67   5       29466   Veillonella parvula     Veillonella     29466:emu_db:36875
final_consensus_1_depth_3081    3081    99.27   11      29466   Veillonella parvula     Veillonella     29466:emu_db:36873
final_consensus_2_depth_2927    2927    99.40   9       29466   Veillonella parvula     Veillonella     29466:emu_db:36869
```

The best mapping references and their corresponding species/genus are denoted. 

## Algorithm Overview

### ASV Generation Pipeline

Savont does the following: cluster reads --> polish and get consensus ASVs --> remove chimeric ASVs. 

Savont uses novel algorithms for clustering with polymorphic markers. We also use some heuristics and statistics to deal with errors, inspired by existing ASV approaches (but adapted to long reads). 

### Taxonomic Classification

Savont has built-in taxonomic profiling for prokaryotic 16S sequences. The `classify` command does the following:

1. **Map ASVs to database using minimap2**
2. **Deal with ambiguous ASV alignments using EM algorithm**
3. **Filter low-identity mappings for species and genus-level classification**

## Database Information

### EMU Database

From [Emu](https://github.com/treangenlab/emu) by Curry et al. (2022, Nature Methods)

- Has more "focused" species classifications
- Lacks breadth of SILVA

### SILVA Database (v138.2) - Non-redundant 99%

- More comprehensive than EMU, especially for understudied species
- Species-level classifications may be split over multiple distinct species

### Notes about Quality Control

1. Check `asv_mappings.tsv` for ASV depth distribution
2. Low-depth ASVs (<20 reads) may be artifacts or rare taxa
3. Examine unmapped ASVs in `asv_mappings.tsv` or the log.

### CHANGELOG

See [the changelog.](CHANGELOG.md)

## Citation

FORTHCOMING. 

If you use SILVA or EMU databases, cite:

1. Quast, Christian, et al. "The SILVA ribosomal RNA gene database project: improved data processing and web-based tools." Nucleic acids research 41.D1 (2012): D590-D596.

2. Curry, Kristen D., et al. "Emu: species-level microbial community profiling of full-length 16S rRNA Oxford Nanopore sequencing data." Nature methods 19.7 (2022): 845-853.


## License

MIT
