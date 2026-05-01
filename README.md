# savont - Amplicon Sequence Variants (ASVs) and taxonomic profiling for long read amplicons

**Savont** generates [**Amplicon Sequence Variants (ASVs)**](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) at **single-nucleotide resolution** from long-read amplicon sequencing data such as

- Oxford Nanopore (ONT) R10.4 sequencing (preferably with SUP basecalling)
- PacBio HiFi sequencing

Savont differs from mapping-based approaches (e.g. Emu or ONT's epi2me workflow). Savont instead follows the Reads -> ASV -> Classification paradigm (just like DADA2, but for noisier long reads).

<p align="center">
    <img alt="github-diagram" src="https://github.com/user-attachments/assets/bf200337-2208-4336-9c5c-6a164e57b3ec" width="75%" />
</p>

## Why savont?

- Savont can resolve ASVs that differ by a single nucleotide. This differs from existing long-read workflows that do fuzzy OTU-like clustering.
- For ONT amplicons, savont requires ~10x less depth to generate ASVs compared to DADA2 / UNOISE.
- Savont also has built-in support for full taxonomic profiling (fastq -> abundance table) for several rRNA databases. 

> [!NOTE]
> Savont is optimized for long reads with >98% accuracy. ONT's R10.4 reads with SUP basecalling or PacBio HiFi are preferred.
> For lower quality reads (e.g. R9.4 ONT data) savont may **not** be useful.

## Preliminary results

We have compiled some [very preliminary results available here](https://github.com/bluenote-1577/savont/wiki/Savont-Preliminary-Results).

<p align="center">
<img width="60%" alt="detection_probability_vs_depth_nm_0" src="https://github.com/user-attachments/assets/19dbacad-2856-4888-a75e-0f1406f73265" />
</p>
<p align="center">
Savont is an order of magnitude (or two) more sensitive for ASV retrieval compared to denoising methods for ONT R10.4 sup reads. Yet, it can retrieve most ASVs (i.e., exact, multi-copy 16s sequences) in this dataset (Zymo Microbial Community Standard). 
</p>


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
[![Anaconda-Server Badge](https://anaconda.org/bioconda/savont/badges/version.svg)](https://anaconda.org/bioconda/savont)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/savont/badges/latest_release_date.svg)](https://anaconda.org/bioconda/savont)

```sh
mamba install -c bioconda savont
# or use conda instead of mamba
```

## Quick start

### Step 1: Generate ASVs from reads

> [!NOTE]
> Savont filters reads based on length and quality. However, savont does not do adapter/primer trimming. Please trim your reads with e.g. [cutadapt](https://cutadapt.readthedocs.io/en/stable/) first. 
```sh
# Full-length 16S rRNA reads -> ASVs
savont asv 16s_full-length.fastq.gz -o savont-out -t 20 

# Full bacterial rRNA operon amplicons -> ASVs
savont asv operon_reads.fastq.gz -o savont-out -t 20 --rrna-operon

# For single-stranded protocols
savont asv 16s_single_strand.fq --single-strand -o savont-out -t 20

# Other types of amplicons with known length
savont asv amplicons.fastq.gz -o savont-out -t 20 --min-read-length 1600 --max-read-length 2100 

# Resulting ASVs
ls savont-out/final_asvs.fasta
```

### Importing ASVs into QIIME2

If you prefer to use QIIME: `savont asv` writes a `feature-table.tsv` alongside `final_asvs.fasta`. Use these two files to create QIIME2 artifacts:

```sh
# Convert TSV feature table to BIOM format (requires the biom-format package)
biom convert -i savont-out/feature-table.tsv -o feature-table.biom \
    --table-type="OTU table" --to-hdf5

# Import feature table
qiime tools import \
    --type 'FeatureTable[Frequency]' \
    --input-path feature-table.biom \
    --output-path feature-table.qza

# Import representative sequences
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path savont-out/final_asvs.fasta \
    --output-path rep-seqs.qza
```

`feature-table.qza` and `rep-seqs.qza` are the standard inputs for QIIME2 diversity, taxonomy, and differential abundance plugins.

## Taxonomic profiling against a reference database

Savont can also classify ASVs and generate a taxonomic profile with abundances. Savont supports two classification approaches:

- **`savont classify`** — minimap2 alignment against database with species- and genus-level output (better for species level)
- **`savont sintax`** — SINTAX k-mer bootstrap; genus-level only (better for unknown taxa)

### Step 2: Download a reference database

```sh
# Download one or more databases (emu-1, silva-138.2, greengenes2-2024.09)
savont download --location /path/databases --dbs emu-1
savont download --location /path/databases --dbs silva-138.2
savont download --location /path/databases --dbs greengenes2-2024.09

# Or download multiple at once
savont download --location /path/databases --dbs emu-1 silva-138.2 
```

### Step 3a: Classify ASVs with alignment (`savont classify`)

```sh
# Classify using any downloaded database (type is auto-detected)
savont classify -i savont-out -d databases/emu-1 -t 20
savont classify -i savont-out -d databases/silva-138.2 -t 20

# Write to a separate output directory
savont classify -i savont-out -d databases/emu-1 -o classification-out -t 20

# Adjust identity thresholds
savont classify -i savont-out -d databases/emu-1 \
    --species-threshold 99.9 --genus-threshold 90.0
```

### Step 3b: Classify ASVs with SINTAX (`savont sintax`)

`savont sintax` uses 12-mer bootstrap resampling (100 iterations by default) to assign genus-level taxonomy with per-rank confidence scores.

```sh
savont sintax -i savont-out -d databases/emu-1 -t 20

# Adjust bootstrap threshold (default: 0.80)
savont sintax -i savont-out -d databases/silva-138.2 --min-bootstrap 0.70
```

## Database Information

### EMU (`emu-1`)

From [Emu](https://github.com/treangenlab/emu) by Curry et al. (2022, Nature Methods). Curated 16S rRNA database with focused species-level classifications.

### SILVA (`silva-138.2`) — SSU Ref NR99 v138.2

More comprehensive than EMU, especially for understudied taxa. Species-level classifications are often split across multiple distinct strains.

### GreenGenes2 (`greengenes2-2024.09`)

GreenGenes2 2024.09 species-level trainset (DADA2 format). Unannotated ranks are reported as `Greengenes_unannotated`.

## Output

### ASV Clustering Output

The `savont asv` command produces:

1. **final_asvs.fasta** - Final ASV sequences (high-quality, chimera-filtered)
2. **feature-table.tsv** - QIIME2-compatible feature table (ASV × sample read counts)
3. **final_clusters.tsv** - Cluster assignments mapping reads to ASVs
4. **temp/** - Directory containing intermediate files

### Classification Output (`savont classify`)

The `savont classify` command produces three output files similar to Emu:

#### 1. species_abundance.tsv / genus_abundance.tsv

Species-level or genus-level taxonomic abundance table:

```
abundance       species         genus   family  order   class   phylum  clade   superkingdom
0.45123         Escherichia_coli        Escherichia     Enterobacteriaceae      ...
0.23456         Staphylococcus_aureus   Staphylococcus  Staphylococcaceae       ...
```

- `abundance` - Relative abundance estimated by EM algorithm
- `tax_id` - Taxonomic identifier
- Full taxonomic lineage from species to superkingdom

#### 2. asv_mappings.tsv

Individual ASV mapping details:

```
asv_header      depth     alignment_identity      number_mismatches       tax_id  species genus   reference
final_consensus_0_depth_5936    5936    99.67   5       29466   Veillonella parvula     Veillonella     29466:emu_db:36875
final_consensus_1_depth_3081    3081    99.27   11      29466   Veillonella parvula     Veillonella     29466:emu_db:36873
final_consensus_2_depth_2927    2927    99.40   9       29466   Veillonella parvula     Veillonella     29466:emu_db:36869
```

The best mapping references and their corresponding species/genus are denoted.

### Classification Output (`savont sintax`)

`savont sintax` produces **genus_abundance.tsv** (same format as above) and **asv_mappings.tsv** with bootstrap confidence scores:

```
asv_header      depth   genus_bootstrap family_bootstrap  genus         family            order   ...
final_consensus_0_depth_5936    5936    0.980   0.990   Veillonella   Veillonellaceae   ...
```

Ranks below `--min-bootstrap` are reported as `UNCLASSIFIED`.

## Algorithm Overview

### ASV Generation Pipeline

Savont clusters reads using novel polymorphic marker (SNPmer) algorithms, polishes to get consensus ASVs, removes chimeras, and refines depths with an EM algorithm over read-level alignments.

### Taxonomic Classification

**`savont classify`**: maps ASVs to the database with minimap2, resolves multi-mappers with an EM algorithm, then filters by identity for species/genus-level calls.

**`savont sintax`**: re-implements [Edgar's sintax algorithm](https://www.biorxiv.org/content/10.1101/074161v1). Subsamples 32 canonical 12-mers per ASV per iteration (100 iterations), finds the best-matching database entry for each subsample, and reports the fraction of iterations supporting each rank as the bootstrap confidence.



### CHANGELOG

See [the changelog.](CHANGELOG.md)

## Citation

FORTHCOMING WORK.

**If you use any provided database, cite:**

1. Quast, Christian, et al. "The SILVA ribosomal RNA gene database project: improved data processing and web-based tools." Nucleic acids research 41.D1 (2012): D590-D596.

2. Curry, Kristen D., et al. "Emu: species-level microbial community profiling of full-length 16S rRNA Oxford Nanopore sequencing data." Nature methods 19.7 (2022): 845-853.

3. McDonald, D., Jiang, Y., Balaban, M. et al. Greengenes2 unifies microbial data in a single reference tree. Nat Biotechnol 42, 715–718 (2024). 

## License

MIT
