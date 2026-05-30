# savont - Amplicon Sequence Variants (ASVs) and taxonomic profiling for long read amplicons


**Savont** generates [**Amplicon Sequence Variants (ASVs)**](https://en.wikipedia.org/wiki/Amplicon_sequence_variant) at **single-nucleotide resolution** from long-read amplicon sequencing data such as

- Oxford Nanopore (ONT) R10.4 sequencing (preferably with SUP basecalling)
- PacBio HiFi sequencing

Savont differs from mapping-based approaches (e.g. Emu or ONT's epi2me workflow). Savont instead follows the Reads -> ASV -> Classification paradigm (just like DADA2, but for noisier long reads).

## Why savont?

- Savont can resolve ASVs that differ by a single nucleotide, **even for nanopore reads.** This differs from existing nanopore workflows for OTU-like clustering.
- For ONT amplicons, savont requires ~10x less depth to generate ASVs compared to DADA2 / UNOISE.
- Savont also has built-in support for full taxonomic profiling (fastq -> abundance table -> QIIME) for several rRNA databases. 

> [!NOTE]
> Savont is optimized for long reads with >98% accuracy. ONT's R10.4 reads with SUP basecalling or PacBio HiFi are preferred.
> For lower quality reads (e.g. R9.4 ONT data) savont may **not** be useful.

<p align="center">
    <img width="90%" alt="github-diagram" src="https://github.com/user-attachments/assets/c0d9e356-ee1d-4d60-a217-c050e5abd0dc" />
</p>

## Results

> [!tip]
> [The savont preprint](https://www.biorxiv.org/content/10.64898/2026.05.26.727271v1) is now out! Please see the preprint for results and methodological details. 

<p align="center">
<img width="90%" alt="image" src="https://github.com/user-attachments/assets/e4da4032-e080-48a5-b04a-e8bbd2ad1976" />
    <p align="center"><i>ONT and HiFi results on 8-species community for full-length 16S sequences.</i></p>
</p>

<p align="center">
<img width="90%" alt="image" src="https://github.com/user-attachments/assets/0ce3f9d9-3040-4e87-aa35-f798767acdc7" />
    <p align="center"><i>Diversity and ASV metrics on complex full-length ONT 16S amplicons.</i> </p>
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

## Quick start for trimmed full-length 16S amplicons

```sh
# generate 16S ASVs for ONT / HiFi after cutadapt
savont asv 16s_full-length_trimmed.fastq.gz -o savont-out -t 20 (optional: --hifi)

# download Greengenes, EMU, or silva databases
savont download --location /path/databases --dbs greengenes2-2024.09 emu-1 silva-138.2

# classify against any database
savont classify -i savont-out -d /path/databases/emu-1 -t 20

# merge multiple savont outputs --> feeds into QIIME2 (optional)
savont merge -i savont-out1 savont-out2 -o merged_output
```

# Full subcommand guide

## Generate ASVs from reads (`savont asv`)

> [!NOTE]
> Savont filters reads based on length and quality. However, savont does not do adapter/primer trimming. Please trim your reads with e.g. [cutadapt](https://cutadapt.readthedocs.io/en/stable/) first.

Default savont parameters assume ONT / HiFi full-length 16S rRNA sequencing with trimmed reads from cutadapt. See below for parameter choices if this is not the case. 
 
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

### ASV output

The `savont asv` command produces:

1. **final_asvs.fasta** - Final ASV sequences (high-quality, chimera-filtered)
2. **feature-table.tsv** - QIIME2-compatible feature table (ASV × sample read counts)
3. **final_clusters.tsv** - Cluster assignments mapping reads to ASVs
4. **temp/** - Directory containing intermediate files


## Taxonomic profiling against a reference database

Savont can also classify ASVs and generate a taxonomic profile with abundances. Savont supports two classification approaches:

- **`savont classify`** — minimap2 alignment against database with species- and genus-level output (better for species level). Uses identity thresholds for taxonomic assignment (thresholds from [Yarza et al.](https://www.nature.com/articles/nrmicro3330))
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


### Classification Output (`savont classify`)

The `savont classify` command produces three output files:

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
asv_header      depth   alignment_identity  number_mismatches   tax_id  species             genus       family              order               class                   phylum          clade   superkingdom    reference
final_consensus_0_depth_5936    5936    99.67   5   29466   Veillonella parvula Veillonella Veillonellaceae Veillonellales  Negativicutes   Bacillota           Bacteria    29466:emu_db:36875
final_consensus_1_depth_3081    3081    99.27   11  29466   Veillonella parvula Veillonella Veillonellaceae Veillonellales  Negativicutes   Bacillota           Bacteria    29466:emu_db:36873
```

All taxonomic ranks from species to superkingdom are included. Unclassified ASVs have `UNCLASSIFIED` in every rank column.

### Classification Output (`savont sintax`)

`savont sintax` produces **genus_abundance.tsv** (same format as above) and **asv_mappings.tsv** with bootstrap confidence scores:

```
asv_header      depth   species_bootstrap   genus_bootstrap family_bootstrap    order_bootstrap class_bootstrap phylum_bootstrap    superkingdom_bootstrap  species         genus       family          ...
final_consensus_0_depth_5936    5936    0.000   0.980   0.990   0.995   0.999   1.000   1.000   UNCLASSIFIED    Veillonella Veillonellaceae  ...
```

Ranks below `--min-bootstrap` are reported as `UNCLASSIFIED`.


## Merge multiple samples (`savont merge`)

When you have run `savont asv` (and optionally `classify`/`sintax`) on multiple samples separately, `savont merge` combines them into a single set of outputs. 

ASVs across samples are matched by exact sequence hash. By default, savont also **fuzzy-merges** ASVs with identical sequences but trimmed to slightly different lengths across runs. 

```sh
# Merge two samples
savont merge -i sample1-out sample2-out -o merged-out

# Merge many samples and give them meaningful names
savont merge \
    -i run1/savont-out run2/savont-out run3/savont-out \
    -o merged-out \
    --relabel SampleA SampleB SampleC
```

> [!NOTE]
> Make sure to use `--relabel` if you have duplicate sample names. 

#### Merge outputs

1. **merged_feature_table.tsv** — QIIME2-compatible feature table; rows are hash-keyed ASVs, columns are samples
2. **merged_rep_seqs.fasta** — representative sequences for all merged ASVs
3. **merged_asv_taxonomy.tsv** — ASV-level taxonomy (Feature ID → lineage); **only written if `savont classify` or `savont sintax` was run on any input directory**
4. **merged_taxon_counts.tsv** — human-readable taxon count table; rows are lineage strings (`Bacteria;Firmicutes;...`), columns are sample counts; useful for quick inspection without QIIME2

#### Importing merged outputs into QIIME2 + stacked bar plot

```sh
# Feature table
biom convert -i merged-out/merged_feature_table.tsv \
    -o feature-table.biom --table-type='OTU table' --to-hdf5
qiime tools import --type 'FeatureTable[Frequency]' \
    --input-path feature-table.biom --output-path feature-table.qza

# Representative sequences
qiime tools import --type 'FeatureData[Sequence]' \
    --input-path merged-out/merged_rep_seqs.fasta --output-path rep-seqs.qza

# Taxonomy (if classify/sintax was run)
qiime tools import --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path merged-out/merged_asv_taxonomy.tsv --output-path taxonomy.qza

qiime taxa barplot --i-table feature-table.qza --i-taxonomy taxonomy.qza \
    --o-visualization taxa-bar-plots.qzv
```


## Database Information

### EMU (`emu-1`)

From [Emu](https://github.com/treangenlab/emu) by Curry et al. (2022, Nature Methods). Curated 16S rRNA database with focused species-level classifications. Less breadth than SILVA or Greengenes2 for complex microbiomes. 

### SILVA (`silva-138.2`) — SSU Ref NR99 v138.2

More comprehensive than EMU, especially for understudied taxa. However, species-level classifications are often split across multiple distinct strains.

### GreenGenes2 (`greengenes2-2024.09`)

GreenGenes2 2024.09 species-level trainset (DADA2 format). Unannotated ranks are reported as `Greengenes_unannotated`.

### CHANGELOG

See [the changelog.](CHANGELOG.md)

## Citation

Jim Shaw, Marie Riisgaard-Jensen, Kasper Skytte Andersen, Rasmus Kirkegaard, Morten Kam Dahl Dueholm, Heng Li. bioRxiv 2026.05.26.727271; doi: https://doi.org/10.64898/2026.05.26.727271

**If you use any provided database, cite:**

1. Quast, Christian, et al. "The SILVA ribosomal RNA gene database project: improved data processing and web-based tools." Nucleic acids research 41.D1 (2012): D590-D596.

2. Curry, Kristen D., et al. "Emu: species-level microbial community profiling of full-length 16S rRNA Oxford Nanopore sequencing data." Nature methods 19.7 (2022): 845-853.

3. McDonald, D., Jiang, Y., Balaban, M. et al. Greengenes2 unifies microbial data in a single reference tree. Nat Biotechnol 42, 715–718 (2024).

**If you use the SINTAX algorithm please cite**: Edgar, Robert C. "SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences." biorxiv (2016): 074161.

## License

MIT
