# Full-length 16S rRNA (Oxford Nanopore) bioinformatic pipeline

Reproducible bioinformatic workflow for taxonomic profiling of a bacterial microbiota
from full-length 16S rRNA (V1 to V9) Oxford Nanopore sequencing. It takes basecalled,
quality-filtered reads and produces combined abundance tables at six taxonomic ranks.

The downstream community-level statistical analysis (contaminant removal, diversity,
PERMANOVA and related tests) is documented in a separate repository/pipeline.

Steps:

1. Read subsampling (depth standardization) with seqkit
2. Chimera detection and removal with VSEARCH (UCHIME)
3. Taxonomic classification with EMU (minimap2 plus expectation-maximization)
4. Combined abundance tables at six taxonomic ranks (EMU combine-outputs)

```
basecalled FASTQ (per barcode)
      |  quality and length filtering (Dorado: 1,300 to 1,800 bp, Q >= 12)
      v
[1] subsample to <= 50,000 reads/sample     (seqkit, fixed seed)          -> chimera_removal.sh
[2] chimera removal                         (VSEARCH: cluster 97%, uchime_denovo) -> chimera_removal.sh
      v
[3] taxonomic classification                (EMU: map-ont, default 16S database)  -> emu_classification.sh
[4] combined abundance tables               (phylum, class, order, family, genus, species) -> emu_classification.sh
      v
    six .tsv tables  ->  statistical analysis (separate pipeline)
```

---

## 1. Study data

- Assay: full-length 16S rRNA gene (V1 to V9), Oxford Nanopore Technologies.
- Instrument and kit: MinION Mk1D, Rapid Sequencing DNA 16S Barcoding Kit V14
  (SQK-16S114.24), FLO-MIN114 (R10.4.1) flow cell, up to 24 barcoded samples per run.
- Basecalling: Dorado (sup model), demultiplexing and adapter/barcode trimming,
  reads retained at 1,300 to 1,800 bp and Q >= 12.
- Input to this pipeline: one FASTQ (or a folder of FASTQ files) per barcode.

---

## 2. Requirements

| Software     | Version   | Purpose                              |
|--------------|-----------|--------------------------------------|
| Linux / WSL2 | -         | Runtime environment                  |
| conda (Miniforge) | -    | Environment management               |
| seqkit       | 2.13.0    | Read subsampling and FASTQ handling  |
| VSEARCH      | 2.32.0    | Chimera detection and removal        |
| EMU          | 3.6.2     | Taxonomic classification             |
| minimap2     | >= 2.22   | Long-read alignment (used by EMU)    |
| osfclient    | -         | Download of the EMU reference DB     |

EMU reference database: default build combining rrnDB and NCBI 16S RefSeq
(49,301 sequences from 17,555 bacterial and archaeal species).

---

## 3. Installation

```bash
# Create the environment
mamba create -n emu -c conda-forge -c bioconda \
    emu "minimap2>=2.22" vsearch seqkit osfclient -y
conda activate emu

# Download the EMU default database
mkdir -p ~/emu_db && export EMU_DATABASE_DIR=~/emu_db
cd ~/emu_db
osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar
tar -xvf emu.tar          # yields species_taxid.fasta and taxonomy.tsv
echo 'export EMU_DATABASE_DIR=~/emu_db' >> ~/.bashrc
```

---

## 4. Repository structure

```
.
├── README.md
├── chimera_removal.sh      # Step 1 and 2 (subsampling + chimera removal)
├── emu_classification.sh   # Step 3 and 4 (EMU + combine-outputs)
├── data/                   # input FASTQ, per barcode  [not tracked]
└── results/                # pipeline outputs           [not tracked]
```

Both scripts read their paths from environment variables or arguments, so no absolute
path is hard-coded. The examples below use generic placeholders such as
`/path/to/data/barcodes`.

---

## 5. Usage

### Step 1 and 2 - Subsampling and chimera removal (chimera_removal.sh)

Standardizes sequencing depth and then removes chimeras, keeping quality scores intact.

Subsampling: each sample is randomly subsampled to a maximum of 50,000 reads (seqkit,
fixed seed = 100); samples below this threshold are kept in full.

Chimera removal: dereplication of high-error long reads yields almost exclusively
singleton sequences, which precludes abundance-based de novo detection. Reads are
clustered at 97% identity (--cluster_size) to recover abundance structure,
--uchime_denovo (abskew 2.0, minh 0.28) is applied to the cluster centroids, and
chimera calls are propagated to every read within a chimeric cluster. Borderline
sequences are retained. The script validates each sample (exact read balance, preserved
quality scores, no surviving chimeric read).

```bash
INPUT_DIR=/path/to/data/barcodes \
OUTPUT_DIR=/path/to/results/chimera \
MAX_READS=50000 SEED=100 CLUSTER_ID=0.97 ABSKEW=2.0 MINH=0.28 \
THREADS=8 \
bash chimera_removal.sh
```

Outputs: `results/chimera/final/<barcode>.sin_quimeras.fastq.gz` (EMU input),
`results/chimera/resumen/resumen_quimeras.csv` (per-sample summary with the percentage
of chimeric reads), plus logs and the full VSEARCH audit trail.

### Step 3 and 4 - Classification and abundance tables (emu_classification.sh)

Classifies each chimera-free sample with EMU (map-ont preset, --keep-counts) and
combines the per-sample profiles into tables at every taxonomic rank. Set BASE to the
directory holding the run, which must contain a final/ subfolder with the chimera-free
FASTQ files.

```bash
export EMU_DATABASE_DIR=~/emu_db
BASE=/path/to/results bash emu_classification.sh RUN_NAME 8
```

Output, the six taxonomic-level tables (in RUN_NAME/emu/):

```
emu-combined-phylum-counts.tsv
emu-combined-class-counts.tsv
emu-combined-order-counts.tsv
emu-combined-family-counts.tsv
emu-combined-genus-counts.tsv
emu-combined-species-counts.tsv
```

Each is a taxa by samples matrix of estimated read counts (and a matching relative
abundance table without the -counts suffix), with the full taxonomic lineage as leading
columns. These tables are the input for the downstream statistical analysis.

---

## 6. Parameter summary

| Parameter | Value | Step |
|-----------|-------|------|
| Read length filter | 1,300 to 1,800 bp | upstream |
| Quality filter | Q >= 12 | upstream |
| Subsampling maximum | 50,000 reads/sample | 1 |
| Subsampling seed | 100 | 1 |
| Chimera clustering identity | 0.97 | 2 |
| UCHIME abskew | 2.0 | 2 |
| UCHIME minh | 0.28 | 2 |
| EMU preset | map-ont | 3 |
| EMU database | rrnDB + NCBI 16S RefSeq (default) | 3 |
| Taxonomic ranks combined | phylum to species | 4 |

---

## 7. Methodological notes

- Subsampling is random and uniform (fixed seed, reproducible). It standardizes depth
  and reduces runtime without altering community proportions. Raw reads are deposited
  unmodified in a public repository, and subsampling is an analysis step only.
- Clustering before UCHIME is required for error-prone long reads, whose exact
  dereplication would otherwise yield only singletons and prevent detection.
- Chimera calls are made at the cluster level (97% identity) and propagated to member
  reads, stated explicitly for transparency.
- Species-level assignments are presumptive, consistent with the known limits of 16S
  for resolving closely related species even with full-length reads.

---

## 8. Citations

- Rognes T, Flouri T, Nichols B, Quince C, Mahe F (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584.
- Curry KD, et al. (2022). Emu: species-level microbial community profiling of full-length 16S rRNA Oxford Nanopore sequencing data. Nature Methods 19:845-853.
- Li H (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34:3094-3100.
- Shen W, Le S, Li Y, Hu F (2016). SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS ONE 11:e0163962.
- Hakimzadeh A, et al. (2025). Are we throwing away good data? Evaluation of chimera detection algorithms on long-read amplicons. PeerJ 13:e20456.

---

## 9. Reproducibility

All stochastic steps use fixed seeds, and tool versions are pinned in Section 2. Raw
sequencing reads are available under [accession, to be added]. This repository contains
the bioinformatic analysis code and the derived abundance tables. The community-level
statistical analysis is provided as a separate pipeline.
