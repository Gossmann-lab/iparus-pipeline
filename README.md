# IPARUS pipeline

This repository contains the scripts used to build the **IPARUS** great tit (*Parus major*) genomic database, from raw sequencing reads to the BAM files and sample information shown on the website.

**This repository matches IPARUS Database Version 1.0.**

- Website: [IPARUS](https://gossmann-lab.github.io/iparus/)
- Manuscript: TODO add citation

---

## Repository contents

```
iparus-pipeline/
├── README.md                     ← this file
├── scripts/
│   ├── 01_make_bams.sh           ← step 1: raw reads → BAM files
│   └── 02_sex_check.sh           ← step 2: genetic sex of each sample
├── data/
│   ├── samples.txt               ← example sample list (input to steps 1 and 2)
│   └── bam_checksums.txt         ← checksums of the Version 1.0 BAMs
└── example_analysis/             ← example analysis from manuscript 
    ├── README.md                 ← description, tools, steps, outputs
    ├── fst_chr1A.sh              ← FST P. minor vs. P. major
    └── pca_chr1A.sh              ← PCA of 20 individuals
```

### Files in `data/`

Where each file in `data/` comes from:

| file | made by | contents |
|---|---|---|
| `samples.txt` | written by hand | one sample name per line |
| `sex_results.tsv` | `scripts/02_sex_check.sh` | Z/autosome depth ratio and sex of each sample |
| `bam_checksums.txt` | `md5sum` on the Version 1.0 BAM files made by `scripts/01_make_bams.sh` | MD5 checksum of each BAM file on IPARUS |

Large files (raw reads and BAMs) are not stored on GitHub; see [Data storage](#data-storage).

---

## Workflow

Run the scripts in this order, from the repository folder:

| step | script | input | output |
|---|---|---|---|
| 1 | `scripts/01_make_bams.sh` | raw paired-end reads (`SAMPLE_1.fq.gz`, `SAMPLE_2.fq.gz`), reference genome, `data/samples.txt` | `bams/SAMPLE/SAMPLE_sorted_dedup.bam` (+ index, QC reports, metrics) |
| 2 | `scripts/02_sex_check.sh` | BAMs from step 1, `data/samples.txt` | `data/sex_results.tsv` |

```bash
git clone https://github.com/gossmann-lab/iparus-pipeline.git
cd iparus-pipeline

bash scripts/01_make_bams.sh data/samples.txt
bash scripts/02_sex_check.sh data/samples.txt
```

Each script has the paths it uses (reference, reads, BAM folder) at the top. If you follow [Installation](#installation), the defaults already fit.

---

## Software

| tool | version | used in | what it is for |
|---|---|---|---|
| Java | 17 or newer | step 1 | runs GATK and Trimmomatic |
| FastQC | 0.12.1 | step 1 | read quality reports |
| Trimmomatic | 0.39 | step 1 | adapter and quality trimming |
| bwa-mem2 | 2.2.1 | step 1 | mapping |
| samtools | 1.10 or newer | steps 1, 2, example | sorting, indexing, read depth |
| GATK | 4.6 | step 1 | read groups, merging, duplicate marking, metrics |
| R | any | step 1, example | insert-size histogram, FST tables |
| Python 3 | any | step 1, example | runs the `gatk` command and PCAngsd |
| md5sum | any (part of Linux) | checking downloads | checksums |
| ANGSD | 0.941 or newer | example | allele frequencies, FST, genotype likelihoods |
| PCAngsd | any recent | example | PCA |

All scripts are written for Linux and run from the terminal.

---

## Installation

**Java, R, Python and samtools** come from the system package manager. This needs admin rights; on a cluster, ask your admin or use `module load`:

```bash
sudo apt install openjdk-17-jre r-base python3 samtools
```

**The rest** are downloads that you unpack into a `tools` folder:

```bash
mkdir -p ~/tools && cd ~/tools

# FastQC
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip && chmod +x FastQC/fastqc

# Trimmomatic (comes with the adapter file step 1 needs)
wget https://github.com/usadellab/Trimmomatic/releases/download/v0.39/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip

# bwa-mem2
wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
tar -xjf bwa-mem2-2.2.1_x64-linux.tar.bz2

# GATK
wget https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip
unzip gatk-4.6.2.0.zip
```

Then make the programs findable by adding this line to the end of `~/.bashrc`, and open a new terminal:

```bash
export PATH=~/tools/FastQC:~/tools/bwa-mem2-2.2.1_x64-linux:~/tools/gatk-4.6.2.0:$PATH
```

If you can't use `apt` for samtools, build it from source instead (download page: https://www.htslib.org/download/):

```bash
cd ~/tools
wget https://github.com/samtools/samtools/releases/download/1.22.1/samtools-1.22.1.tar.bz2
tar -xjf samtools-1.22.1.tar.bz2 && cd samtools-1.22.1
./configure --prefix=$HOME/tools/samtools && make && make install
# then add ~/tools/samtools/bin to the PATH line above
```

ANGSD and PCAngsd are only needed for the example analysis; see [example_analysis/README.md](example_analysis/README.md#tools).

**Reference genome:** download *Parus major* 1.1 (GCF_001522545.3):

```bash
mkdir -p ~/iparus && cd ~/iparus
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/522/545/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.fna.gz
gunzip GCF_001522545.3_Parus_major1.1_genomic.fna.gz
```

---

## Data storage

The raw reads and BAM files are too large for GitHub.

| data | where |
|---|---|
| BAM files of Database Version 1.0 (+ index) | [IPARUS Database](https://gossmann-lab.github.io/iparus/) |
| raw reads | public archives (NCBI SRA / ENA), under the accessions as in `data/samples.txt` |
| reference genome and annotation | NCBI, GCF_001522545.3 (see [Installation](#installation)) |

Before using downloaded BAM files, check them against `data/bam_checksums.txt` (see [Checking downloaded BAM files](#checking-downloaded-bam-files)).

---

## Step 1: Making BAM files

`scripts/01_make_bams.sh` turns raw paired-end Illumina reads into one duplicate-marked BAM file per sample, aligned to *Parus major* 1.1. It runs:

1. a quality check (FastQC)
2. trimming (Trimmomatic)
3. mapping (bwa-mem2)
4. duplicate marking and quality metrics (GATK)

**Input.** Put the raw reads in `~/iparus/raw_reads`, either directly or one folder per sample, the way sequencing companies deliver them. The files must be named `SAMPLE_1.fq.gz` and `SAMPLE_2.fq.gz`; the folder names don't matter. For SRA runs, `fasterq-dump` followed by `gzip` gives these names.

```
~/iparus/
├── GCF_001522545.3_Parus_major1.1_genomic.fna     ← reference (+ index files made on the first run)
├── raw_reads/
│   ├── SRR1793405/
│   │   ├── SRR1793405_1.fq.gz
│   │   └── SRR1793405_2.fq.gz
│   └── SRR1793407/
│       ├── SRR1793407_1.fq.gz
│       └── SRR1793407_2.fq.gz
└── bams/                                          ← created by step 1
```

The sample list `data/samples.txt` has one name per line: the file name without `_1.fq.gz`, e.g. `SRR1793405`.

**Run:**

```bash
bash scripts/01_make_bams.sh data/samples.txt    # all samples in the list
bash scripts/01_make_bams.sh SRR1793405          # or just one sample
```

The script first checks that every tool is installed and that the reads can be found, and stops with a message if something is missing. The first run also indexes the reference. That happens only once, takes a few hours and needs about 40 GB of RAM. Mapping needs about 48 GB.

If a run stops partway through a sample, delete that sample's folder in `bams/` before running it again.

**Output:**

```
~/iparus/bams/SAMPLE/
├── SAMPLE_sorted_dedup.bam (+ .bai)    ← the final BAM
├── SAMPLE_dedup_metrics.txt            duplication rate
├── SAMPLE_alignment_metrics.txt        mapping rate etc.
├── SAMPLE_insert_metrics.txt
├── SAMPLE_insert_size_histogram.pdf
├── *.log                               log of each step
├── qc1/                                FastQC report of the raw reads
├── trim/                               trimmed reads
└── qc2/                                FastQC report of the trimmed reads
```

---

## Step 2: Sex check

`scripts/02_sex_check.sh` works out the genetic sex of each sample from its BAM file. It compares the mean read depth on the Z chromosome (`NC_031799.1`) with the depth on an autosome (`NC_031770.1`, chromosome 3).

In birds, females are ZW and so have one copy of Z, which gives a ratio of about **0.5**. Males are ZZ, which gives a ratio of about **1.0**.

**Input:** the BAMs from step 1, or BAMs downloaded from the website. The script finds them in either layout:

- `~/iparus/bams/SAMPLE/SAMPLE_sorted_dedup.bam` (made by step 1)
- `~/iparus/bams/SAMPLE_sorted_dedup.bam` (downloaded)

For a different folder, change `BAM_DIR` at the top of the script.

**Run:**

```bash
bash scripts/02_sex_check.sh data/samples.txt
bash scripts/02_sex_check.sh SRR2961674          # or just one sample
```

If a BAM has no index, the script creates one.

**Output:** `data/sex_results.tsv`, with one line per sample:

| column | meaning |
|---|---|
| `SampleID` | sample name |
| `Autosome_Depth` | mean depth on chromosome 3 |
| `Z_Depth` | mean depth on the Z chromosome |
| `Ratio_Z_Autosome` | Z depth / autosome depth |
| `Sex` | `female` (ratio ≤ 0.65), `male` (ratio ≥ 0.85), `unclear` (in between), `NoCoverage`, `FileNotFound`, or `Error` (the chromosomes are not in the BAM, e.g. it was aligned to a different reference) |

Example output. The file itself is tab-separated; it's shown aligned here for readability. The three females have a ratio of about 0.55, and the males about 1.0:

```
SampleID    Autosome_Depth  Z_Depth  Ratio_Z_Autosome  Sex
SRR765717   13.3287         7.32849  0.5498            female
SRR5423301  49.7291         50.5102  1.0157            male
SRR5423302  48.235          50.3875  1.0446            male
SRR9656355  15.7751         16.2848  1.0323            male
SRR9656356  11.9793         6.8109   0.5686            female
SRR9656357  13.0044         13.4324  1.0329            male
SRR9656367  12.5473         12.9261  1.0302            male
SRR9656368  16.2182         9.0897   0.5605            female
SRR9656369  14.8677         15.3608  1.0332            male
```

---

## Checking downloaded BAM files

Large downloads can be cut off or damaged without any error message. An MD5 checksum works like a fingerprint of a file: if even one byte differs, the checksum is different.

`data/bam_checksums.txt` lists the MD5 checksum of every BAM file in Database Version 1.0, as stored on sciebo. You don't need to create a checksum file yourself; you only compare your downloaded files against this list. Each line holds a checksum, two spaces, and then the file name:

```
3f2a9c1e5b7d...  SRR2961674_sorted_dedup.bam
9b41d0c7a2e8...  SRR2961676_sorted_dedup.bam
```

Go to the folder with your downloaded BAMs and run:

```bash
cd ~/iparus/bams
md5sum -c --ignore-missing /path/to/iparus-pipeline/data/bam_checksums.txt
```

`md5sum` computes the checksum of each of your files and compares it with the list. Each downloaded file is reported as `OK` or `FAILED`; download any file marked `FAILED` again. `--ignore-missing` skips the BAMs on the list that you haven't downloaded.

On macOS, which has no `md5sum`, run `md5 -r FILE` and compare the result with the line for that file in `data/bam_checksums.txt`.

As an extra check, `samtools quickcheck` confirms that each BAM is complete and readable. It prints nothing if all files are fine, and otherwise lists the broken ones:

```bash
samtools quickcheck -v *.bam
```

---

## Example analysis

[example_analysis/](example_analysis/README.md) shows what the IPARUS BAM files can be used for as in the manuscript. It compares 8 *Parus minor* from China with 8 European *Parus major* along chromosome 1A, using:

- FST in windows and per gene (`example_analysis/fst_chr1A.sh`)
- a PCA of 20 individuals, including 4 Asian *P. major* (`example_analysis/pca_chr1A.sh`)

Both use ANGSD genotype likelihoods. Its README covers the extra tools needed, every step, and every output file.
