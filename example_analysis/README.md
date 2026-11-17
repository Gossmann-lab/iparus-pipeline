# Example analysis: *Parus minor* vs. *Parus major* on chromosome 1A

[← back to the main README](../README.md)

This example shows what the iParus alignments can be used for as in the manuscript. It takes BAM files from the database and compares Japanese tits (*Parus minor*, China) with European great tits (*Parus major*) along chromosome 1A, in two ways:

- **Genetic divergence (FST)** across the whole chromosome, in windows, and per gene. Script: [`fst_chr1A.sh`](fst_chr1A.sh).
- **Population structure (PCA)** of 20 individuals: the two groups above plus four Asian *P. major*. Script: [`pca_chr1A.sh`](pca_chr1A.sh).

Both analyses work on genotype likelihoods (ANGSD) instead of called genotypes, so moderate- and low-coverage samples can be used directly.

## Samples

| group | n | accessions |
|---|---|---|
| *P. minor*, China | 8 | SRR1793405, SRR1793407, SRR1793437, SRR1793438, SRR1793439, SRR1793445, SRR9656357, SRR9656369 |
| *P. major*, Europe | 8 | SRR2961674, SRR2961676, SRR2961680, SRR2961689, SRR2961699, SRR2961701, SRR2961712, SRR2961728 |
| *P. major*, Asia (PCA only) | 4 | SRR1793428, SRR1793430, SRR1793432, SRR1793435 |

The two FST groups each have eight males, one per sampling locality. If you use your own samples, the [sex check](../README.md#step-2-sex-check) tells you which are male and which are female.

## Tools

You need samtools and R (see [Installation](../README.md#installation) in the main README), plus:

| tool | version | used for |
|---|---|---|
| ANGSD (with `realSFS`) | 0.941 or newer | site allele frequencies, 2D-SFS, FST, genotype likelihoods |
| PCAngsd | any recent | covariance matrix for the PCA |
| bc | any | small calculations in the scripts (usually already installed) |

Install them like this:

```bash
cd ~/tools

# ANGSD (realSFS is built with it)
git clone --recursive https://github.com/samtools/htslib.git
git clone https://github.com/ANGSD/angsd.git
cd htslib && make && cd ../angsd && make HTSSRC=../htslib && cd ..

# PCAngsd (needs Python 3 and pip)
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd && pip3 install . && cd ..
```

Then add ANGSD to your PATH by adding this line to the end of `~/.bashrc`:

```bash
export PATH=~/tools/angsd:~/tools/angsd/misc:$PATH
```

## Input

Put everything in one working folder (`WORKING_DIR` at the top of both scripts):

```
iparus/
├── GCF_001522545.3_Parus_major1.1_genomic.fna   ← reference
├── GCF_001522545.3_Parus_major1.1_genomic.gff   ← RefSeq annotation (same FTP folder as the reference)
└── bams/
    ├── SRR1793405_sorted_dedup.bam   (+ .bai)
    └── ...                           ← the 20 BAMs from the iParus database
```

After downloading the BAMs, check that they arrived intact (see [Step 3: BAM checksums](../README.md#step-3-bam-checksums)).

The annotation is in the same NCBI folder as the reference:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/522/545/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.gff.gz
gunzip GCF_001522545.3_Parus_major1.1_genomic.gff.gz
```

## Running

```bash
bash example_analysis/fst_chr1A.sh
bash example_analysis/pca_chr1A.sh
```

Both scripts use 16 threads (`THREADS` at the top). Neither makes figures: they write tables, and plotting is done separately.

**Other chromosomes:** change `CHROMS` at the top of both scripts. The accessions of the main chromosomes are listed in `fst_chr1A.sh`, e.g. `NC_031768.1` is chromosome 1 and `NC_031799.1` is the Z. Keep the two scripts on the same chromosomes. If you change chromosome, also update `REGION` and `CONTROL` in `fst_chr1A.sh`, since they must lie on an analysed chromosome. A single chromosome gives far fewer SNPs than the whole genome: `pca_chr1A.sh` warns when there are fewer than 50,000, because the components can then be unstable.

## What the FST script does

1. **Site allele frequency likelihoods** for each group (`angsd -doSaf`). Filters: mapping quality ≥ 30, base quality ≥ 20, data in at least 6 of the 8 individuals, only uniquely mapped reads in proper pairs.
2. **Joint site frequency spectrum** of the two groups (`realSFS`), used as the prior for FST.
3. **FST** (`realSFS fst`):
   - over the whole region; the weighted value (ratio of averages) is the one reported
   - in 50 kb windows with a 10 kb step, for an overview along the chromosome
   - in 5 kb windows with a 1 kb step, for gene-level detail
4. **Mean FST per chromosome**, from non-overlapping 50 kb windows with at least 20,000 sites.
5. **FST per gene**: sum(α)/sum(β) over all sites inside each annotated gene (RefSeq GFF). This is the same estimator as the global value, so no windowed value has to be assigned to a gene.
6. **Ranking genes against their background.** Raw FST favours genes that lie in already-diverged regions, so each gene is compared with the mean FST of genes within 1 Mb of it. Genes are ranked by this excess. Only genes with at least 2,000 sites are kept.
7. **Checks that the signal is not an artefact:**
   - read depth in the region of interest (`REGION`, around *MOV10L1*) relative to a nearby control region (`CONTROL`), for every individual; a copy-number or mapping difference between the groups would show up as different ratios
   - the share of reads in the region that map uniquely (MAPQ ≥ 30), for one individual per group

## What the PCA script does

1. **Genotype likelihoods** for all 20 individuals in a single joint ANGSD run, so that major and minor alleles are assigned the same way in every sample. Only SNPs are kept (p < 1e-6, minor allele frequency ≥ 0.05, data in at least 15 of the 20 individuals); the other filters are as for FST.
2. **Covariance matrix** with PCAngsd. Eigen-decomposing this matrix gives the principal components.

## Output

`fst_chr1A/`:

| file | contents |
|---|---|
| `global.fst` | unweighted and weighted FST over the analysed region |
| `windows_50kb.txt`, `windows_5kb.txt` | windowed FST (region, chromosome, midpoint, number of sites, FST) |
| `chr_means.txt` | mean FST per chromosome |
| `fst_per_gene.txt` | FST per annotated gene, with number of sites |
| `genes_ranked.csv` | genes with their FST, local background FST and excess, sorted by excess |
| `depth.tsv` | per-individual depth in control and target region, and their ratio |
| `mapq.txt` | number of reads in the region at MAPQ ≥ 0 and ≥ 30 |

`pca_chr1A/`:

| file | contents |
|---|---|
| `all20_pca.cov` | covariance matrix (20 × 20) |
| `all20.bamlist` | sample order; row *i* of the matrix is line *i* of this file |
| `all20.beagle.gz` | genotype likelihoods of the SNPs used |
