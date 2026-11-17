#!/usr/bin/env bash
#
# Principal component analysis of twenty individuals, using
# alignments retrieved from IPARUS.
#
# The chromosomes analysed are set by CHROMS below, and should
# match those used in fst_chr1A.sh.
#
# Output, in the directory named by OUT
#   all20_pca.cov     covariance matrix
#   all20.bamlist     sample order, needed when plotting
#
# Requires: ANGSD >= 0.941, PCAngsd

set -euo pipefail

THREADS=16
WORKING_DIR="/vol/storage/iparus"
REF="GCF_001522545.3_Parus_major1.1_genomic.fna"
OUT="pca_chr1A"

# RefSeq accessions; see fst_chr1A.sh for the full list.
CHROMS=( "NC_031773.1" )        # chromosome 1A

# Require data in at least this many of the 20 individuals.
MININD=15

cd "$WORKING_DIR"
mkdir -p "$OUT"

for prog in angsd pcangsd samtools; do
    command -v "$prog" >/dev/null 2>&1 || {
        echo "ERROR: $prog not in PATH"; exit 1; }
done

# ------------------------------------------------------------
# Samples
#
# Order matters: the plotting script matches each row of the
# covariance matrix to this list. Eight P. minor, then eight
# European P. major, then four Asian P. major.
# ------------------------------------------------------------

cat > "$OUT/all20.bamlist" <<'EOF'
bams/SRR1793405_sorted_dedup.bam
bams/SRR1793407_sorted_dedup.bam
bams/SRR1793437_sorted_dedup.bam
bams/SRR1793438_sorted_dedup.bam
bams/SRR1793439_sorted_dedup.bam
bams/SRR1793445_sorted_dedup.bam
bams/SRR9656357_sorted_dedup.bam
bams/SRR9656369_sorted_dedup.bam
bams/SRR2961674_sorted_dedup.bam
bams/SRR2961676_sorted_dedup.bam
bams/SRR2961680_sorted_dedup.bam
bams/SRR2961689_sorted_dedup.bam
bams/SRR2961699_sorted_dedup.bam
bams/SRR2961701_sorted_dedup.bam
bams/SRR2961712_sorted_dedup.bam
bams/SRR2961728_sorted_dedup.bam
bams/SRR1793428_sorted_dedup.bam
bams/SRR1793430_sorted_dedup.bam
bams/SRR1793432_sorted_dedup.bam
bams/SRR1793435_sorted_dedup.bam
EOF

printf '%s\n' "${CHROMS[@]}" > "$OUT/regions.txt"

[[ -f "${REF}.fai" ]] || samtools faidx "$REF"
echo "analysing: $(paste -sd' ' "$OUT/regions.txt")"

# ------------------------------------------------------------
# Genotype likelihoods
#
# A single joint run over all twenty individuals: major and
# minor alleles must be assigned consistently across samples,
# so per-group runs cannot be combined afterwards. Restricted
# to variable sites, unlike the site allele frequency
# likelihoods used for FST, which need the full spectrum.
# ------------------------------------------------------------

angsd \
    -bam "$OUT/all20.bamlist" \
    -ref "$REF" \
    -rf "$OUT/regions.txt" \
    -out "$OUT/all20" \
    -nThreads "$THREADS" \
    -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doCounts 1 \
    -SNP_pval 1e-6 -minMaf 0.05 \
    -minMapQ 30 -minQ 20 -minInd "$MININD" \
    -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
    -skipTriallelic 1

NSNP=$(( $(zcat "$OUT/all20.beagle.gz" | wc -l) - 1 ))
echo "SNPs: $NSNP"

# A single chromosome yields far fewer SNPs than a genome-wide
# run. Below roughly 50,000 the components become unstable with
# twenty individuals.
if [[ "$NSNP" -lt 50000 ]]; then
    echo "WARNING: only $NSNP SNPs - components may be unstable"
fi

# ------------------------------------------------------------
# Covariance matrix
# ------------------------------------------------------------

pcangsd -b "$OUT/all20.beagle.gz" -o "$OUT/all20_pca" -t "$THREADS"

echo
echo "Done. Covariance matrix: $OUT/all20_pca.cov"
echo "Sample order:           $OUT/all20.bamlist"
