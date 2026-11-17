#!/usr/bin/env bash
#
# Divergence between Parus minor (China) and European P. major,
# using alignments retrieved from IPARUS.
#
# The chromosomes analysed are set by CHROMS below. As written
# the analysis covers chromosome 1A; add accessions to the list
# to extend it, or replace them to analyse a different region.
#
# The principal component analysis is in pca_chr1A.sh.
# Figures are produced separately; this script does no plotting.
#
# Input
#   bams/          alignments, indexed
#   reference FASTA and RefSeq GFF (paths below)
#
# Output, all in the directory named by OUT
#   global.fst             mean FST over the analysed region
#   windows_50kb.txt       50 kb windows, 10 kb step
#   windows_5kb.txt         5 kb windows,  1 kb step
#   chr_means.txt          mean FST per chromosome analysed
#   fst_per_gene.txt       FST within annotated gene bounds
#   genes_ranked.csv       genes ranked by excess over local
#   depth.tsv              per-individual depth, region/control
#   mapq.txt               unique-mapping check
#
# Requires: ANGSD >= 0.941 (angsd, realSFS), samtools, R

set -euo pipefail

THREADS=16
WORKING_DIR="/vol/storage/iparus"
REF="GCF_001522545.3_Parus_major1.1_genomic.fna"
GFF="GCF_001522545.3_Parus_major1.1_genomic.gff"
OUT="fst_chr1A"

# ------------------------------------------------------------
# Chromosomes to analyse
#
# RefSeq accessions in Parus_major1.1:
#   NC_031768.1  chromosome 1     114.1 Mb
#   NC_031769.1  chromosome 2     150.3 Mb
#   NC_031770.1  chromosome 3     111.6 Mb
#   NC_031773.1  chromosome 1A     71.4 Mb
#   NC_031792.1  chromosome 24      6.4 Mb
#   NC_031799.1  chromosome Z      74.5 Mb
#
# Names used in output and figures.
# ------------------------------------------------------------

CHROMS=( "NC_031773.1" )
declare -A CHRNAME=(
  [NC_031768.1]="1"  [NC_031769.1]="2" [NC_031770.1]="3"
  [NC_031773.1]="1A" [NC_031792.1]="24" [NC_031799.1]="Z"
)

# Region used for the depth comparison, and a nearby control.
# Both must lie within the analysed chromosomes.
REGION="NC_031773.1:18511000-18528000"    # MOV10L1 plateau
CONTROL="NC_031773.1:18400000-18432000"   # background

cd "$WORKING_DIR"
mkdir -p "$OUT"

for prog in angsd realSFS samtools Rscript; do
    command -v "$prog" >/dev/null 2>&1 || {
        echo "ERROR: $prog not in PATH"; exit 1; }
done

# ------------------------------------------------------------
# 1. Sample lists and regions
#
# Eight males per group, one individual per sampling locality.
# ------------------------------------------------------------

cat > "$OUT/minor.bamlist" <<'EOF'
bams/SRR1793405_sorted_dedup.bam
bams/SRR1793407_sorted_dedup.bam
bams/SRR1793437_sorted_dedup.bam
bams/SRR1793438_sorted_dedup.bam
bams/SRR1793439_sorted_dedup.bam
bams/SRR1793445_sorted_dedup.bam
bams/SRR9656357_sorted_dedup.bam
bams/SRR9656369_sorted_dedup.bam
EOF

cat > "$OUT/major.bamlist" <<'EOF'
bams/SRR2961674_sorted_dedup.bam
bams/SRR2961676_sorted_dedup.bam
bams/SRR2961680_sorted_dedup.bam
bams/SRR2961689_sorted_dedup.bam
bams/SRR2961699_sorted_dedup.bam
bams/SRR2961701_sorted_dedup.bam
bams/SRR2961712_sorted_dedup.bam
bams/SRR2961728_sorted_dedup.bam
EOF

printf '%s\n' "${CHROMS[@]}" > "$OUT/regions.txt"

[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

echo "analysing:"
while read -r c; do
    len=$(awk -v c="$c" '$1==c {print $2}' "${REF}.fai")
    [[ -n "$len" ]] || { echo "ERROR: $c not in ${REF}.fai"; exit 1; }
    printf '  %-14s chromosome %-3s %6.1f Mb\n' \
           "$c" "${CHRNAME[$c]:-?}" "$(echo "$len / 1000000" | bc -l)"
done < "$OUT/regions.txt"

# ------------------------------------------------------------
# 2. Site allele frequency likelihoods, per group
# ------------------------------------------------------------

for grp in minor major; do
    angsd \
        -bam "$OUT/${grp}.bamlist" \
        -ref "$REF" -anc "$REF" \
        -rf "$OUT/regions.txt" \
        -out "$OUT/$grp" \
        -nThreads "$THREADS" \
        -GL 2 -doSaf 1 -doMajorMinor 1 -doMaf 1 -doCounts 1 \
        -minMapQ 30 -minQ 20 -minInd 6 \
        -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
        -skipTriallelic 1
done

# ------------------------------------------------------------
# 3. Two-dimensional SFS and FST
# ------------------------------------------------------------

realSFS "$OUT/minor.saf.idx" "$OUT/major.saf.idx" \
    -P "$THREADS" > "$OUT/minor_major.2dsfs"

realSFS fst index "$OUT/minor.saf.idx" "$OUT/major.saf.idx" \
    -sfs "$OUT/minor_major.2dsfs" -fstout "$OUT/minor_major"

# Unweighted and weighted. Weighted, the ratio of averages, is
# the value reported.
realSFS fst stats "$OUT/minor_major.fst.idx" | tee "$OUT/global.fst"

# 50 kb windows for the overview, 5 kb for resolving genes.
realSFS fst stats2 "$OUT/minor_major.fst.idx" -win 50000 -step 10000 \
    > "$OUT/windows_50kb.txt"

realSFS fst stats2 "$OUT/minor_major.fst.idx" -win 5000 -step 1000 \
    > "$OUT/windows_5kb.txt"

# ------------------------------------------------------------
# 4. Mean FST per chromosome analysed
# ------------------------------------------------------------

Rscript - "$OUT" <<'RSCRIPT' | tee "$OUT/chr_means.txt"
args <- commandArgs(trailingOnly = TRUE)
OUT  <- args[1]

w <- read.table(file.path(OUT, "windows_50kb.txt"), header = FALSE,
                skip = 1, stringsAsFactors = FALSE)
names(w)[1:5] <- c("region","chr","midPos","Nsites","Fst")
w$Fst <- pmax(as.numeric(w$Fst), 0)
w <- w[is.finite(w$Fst) & w$Nsites >= 20000, ]
w <- w[seq(1, nrow(w), by = 5), ]          # non-overlapping

nm <- c(NC_031768.1 = "1",  NC_031769.1 = "2",  NC_031770.1 = "3",
        NC_031773.1 = "1A", NC_031792.1 = "24", NC_031799.1 = "Z")

m <- tapply(w$Fst, nm[w$chr], mean)
for (k in names(m)) cat(sprintf("chr %-3s %.3f  (%d windows)\n",
                                k, m[[k]], sum(nm[w$chr] == k)))
RSCRIPT

# ------------------------------------------------------------
# 5. Per-gene FST
#
# FST for a region is sum(alpha)/sum(beta) over its sites --
# the same ratio-of-averages estimator used for the overall
# value, applied within annotated gene boundaries. This avoids
# the ambiguity of assigning a windowed value to a gene.
# ------------------------------------------------------------

awk -F'\t' -v CHROMS="$(printf '%s,' "${CHROMS[@]}")" '
BEGIN { n = split(CHROMS, a, ","); for (i=1; i<=n; i++) keep[a[i]] = 1 }
$3=="gene" && ($1 in keep) {
    match($9,/gene=[^;]*/);        g = RSTART ? substr($9,RSTART+5,RLENGTH-5)  : "."
    match($9,/gene_biotype=[^;]*/); b = RSTART ? substr($9,RSTART+13,RLENGTH-13) : "."
    print $1"\t"$4"\t"$5"\t"g"\t"b
}' "$GFF" > "$OUT/genes.bed"

echo "annotated genes in the analysed region: $(wc -l < "$OUT/genes.bed")"

realSFS fst print "$OUT/minor_major.fst.idx" 2>/dev/null \
  | awk -v GENES="$OUT/genes.bed" '
    BEGIN {
        FS="\t"; BIN=100000
        while ((getline line < GENES) > 0) {
            split(line, f, "\t"); n++
            chr[n]=f[1]; s[n]=f[2]; e[n]=f[3]; nm[n]=f[4]; bt[n]=f[5]
            for (b = int(f[2]/BIN); b <= int(f[3]/BIN); b++)
                idx[f[1], b] = idx[f[1], b] " " n
        }
        close(GENES)
    }
    {
        key = $1 SUBSEP int($2/BIN)
        if (!(key in idx)) next
        m = split(idx[key], c, " ")
        for (i=1; i<=m; i++) {
            j = c[i]
            if (j != "" && $2 >= s[j] && $2 <= e[j]) {
                A[j] += $3; B[j] += $4; N[j]++
            }
        }
    }
    END {
        print "gene\tbiotype\tchr\tstart\tend\tnsites\tfst"
        for (j=1; j<=n; j++)
            if (N[j] > 0 && B[j] > 0)
                printf "%s\t%s\t%s\t%d\t%d\t%d\t%.6f\n",
                       nm[j], bt[j], chr[j], s[j], e[j], N[j], A[j]/B[j]
    }' > "$OUT/fst_per_gene.txt"

# ------------------------------------------------------------
# 6. Rank genes against their local background
#
# Raw FST favours genes lying in already elevated regions, so
# each gene is also compared with the mean FST of genes within
# 1 Mb, excluding itself.
# ------------------------------------------------------------

Rscript - "$OUT" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
OUT  <- args[1]

g <- read.delim(file.path(OUT, "fst_per_gene.txt"),
                stringsAsFactors = FALSE)
g <- g[g$nsites >= 2000, ]
g$mid <- (g$start + g$end) / 2

g$local <- NA_real_
for (ch in unique(g$chr)) {
  i <- which(g$chr == ch); sub <- g[i, ]
  for (k in seq_along(i)) {
    near <- abs(sub$mid - sub$mid[k]) <= 1e6; near[k] <- FALSE
    if (sum(near) >= 3) g$local[i[k]] <- mean(sub$fst[near])
  }
}
g$excess <- g$fst - g$local

write.csv(g[order(-g$excess), ], file.path(OUT, "genes_ranked.csv"),
          row.names = FALSE)

cat(sprintf("genes with >= 2000 sites: %d\n\n", nrow(g)))
print(head(g[order(-g$excess),
            c("gene","biotype","chr","start","end","nsites",
              "fst","local","excess")], 20),
      row.names = FALSE, digits = 3)

cat("\n")
print(g[g$gene %in% c("MOV10L1","MLC1","PANX2"),
        c("gene","start","end","nsites","fst","local","excess")],
      row.names = FALSE, digits = 3)
RSCRIPT

# ------------------------------------------------------------
# 7. Read depth across the signal, normalised per individual
#
# Confirms the signal is not driven by differential mapping or
# by a copy-number difference between groups.
# ------------------------------------------------------------

meandepth() {
    samtools depth -a -q 20 -Q 30 -r "$1" "$2" \
      | awk '{t+=$3} END{if (NR>0) printf "%.2f", t/NR; else printf "NA"}'
}

printf 'group\tsample\tcontrol\tregion\tratio\n' > "$OUT/depth.tsv"
for grp in minor major; do
    while read -r bam; do
        [[ -n "$bam" ]] || continue
        n=$(basename "$bam" _sorted_dedup.bam)
        c=$(meandepth "$CONTROL" "$bam")
        p=$(meandepth "$REGION"  "$bam")
        printf '%s\t%s\t%s\t%s\t%.2f\n' "$grp" "$n" "$c" "$p" \
               "$(echo "$p / $c" | bc -l)" >> "$OUT/depth.tsv"
    done < "$OUT/${grp}.bamlist"
done
column -t "$OUT/depth.tsv"

# Proportion of reads mapping uniquely over the region,
# one individual per group.
{
  for s in SRR1793439 SRR2961674; do
      printf '%-12s Q0=%s Q30=%s\n' "$s" \
        "$(samtools view -c      bams/${s}_sorted_dedup.bam $REGION)" \
        "$(samtools view -c -q 30 bams/${s}_sorted_dedup.bam $REGION)"
  done
} | tee "$OUT/mapq.txt"

echo
echo "Done. Figures are produced separately from the files in $OUT/."
