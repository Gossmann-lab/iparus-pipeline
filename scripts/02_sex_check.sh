#!/bin/bash
# IPARUS step 2: genetic sex from read depth on the Z chromosome (NC_031799.1)
# divided by read depth on an autosome (NC_031770.1).
# Birds: ~0.5 = female (ZW), ~1.0 = male (ZZ).
#
# Run from the repository folder:
#       bash scripts/02_sex_check.sh data/samples.txt   (same list as step 1)
#   or: bash scripts/02_sex_check.sh SAMPLE_NAME
# Result: data/sex_results.tsv
#
# Reads BAM_DIR/SAMPLE/SAMPLE_sorted_dedup.bam or BAM_DIR/SAMPLE_sorted_dedup.bam

# ---- edit these ----
BAM_DIR=$HOME/iparus/bams
OUTPUT=sex_results.tsv
# --------------------

command -v samtools > /dev/null || { echo "samtools is not installed (see README)"; exit 1; }

if [ -f "${1:?give a sample list file or a sample name}" ]; then SAMPLES=$(grep -v '^#' "$1"); else SAMPLES=$1; fi

echo -e "SampleID\tAutosome_Depth\tZ_Depth\tRatio_Z_Autosome\tSex" > "$OUTPUT"

for S in $SAMPLES; do
    BAM=$BAM_DIR/$S/${S}_sorted_dedup.bam                         # made with iparus_bam.sh
    [ -f "$BAM" ] || BAM=$BAM_DIR/${S}_sorted_dedup.bam            # downloaded from the website

    if [ ! -f "$BAM" ]; then
        echo "$S: BAM not found"
        echo -e "$S\tNA\tNA\tNA\tFileNotFound" >> "$OUTPUT"
        continue
    fi
    [ -f "$BAM.bai" ] || [ -f "${BAM%.bam}.bai" ] || samtools index "$BAM"

    AUTO=$(samtools coverage -r NC_031770.1 "$BAM" 2>/dev/null | awk 'NR==2 {print $7}')
    Z=$(samtools coverage -r NC_031799.1 "$BAM" 2>/dev/null | awk 'NR==2 {print $7}')

    if [ -z "$AUTO" ] || [ -z "$Z" ]; then
        echo "$S: chromosomes not found in BAM"
        echo -e "$S\tNA\tNA\tNA\tError" >> "$OUTPUT"
        continue
    fi

    # female if ratio <= 0.65, male if >= 0.85, otherwise unclear
    awk -v s="$S" -v a="$AUTO" -v z="$Z" 'BEGIN {
        r = (a > 0) ? z / a : 0
        sex = (a == 0) ? "NoCoverage" : (r <= 0.65) ? "female" : (r >= 0.85) ? "male" : "unclear"
        printf "%s\t%s\t%s\t%.4f\t%s\n", s, a, z, r, sex
    }' | tee -a "$OUTPUT"
done

echo "Results saved in $OUTPUT"
