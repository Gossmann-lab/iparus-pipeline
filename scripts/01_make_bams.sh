#!/bin/bash
# IPARUS step 1: raw paired-end reads -> duplicate-marked BAM, aligned to Parus major 1.1
#
# Run from the repository folder:
#       bash scripts/01_make_bams.sh data/samples.txt   (one sample name per line)
#   or: bash scripts/01_make_bams.sh SAMPLE_NAME        (a single sample)
#
# Raw reads must be named SAMPLE_1.fq.gz and SAMPLE_2.fq.gz and be somewhere in
# FASTQ_DIR, either directly or in a subfolder (e.g. one folder per sample).
# Result: OUT_DIR/SAMPLE/SAMPLE_sorted_dedup.bam

set -euo pipefail

# ---- edit these paths ----
REF=$HOME/iparus/GCF_001522545.3_Parus_major1.1_genomic.fna
FASTQ_DIR=$HOME/iparus/raw_reads
OUT_DIR=$HOME/iparus/bams
TRIMMOMATIC=$HOME/tools/Trimmomatic-0.39/trimmomatic-0.39.jar
ADAPTERS=$HOME/tools/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa
# --------------------------

# check that all tools are installed
for tool in java fastqc bwa-mem2 samtools gatk Rscript; do
    command -v $tool > /dev/null || { echo "$tool is not installed (see README)"; exit 1; }
done
for file in "$TRIMMOMATIC" "$ADAPTERS" "$REF"; do
    [ -f "$file" ] || { echo "$file not found (see README)"; exit 1; }
done

# sample list file or single sample name
if [ -f "${1:?give a sample list file or a sample name}" ]; then SAMPLES=$(grep -v '^#' "$1"); else SAMPLES=$1; fi

# index the reference (only the first time)
[ -f "$REF.bwt.2bit.64" ] || bwa-mem2 index "$REF"
[ -f "$REF.fai" ] || samtools faidx "$REF"
[ -f "${REF%.*}.dict" ] || gatk CreateSequenceDictionary -R "$REF" -O "${REF%.*}.dict"

TMP=$OUT_DIR/tmp
mkdir -p "$TMP"
JAVA="-Xmx40G -Xms40G -Djava.io.tmpdir=$TMP"

for S in $SAMPLES; do
    echo "=== $S"
    D=$OUT_DIR/$S
    mkdir -p "$D/qc1" "$D/trim" "$D/qc2"
    # find the reads anywhere under FASTQ_DIR (e.g. in a folder per sample)
    R1=$(find "$FASTQ_DIR" -name "${S}_1.fq.gz" | head -n 1)
    R2=$(find "$FASTQ_DIR" -name "${S}_2.fq.gz" | head -n 1)
    [ -n "$R1" ] && [ -n "$R2" ] || { echo "Reads for $S not found in $FASTQ_DIR"; exit 1; }

    # 1. quality check and trimming
    fastqc -o "$D/qc1" -f fastq "$R1" "$R2" > "$D/qc1/$S.log" 2>&1
    java -jar "$TRIMMOMATIC" PE -threads 5 -phred33 -trimlog "$D/trim/$S.txt" "$R1" "$R2" \
        -baseout "$D/trim/trim_$S.fq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        > "$D/trim/trim_$S.log" 2>&1
    rm "$D/trim/$S.txt"
    fastqc -o "$D/qc2" -f fastq "$D/trim/trim_${S}_1P.fq.gz" "$D/trim/trim_${S}_2P.fq.gz" > "$D/qc2/$S.log" 2>&1

    # 2. mapping
    gatk --java-options "$JAVA" FastqToSam --FASTQ "$D/trim/trim_${S}_1P.fq.gz" --FASTQ2 "$D/trim/trim_${S}_2P.fq.gz" \
        --OUTPUT "$D/${S}_fq2sam.bam" --READ_GROUP_NAME "$S" --SAMPLE_NAME "$S" --LIBRARY_NAME "$S" \
        --PLATFORM ILLUMINA --PREDICTED_INSERT_SIZE 100 --SORT_ORDER queryname > "$D/${S}_rg.log" 2>&1
    gatk --java-options "$JAVA" MarkIlluminaAdapters --INPUT "$D/${S}_fq2sam.bam" --OUTPUT "$D/${S}_nwxt.bam" \
        --METRICS "$D/${S}_nwxt_metrics.txt" >> "$D/${S}_rg.log" 2>&1
    gatk --java-options "$JAVA" SamToFastq --INPUT "$D/${S}_nwxt.bam" --FASTQ "$D/${S}_sam2fq.fastq" \
        --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --INTERLEAVE true --INCLUDE_NON_PF_READS true >> "$D/${S}_rg.log" 2>&1
    rm "$D/${S}_nwxt.bam"

    bwa-mem2 mem -t 1 -p "$REF" "$D/${S}_sam2fq.fastq" -o "$D/$S.sam" > "$D/${S}_bwa.log" 2>&1
    rm "$D/${S}_sam2fq.fastq"
    samtools sort -@ 1 -T "$TMP/$S" -o "$D/${S}_alned.bam" "$D/$S.sam"
    rm "$D/$S.sam"

    gatk --java-options "$JAVA" MergeBamAlignment --REFERENCE_SEQUENCE "$REF" --UNMAPPED_BAM "$D/${S}_fq2sam.bam" \
        --ALIGNED_BAM "$D/${S}_alned.bam" --OUTPUT "$D/${S}_merged.bam" --CREATE_INDEX false --ADD_MATE_CIGAR true \
        --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true --INCLUDE_SECONDARY_ALIGNMENTS true \
        --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS \
        > "$D/${S}_merge.log" 2>&1
    rm "$D/${S}_fq2sam.bam" "$D/${S}_alned.bam"
    samtools sort -@ 1 -n -T "$TMP/$S" "$D/${S}_merged.bam" -o "$D/${S}_merged_sorted.bam"
    rm "$D/${S}_merged.bam"

    # 3. mark duplicates and collect metrics
    gatk --java-options "-Xmx32G -Xms32G -Djava.io.tmpdir=$TMP" MarkDuplicatesSpark -I "$D/${S}_merged_sorted.bam" \
        -M "$D/${S}_dedup_metrics.txt" -O "$D/${S}_sorted_dedup.bam" --tmp-dir "$TMP" > "$D/${S}_mrkdup.log" 2>&1
    rm "$D/${S}_merged_sorted.bam"
    gatk CollectAlignmentSummaryMetrics -R "$REF" -I "$D/${S}_sorted_dedup.bam" \
        -O "$D/${S}_alignment_metrics.txt" > "$D/${S}_colmx.log" 2>&1
    gatk CollectInsertSizeMetrics --INPUT "$D/${S}_sorted_dedup.bam" --OUTPUT "$D/${S}_insert_metrics.txt" \
        --Histogram_FILE "$D/${S}_insert_size_histogram.pdf" >> "$D/${S}_colmx.log" 2>&1

    echo "Done: $D/${S}_sorted_dedup.bam"
done
