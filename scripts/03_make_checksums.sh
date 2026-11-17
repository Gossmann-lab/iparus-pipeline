#!/bin/bash
# IPARUS step 3: MD5 checksum of every BAM file, so that people who download
# the BAMs can check that their copies are complete and unchanged.
#
# Run from the repository folder:
#       bash scripts/03_make_checksums.sh
# Result: data/bam_checksums.txt   (format of md5sum: checksum, two spaces, file name)

# ---- edit these ----
BAM_DIR=$HOME/iparus/bams
OUTPUT=data/bam_checksums.txt
# --------------------

OUTPUT=$(realpath "$OUTPUT")
: > "$OUTPUT"

# file names only (no folders), so the list also works for downloaded BAMs
# that are all stored in one folder
find "$BAM_DIR" -name '*_sorted_dedup.bam' | sort | while read -r BAM; do
    (cd "$(dirname "$BAM")" && md5sum "$(basename "$BAM")") | tee -a "$OUTPUT"
done

echo "Checksums of $(wc -l < "$OUTPUT") BAM files saved in $OUTPUT"
