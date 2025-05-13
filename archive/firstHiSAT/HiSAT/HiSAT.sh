#!/bin/bash

# Path to indexed Mycobacterium marinum genome
MMARINUM_INDEX="mmarinum_E11_index"  # Update if your index prefix is named differently

# Loop through paired-end files
for fq1 in *_1.fastq; do
    base=$(basename "$fq1" _1.fastq)
    fq2="${base}_2.fastq"
    echo "Processing $base"

    # Align paired reads
    hisat2 -x $MMARINUM_INDEX -1 "$fq1" -2 "$fq2" -S "${base}_mmarinum.sam"

    # Convert to sorted BAM
    samtools view -bS "${base}_mmarinum.sam" | samtools sort -o "${base}_mmarinum_sorted.bam"

    # Extract mapped reads only
    samtools view -b -F 4 "${base}_mmarinum_sorted.bam" > "${base}_mmarinum_mapped.bam"

    # Convert back to paired FASTQ files
    samtools fastq -1 "${base}_mmarinum_1.fastq" -2 "${base}_mmarinum_2.fastq" -0 /dev/null -s /dev/null "${base}_mmarinum_mapped.bam"

    # Clean up
    rm "${base}_mmarinum.sam" "${base}_mmarinum_sorted.bam" "${base}_mmarinum_mapped.bam"

    echo "Finished processing $base"
done

echo "All paired-end samples processed."
