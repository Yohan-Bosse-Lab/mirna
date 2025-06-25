#!/bin/bash


# Variables
REFERENCE_GENOME='/mnt/sde/renseb01/Documents/Steinberg_Christian/mirna/data/reference/thymine/Hsa_mature_thymine.fa'
REFERENCE_INDEX='/mnt/sde/renseb01/Documents/Steinberg_Christian/mirna/data/reference/bowtie2/Hsa'
FASTQDIR='/mnt/sde/renseb01/Documents/Steinberg_Christian/mirna/data/fastq_trim_REFRAME_2025/'
OUTDIR='/mnt/sde/renseb01/Documents/Steinberg_Christian/mirna/results/bams_REFRAME_2025/'
THREADS=10

#for loop
for READS1 in $FASTQDIR*_trimmed_R1.fq.gz
  do
    SAMPLE1=${READS1//$FASTQDIR/}
    TEMP1=${SAMPLE1//'_trimmed_R1.fq.gz'/''}
    PREFIX=${TEMP1//'NS.X0180.004.NEBNext_dual_i7_'/''}

    READS2=${READS1//'_trimmed_R1.fq.gz'/'_trimmed_R2.fq.gz'}
    SAMPLE2=${READS2//$FASTQDIR/}
    TEMP2=${SAMPLE2//'_trimmed_R1.fq.gz'/''}
    PREFIX2=${TEMP2//'NS.X0180.004.NEBNext_dual_i7_'/''}

    echo $READS1
    echo $READS2
    echo $SAMPLE1
    echo $TEMP1
    echo $PREFIX
    echo $REFERENCE_INDEX
    echo "$OUTDIR${PREFIX}.sam"

    # Perform the alignment (for single-end reads)
    echo "Aligning paired-end reads to reference genome..."
    bowtie2 -x $REFERENCE_INDEX -1 $READS1 -2 $READS2 -N --end-to-end >$OUTDIR${PREFIX}.sam

    # Convert SAM to BAM (optional but recommended for downstream processing)
    echo "Converting SAM to BAM..."
    samtools view -Sb "$OUTDIR${PREFIX}.sam" > "$OUTDIR${PREFIX}.bam"

    # Sort the BAM file (optional)
    echo "Sorting BAM file..."
    samtools sort "$OUTDIR${PREFIX}.bam" -o "$OUTDIR${PREFIX}_sorted.bam"

    # Index the sorted BAM file (optional)
    echo "Indexing the sorted BAM file..."
    samtools index "$OUTDIR${PREFIX}_sorted.bam"

    # Generate the expression data
    samtools idxstats "$OUTDIR${PREFIX}_sorted.bam" >"$OUTDIR${PREFIX}.idxstats"

    # Clean up: Remove the original SAM file if desired
    rm "$OUTDIR${PREFIX}.sam"
    rm "$OUTDIR${PREFIX}.bam"

    echo "Alignment completed. Output files:"
    echo "- ${PREFIX}_sorted.bam"
    echo "- ${PREFIX}_sorted.bam.bai"  # BAM index file
  done
