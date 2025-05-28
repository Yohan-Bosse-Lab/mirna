#!/bin/bash


# Variables
REFERENCE_GENOME='/mnt/sde/renseb01/Documents/Steinberg_Christian/mirna/data/reference/thymine/Hsa_mature_thymine.fa'
REFERENCE_INDEX='/mnt/sde/renseb01/Documents/Steinberg_Christian/mirna/data/reference/thymine/reference_index'
FASTQDIR='/mnt/sde/renseb01/Documents/Steinberg_Christian/mirna/data/fastq_trim/'
OUTDIR='/mnt/sde/renseb01/Documents/Steinberg_Christian/mirna/results/bams/'
THREADS=8

#for loop
for READS in $FASTQDIR*R1_trimmed.fq
  do
    SAMPLE=${READS//$FASTQDIR/}
    TEMP=${SAMPLE//'_R1_trimmed.fq'/''}
    PREFIX=${TEMP//'NX.VH01519_159.001.NEBNext_dual_i7_'/''}

    echo $READS
    echo $SAMPLE
    echo $TEMP
    echo $PREFIX
    echo "$OUTDIR${PREFIX}.sam"

    # Perform the alignment (for single-end reads)
    echo "Aligning single-end reads to reference genome..."
    #bwa mem -t $THREADS $REFERENCE_GENOME $READS > "$OUTDIR${PREFIX}.sam"
    bowtie -n 0 -l 15 -S $REFERENCE_INDEX $READS >$OUTDIR${PREFIX}.sam


    # Convert SAM to BAM (optional but recommended for downstream processing)
    echo "Converting SAM to BAM..."
    samtools view -Sb "$OUTDIR${PREFIX}.sam" > "$OUTDIR${PREFIX}.bam"

    # Sort the BAM file (optional)
    echo "Sorting BAM file..."
    samtools sort "$OUTDIR${PREFIX}.bam" -o "$OUTDIR${PREFIX}_sorted.bam"

    # Index the sorted BAM file (optional)
    echo "Indexing the sorted BAM file..."
    samtools index "$OUTDIR${PREFIX}_sorted.bam"


    # Count the data
    samtools idxstats "$OUTDIR${PREFIX}_sorted.bam" >"$OUTDIR${PREFIX}.idxstats"

    # Clean up: Remove the original SAM file if desired
    rm "$OUTDIR${PREFIX}.sam"
    rm "$OUTDIR${PREFIX}.bam"

    echo "Alignment completed. Output files:"
    echo "- ${PREFIX}_sorted.bam"
    echo "- ${PREFIX}_sorted.bam.bai"  # BAM index file
  done
