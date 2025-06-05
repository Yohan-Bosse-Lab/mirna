#!/bin/bash

#dir
fastqdir='/home/renseb01/Documents/lord/raw_data/REFRAME_miRNA/REFRAME_2025/fastq_files/'
fastq_trim='/mnt/sde/renseb01/Documents/Steinberg_Christian/mirna/data/fastq_trim_REFRAME_2025/'

#for loop
for R1 in $fastqdir/*_R1.fastq.gz
  do
    nameR1=${R1//$fastqdir/}
    trimgal_R1=${nameR1//'_R1.fastq.gz'/'_R1_val_1.fq.gz'}
    fastp_R1=${trimgal_R1//'_R1_val_1.fq.gz'/'_trimmed_R1.fq.gz'}

    R2=${R1//'R1'/'R2'}
    nameR2=${R2//$fastqdir/}
    trimgal_R2=${nameR2//'_R2.fastq.gz'/'_R2_val_2.fq.gz'}
    fastp_R2=${trimgal_R2//'_R2_val_2.fq.gz'/'_trimmed_R2.fq.gz'}

    echo '1st pair'
    echo $nameR1

    echo '2nd pair'
    echo $nameR2

    #trim_galore
    trim_galore --cores 6 \
                --length 16 \
                --fastqc \
                -q 20 \
                -o $fastq_trim \
                --phred33 \
                --paired $fastqdir$nameR1 $fastqdir$nameR2 \
                --path_to_cutadapt cutadapt


    fastp -i $fastq_trim$trimgal_R1 \
          -I $fastq_trim$trimgal_R2 \
          -o $fastq_trim$fastp_R1 \
          -O $fastq_trim$fastp_R2 \
          --qualified_quality_phred 20 \
          --unqualified_percent_limit 0

  done

rm $fastq_trim*val_1*
rm $fastq_trim*val_2*
