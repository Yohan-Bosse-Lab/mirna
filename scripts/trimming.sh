#!/bin/bash

conda activate trimming

#dir
fastqdir='/home/renseb01/Documents/lord/raw_data/REFRAME_miRNA/pilot_2025/'
fastq_trim='/mnt/sde/renseb01/Documents/Steinberg_Christian/mirna/data/fastq_trim3/'

#for loop
for R1 in $fastqdir/*CASE*
  do
    #prepare files
    nameR1=${R1//$fastqdir/}
    trimgal_R1=${nameR1//'_R1.fastq.gz'/'_R1_val_1.fq.gz'}
    cutadapt_R1=${nameR1//'_R1.fastq.gz'/'_cutadapt_R1.fq.gz'}

    echo $R1
    echo $nameR1
    echo $trimgal_R1
    echo $cutadapt_R1

    #trim_galore
    trim_galore --cores 2 \
                --length 16 \
                -q 20 \
                -o $fastq_trim \
                --phred33 \
                $fastqdir$nameR1 \
                --path_to_cutadapt cutadapt

    #remove the specific polyG sequences
#    cutadapt -a GGGGGGGGGGGGGGG \
#             -A GGGGGGGGGGGGGGG \
#             --cores 2 \
#             -o $fastq_trim$cutadapt_R1 \
#             $fastq_trim$trimgal_R1
  done
