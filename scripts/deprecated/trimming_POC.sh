#!/bin/bash

#dir
fastqdir='/home/renseb01/Documents/lord/raw_data/REFRAME_miRNA/REFRAME_2025/fastq_files/'
fastq_trim='/mnt/sde/renseb01/Documents/Steinberg_Christian/mirna/data/fastq_trim2/'

#for loop
for R1 in $fastq_trim/*val_1.fq.gz
  do
    nameR1=${R1//$fastq_trim/}
    trimgal_R1=$nameR1
    cutadapt_R1=${nameR1//'_R1_val_1.fq.gz'/'_cutadapt_R1.fq.gz'}

    #R2=${R1//'R1'/'R2'}
    #nameR2=${R2//$fastq_trim/}
    trimgal_R2=${nameR1//'_R1_val_1.fq.gz'/'_R2_val_2.fq.gz'}
    cutadapt_R2=${nameR1//'_R1_val_1.fq.gz'/'_cutadapt_R2.fq.gz'}

#    echo $R1
 #   echo $nameR1
    echo '1st pair'
    echo $trimgal_R1
    echo $cutadapt_R1

  #  echo $R2
   # echo $nameR2
    echo '2nd pair'
    echo $trimgal_R2
    echo $cutadapt_R2

    #trim_galore
#    trim_galore --cores 2 \
 #               --length 16 \
  #              -q 20 \
   #             -o $fastq_trim \
    #            --phred33 \
     #           --paired $fastqdir$nameR1 $fastqdir$nameR2 \
      #          --path_to_cutadapt cutadapt

    #remove the specific polyG sequences
    cutadapt -a GGGGGGGGGG \
             -A GGGGGGGGGG \
             --cores 2 \
             -o $fastq_trim$cutadapt_R1 \
             -p $fastq_trim$cutadapt_R2 $fastq_trim$trimgal_R1 $fastq_trim$trimgal_R2
  done
