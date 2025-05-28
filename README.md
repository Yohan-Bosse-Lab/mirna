# mirna
  * miRNA analysis pipeline
  * 1. trimming the reads (`trim_galore` &  `cutadapt`)
  * 2. Align against the reference miRNA (`miRBase`) with Urical changed to Thymine with **bowtie**. Then  **samtools idxstats** to count the reads per reference miRNA. OR
  * 3. Align against the reference genome (`GRCh38`) with **bowtie**, then count the reads using a genome annotation file (`gencode.v43.basic.annotation.gtf`) and **htseq-count**
  * 4. Plot results by running `mirna_expression.Rmd`


# Further information
  * sebastien.renaut.1@ulaval.ca
  * spring 2025
