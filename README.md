# miRNA analysis pipeline  
* 1.Trimming the reads (`trim_galore` &  `cutadapt`) with  `trimming.sh` (single-end pilot) or `trimming_PE.sh` (paired-end REFRAME).  
  * Note that the REFRAME dataset are paired-end, the pilot was single end. We will only use the REFRAME dataset. 
* 2.Align the reads   
  * A. against the reference miRNA (`miRBase`) with Uracil (U) changed to Thymine (T) with **bowtie**. Then  **samtools idxstats** to count the reads per reference miRNA `align_miRNA_reference.R` for single-end pilot or `align_miRNA_reference_PE.R` for (REFRAME). OR
  * B. against the reference genome (`GRCh38`) with **bowtie**, then count the reads using a genome annotation file (`gencode.v43.basic.annotation.gtf`) and **htseq-count**. By aligning to the genome you can see what other kind of miRNA are present.
* 3.Plot results by running `mirna_expression.Rmd`  
  * A. Number of miRNA found and how many have an expression value.  
  * B. Histogram of the number of aligned sequences in cases vs. controls.  
  * C. Boxplot for the top 10 most expressed miRNA and how they compare between Case and control (number of reads per sample).   
* 4.Plot results by running `biomarker_discovery.Rmd`.  This produces a figure containing the results of several analyses.  
  * A. A Principal Component Analysis showing how the Case and Control group themselves.  
  * B. A ROC curve showing the result of the final LASSO logistic regression model. (60/40 training/testing, 100 resample to assess final gene model, but read the script to see how the final model was chosen).   
  * C. A plot showing what the final genes (biomarkers) were (i.e. based on 100 resamples, we choose the genes with a now zero median LASSO coefficient).  
  * D. A volcano plot to see the differential expression analysis. The final genes chosen in the LASSO logistic regression are labeled in blue.    


# Further information
  * sebastien.renaut.1@ulaval.ca
  * spring 2025
