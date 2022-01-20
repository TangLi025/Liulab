
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("exomePeak2")

library(exomePeak2)
library(BSgenome)

set.seed(1)

GENE_ANNO_GTF = "/disk1/home/user_09/reference/annotation/mm10/gencode.vM19.annotation.gtf"

f1 = "/disk1/home/user_09/LinLong/04_bam_rmdup/Lysate_IP_rep1.bam"
f2 = "/disk1/home/user_09/LinLong/04_bam_rmdup/Lysate_IP_rep2.bam"
f3 = "/disk1/home/user_09/LinLong/04_bam_rmdup/Result_IP_rep1.bam"
f4 = "/disk1/home/user_09/LinLong/04_bam_rmdup/Result_IP_rep2.bam"
IP_BAM = c(f1,f2)
TREATED_IP_BAM = c(f3,f4)
#IP_BAM = c(f1,f2,f3,f4)

f1 = "/disk1/home/user_09/LinLong/04_bam_rmdup/Lysate_input_rep1.bam"
f2 = "/disk1/home/user_09/LinLong/04_bam_rmdup/Lysate_input_rep2.bam"
f3 = "/disk1/home/user_09/LinLong/04_bam_rmdup/Result_input_rep1.bam"
f4 = "/disk1/home/user_09/LinLong/04_bam_rmdup/Result_input_rep2.bam"
INPUT_BAM = c(f1,f2)
TREATED_INPUT_BAM = c(f3,f4)
#INPUT_BAM = c(f1,f2,f3,f4)

txdb =makeTxDbFromGFF("/disk1/home/user_09/reference/annotation/mm10/gencode.vM19.annotation.gtf")
bsgenome = getBSgenome("mm10")

### 6 Multi-step Functions
# The exomePeak2 package can achieve peak calling and peak statistics calculation with multiple functions.

# 1. Check the BAM files of MeRIP-seq experiments

# The scanMeripBAM() is used to organize the BAM files used in a MeRIP-Seq analysis, it also specify the important BAM reading parameters such as the paired end library type and the FLAG filters.

MeRIP_Seq_Alignment <- scanMeripBAM(
  bam_ip = IP_BAM,
  bam_input = INPUT_BAM,
  paired_end = TRUE
)

# For MeRIP-seq experiment with interactive design (contain control and treatment groups), we can add the treated input and IP files with the arguments bam_treated_input and bam_treated_ip:
  
MeRIP_Seq_Alignment <- scanMeripBAM(
    bam_ip = IP_BAM,
    bam_input = INPUT_BAM,
    bam_treated_input = TREATED_INPUT_BAM,
    bam_treated_ip = TREATED_IP_BAM,
    paired_end = TRUE
  ) 

# 2. Perform peak calling on exons

# After checking the MeRIP-seq BAMs, peak calling can be performed with exomePeakCalling(). The function will first generate sliding windows on exons of the annotation file, and it will then count the read ans it GLMs to detect the significantly modified peaks. If either genome or bsgenome arguments are provided, the GC content bias correction will be performed while peak calling.

SummarizedExomePeaks <- exomePeakCalling(merip_bams = MeRIP_Seq_Alignment,
                                         txdb = txdb,
                                         bsgenome = bsgenome,
                                         fragment_length = 150,
                                         parallel = 30) 

# 3. Estimate correction factors

SummarizedExomePeaks <- estimateSeqDepth(SummarizedExomePeaks)
SummarizedExomePeaks <- normalizeGC(SummarizedExomePeaks)

# P.S. The GC content bias factors can still be estimated by normalizeGC() if no genome sequence information are provided in the previous steps, but genome information need to be provided at normalizeGC() via the bsgenome argument.

# In addition, the sequencing depth and GC content offset can be calculated according to 3 types of feature ranges: “All” (default), “Background” and “Modification”. The impact of different estimation scopes are highly data dependent, and reasonable changes of these may significantly improve the biological outcomes of the downstream analysis. For example, in differential analysis, it is recommended to use estimateSeqDepth(SummarizedExomePeaks, from='Background') so that the direction of differential modification directions can be more consistent with the expectation of the perturbed protein factors.

### 4. Report the GLM statistics

# Using the normalization factors calculated above, the GLM statistics of the IP/input design can be calculated with the function glmM():
  
  SummarizedExomePeaks <- glmM(SummarizedExomePeaks) 

# If the treated IP and input BAM are provided, glmDM() function can be used to conduct differential modification analysis on modification Peaks with interactive GLM:
  
  SummarizedExomePeaks <- glmDM(SummarizedExomePeaks)

# All of the GLM statistics and normalization factors calculated will be contained in the SummarizedExomePeaks object, which can be saved (ex. with saveRDS()) to store the intermediate results in a given pipeline of analysis.

### 5. Scatter plot between GC content and log2 fold changes (LFCs).

# A scatter plot can be produced with plotLfcGC() to visualize the trend of LFCs with GC. It is recommended to plot it twice before and after conducting the GC content normalization using normalizeGC().

plotLfcGC(SummarizedExomePeaks, point_size = 1, xlim = c(0.4, 0.85)) 


### 6. Bar plot of sequencing depth estimates

# The library size factors can be visualized and compared using plotSizeFactors(). The label will allow comparison of 3 types of feature ranges for size factor calculations, including the background, modification and both regions (All).

plotSizeFactors(SummarizedExomePeaks)


### 7. Export the modification peaks and their statistics

# The analysis results can be saved on the system in formats of CSV, BED, or RDS.

exportResults(SummarizedExomePeaks, format = "BED") 

# Please note that exportResults() can also specify the filtering thresholds and style of the result table.

### 7 Contact
# If you encounter any problems during use, please contact the maintainer of exomePeak2:
  
#  Zhen Wei : zhen.wei01@xjtlu.edu.cn

### 8 References
# KD Hansen, RA Irizarry, and Z Wu, Removing technical variability in RNA-seq data using conditional quantile normalization. Biostatistics 2012 vol. 13(2) pp. 204-216.

# Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. https://doi.org/10.1186/s13059-014-0550-8

# Zhu A, Ibrahim JG, Love MI (2018). “Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences.” Bioinformatics. doi: 10.1093/bioinformatics/bty895.

### 9 Session Info
# sessionInfo()
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.3 LTS
## 
## Matrix products: default
## BLAS:   /home/biocbuild/bbs-3.14-bioc/R/lib/libRblas.so
## LAPACK: /home/biocbuild/bbs-3.14-bioc/R/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB              LC_COLLATE=C              
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] splines   stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] BSgenome.Hsapiens.UCSC.hg19_1.4.3 BSgenome_1.62.0                  
##  [3] rtracklayer_1.54.0                Biostrings_2.62.0                
##  [5] XVector_0.34.0                    exomePeak2_1.6.1                 
##  [7] cqn_1.40.0                        quantreg_5.86                    
##  [9] SparseM_1.81                      preprocessCore_1.56.0            
## [11] nor1mix_1.3-0                     mclust_5.4.8                     
## [13] SummarizedExperiment_1.24.0       Biobase_2.54.0                   
## [15] GenomicRanges_1.46.1              GenomeInfoDb_1.30.0              
## [17] IRanges_2.28.0                    S4Vectors_0.32.3                 
## [19] BiocGenerics_0.40.0               MatrixGenerics_1.6.0             
## [21] matrixStats_0.61.0                BiocStyle_2.22.0                 
## 
## loaded via a namespace (and not attached):
##   [1] systemfonts_1.0.3        BiocFileCache_2.2.0      plyr_1.8.6              
##   [4] BiocParallel_1.28.3      listenv_0.8.0            ggplot2_3.3.5           
##   [7] digest_0.6.29            foreach_1.5.1            htmltools_0.5.2         
##  [10] magick_2.7.3             fansi_0.5.0              magrittr_2.0.1          
##  [13] memoise_2.0.1            recipes_0.1.17           globals_0.14.0          
##  [16] annotate_1.72.0          gower_0.2.2              bdsmatrix_1.3-4         
##  [19] prettyunits_1.1.1        colorspace_2.0-2         apeglm_1.16.0           
##  [22] rappdirs_0.3.3           blob_1.2.2               textshaping_0.3.6       
##  [25] xfun_0.28                dplyr_1.0.7              crayon_1.4.2            
##  [28] RCurl_1.98-1.5           jsonlite_1.7.2           genefilter_1.76.0       
##  [31] survival_3.2-13          iterators_1.0.13         glue_1.5.1              
##  [34] gtable_0.3.0             ipred_0.9-12             zlibbioc_1.40.0         
##  [37] MatrixModels_0.5-0       DelayedArray_0.20.0      future.apply_1.8.1      
##  [40] scales_1.1.1             mvtnorm_1.1-3            DBI_1.1.1               
##  [43] Rcpp_1.0.7               emdbook_1.3.12           xtable_1.8-4            
##  [46] progress_1.2.2           bit_4.0.4                lava_1.6.10             
##  [49] prodlim_2019.11.13       httr_1.4.2               RColorBrewer_1.1-2      
##  [52] ellipsis_0.3.2           farver_2.1.0             pkgconfig_2.0.3         
##  [55] XML_3.99-0.8             nnet_7.3-16              sass_0.4.0              
##  [58] dbplyr_2.1.1             locfit_1.5-9.4           utf8_1.2.2              
##  [61] caret_6.0-90             labeling_0.4.2           tidyselect_1.1.1        
##  [64] rlang_0.4.12             reshape2_1.4.4           AnnotationDbi_1.56.2    
##  [67] munsell_0.5.0            tools_4.1.2              cachem_1.0.6            
##  [70] generics_0.1.1           RSQLite_2.2.9            evaluate_0.14           
##  [73] stringr_1.4.0            fastmap_1.1.0            ragg_1.2.1              
##  [76] yaml_2.2.1               ModelMetrics_1.2.2.2     knitr_1.36              
##  [79] bit64_4.0.5              purrr_0.3.4              KEGGREST_1.34.0         
##  [82] future_1.23.0            nlme_3.1-153             xml2_1.3.3              
##  [85] biomaRt_2.50.1           compiler_4.1.2           filelock_1.0.2          
##  [88] curl_4.3.2               png_0.1-7                tibble_3.1.6            
##  [91] geneplotter_1.72.0       bslib_0.3.1              stringi_1.7.6           
##  [94] highr_0.9                GenomicFeatures_1.46.1   lattice_0.20-45         
##  [97] Matrix_1.4-0             vctrs_0.3.8              pillar_1.6.4            
## [100] lifecycle_1.0.1          BiocManager_1.30.16      jquerylib_0.1.4         
## [103] data.table_1.14.2        bitops_1.0-7             conquer_1.2.1           
## [106] R6_2.5.1                 BiocIO_1.4.0             bookdown_0.24           
## [109] parallelly_1.29.0        codetools_0.2-18         MASS_7.3-54             
## [112] assertthat_0.2.1         DESeq2_1.34.0            rjson_0.2.20            
## [115] withr_2.4.3              GenomicAlignments_1.30.0 Rsamtools_2.10.0        
## [118] GenomeInfoDbData_1.2.7   parallel_4.1.2           hms_1.1.1               
## [121] grid_4.1.2               rpart_4.1-15             timeDate_3043.102       
## [124] coda_0.19-4              class_7.3-19             rmarkdown_2.11          
## [127] bbmle_1.0.24             pROC_1.18.0              numDeriv_2016.8-1.1     
## [130] lubridate_1.8.0          restfulr_0.0.13