rm(list = ls())
options(stringsAsFactors = FALSE)
library("RADAR")

# To analyze for MeRIP-seq data, RADAR expect a pair of BAM file for each sample. 
# For each pair of BAM file, the naming convention for RADAR is sample.input.bam for INPUT 
# and sample.<methylation name>.bam for the IP, e.g. sample.m6A.bam for m6A IP sample. 
# They key is INPUT and IP sample have the same prefix.

# countReads() function to :
#1) concatenate exons of each gene to get a “longest isoform” transcript of each gene.
#2) Divide transcripts into consecutive bins of user defined width. 
#3) quantify reads mapped to each bin.

# The default setting for bin width to slice the transcript is 50bp. 
# One can set this by parameter binSize = #bp. 
# For very shallow sequencing depth (e.g. less than 10M mappable reads per library), 
# we recommand setting bin width to larger size such as 100bp to increase the number of reads countable in each bin.

# strandToKeep: According to library preparation protocol, choose which strand to count. 
# Stranded RNA library usually seq the "ooposite" strand. Small RNA library seq the "same" strand.

# The parameter fragmentLength defines the number of nucleotide to shift from the head of a read in BAM record to the center of RNA fragment.
# In the example data, the RNA has been sonicated into ~150nt fragment before IP and library preparation.

samplename <- c("Lysate.rep1","Lysate.rep2","Result.rep1","Result.rep2")
setwd("~/LinLong/20211230")
lovo_radar <- countReads(
  samplenames = samplename[1:4],  # prefix for bam files
  gtf = "/disk1/home/user_09/reference/annotation/mm9/gencode.vM1.annotation.gtf",
  bamFolder = "./04_bam_dedup",
  modification = "IP", 
  strandToKeep = "opposite",
  threads = 50,
  binSize = 100,
  outputDir='./16_diff_peaks/dedup',
  fragmentLength = 150,
  paired = TRUE,
  saveOutput = TRUE
)
hydrogel=lovo_radar
summary(lovo_radar)

# library normalization step
# use top read count bins of IP and corresponding input gene-level read count to compute the estimated enrichment of each sample. 
# This procedure is under the assumption that samples in the same study have the same IP efficiency. 
# We then normalize the IP read counts by the estimated IP efficiency. 
# result: There will be normalized gene-level INPUT count matrix stored in this object, 
# which is essentially RNA-seq read count matrix and can be accessed by geneExression(radar)
lovo_radar <- normalizeLibrary( lovo_radar )
sizeFactors(lovo_radar)

# we can use the geneSum of INPUT to adjust for the variation of expression level of the IP read count. 
# This step aims to account for the variation of IP read count attributed to variation of pre-IP gene expression level.
lovo_radar <- adjustExprLevel( lovo_radar )

# we will need to filter out bins of very low read count 
# because under sampled locus can be strongly affected by technical variation.
# This step will also filter out bins where IP has less coverage than Input 
# because we only care about loci where m6A is enrichment.
variable(lovo_radar) <- data.frame( group = factor(c("Ctl","Ctl","Treated","Treated") ))
lovo_radar <- filterBins( lovo_radar ,minCountsCutOff = 15)

# Now we have the pre-processed read counts matrix for testing differential methylation.
# lovo_radar <- diffIP(lovo_radar)
# if using linux or mac, you can also use multi-thread mode
lovo_radar <- diffIP_parallel(lovo_radar, thread = 40)
head(lovo_radar@test.est)

# The code above run test on the effect of predictor variable on the methylation level, 
# which is fine if all covariates are well balanced. 
# However, in most cases, samples might have been processed in more than one batches (especially when sample size is large)
# where batch effect and other covariate could contributed to the variation.
# In order to check for unwanted variation, 
# we can take the top 1000 bins ranked by count number (basically using the high read count bins) to plot PCA:
top_bins <- extractIP(lovo_radar,filtered = T)[order(rowMeans( extractIP(lovo_radar,filtered = T) ),decreasing = T)[1:1000],]
plotPCAfromMatrix(top_bins,group = unlist(variable(lovo_radar)) )

# Due to the resolution of the MeRIP-seq experiment where RNA molecules are fragmented into 100-300nt
# neighboring bins can usually contain reads from the same locus. 
# Therefore, we do a post-processing to merge significant neighboring bins after the test to obtain a final list of differential peaks.
# We merge the p-value of connecting bins by fisher’s method and report the max beta from neighbouring bins.
# Here, we use FDR<0.1 and log fold change > 0.5 as default cutoff for selecting significant bins. 
res <- reportResult( lovo_radar, cutoff = 0.5, Beta_cutoff = 0, threads = 40 ) # 5分半

result <- results(res)
plotHeatMap(res)



# we can extract the final result as an table (in BED12 format) 
# for the genomic location of the differential methylated peaks with p-values and log fold changes.
lovo_res <- results( res)
head(lovo_res)
length(unique(lovo_res$name))
tmp=lovo_res[abs(lovo_res$logFC)>=1,]
length(unique(tmp$name))
lovo.m6a=tmp
write.csv(lovo.m6a,file='./16_diff_peaks/hytrogel_diff_peaks_log2fc1_fdr0.35.csv',quote = F,row.names = F)
