log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")


library(rtracklayer)
library(Rsubread)
setwd("~/LinLong/")

fracOverlapFeature <- 0.8
fracOverlap <- 0.2

bed <- snakemake@input[[1]]

bed_df <- as.data.frame(rtracklayer::import(bed))
saf <- data.frame(bed_df$name,bed_df[c(1,2,3,5)])
colnames(saf)<-c("GeneID","Chr","Start","End","Strand")

bam <- c(snakemake@input[[2]],snakemake@input[[3]])
count <- featureCounts(files = bam, 
                       annot.ext = saf,
                       isGTFAnnotationFile = F,
                       minMQS = 20, 
                       strandSpecific = 2,
                       countMultiMappingReads = FALSE,
                       isPairedEnd = TRUE,
                       #fracOverlapFeature = fracOverlapFeature,
                       fracOverlap = fracOverlap,
                       maxFragLength = 2000,                         
                       nthreads = 20)
write.csv(count$counts,snakemake@output[[1]])
print("fragments are")
print(colSums(count$stat[,-1]))
peak_to_retain <- NA
print(length(rownames(saf)))

for (j in 1:length(rownames(saf)))
{{
  if (saf$Strand[j]== "-")
  {
    peak_to_retain[j] <- count$counts[j,1] >= 1
  }
  else if (saf$Strand[j]== "+")
  {
    peak_to_retain[j] <- count$counts[j,2] >= 1
  }
}}

export(bed_df[peak_to_retain,],snakemake@output[[2]])
print(paste0(length(peak_to_retain)," peaks, keep ",sum(peak_to_retain)," peaks"))


