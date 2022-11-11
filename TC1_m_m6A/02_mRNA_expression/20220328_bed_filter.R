library(rtracklayer)
library(Rsubread)
setwd("/disk1/home/user_08/Data/TC1-planB/10_mRNA/05_bam_hisat2/01_bam_sorted/")
#setwd("~/Data/TC1-planB/05_bam_change_index_hisat2/02_bam_dedup/")
bed_dir <- "06_bed_modify/" 
bed_file <- dir(bed_dir)[grep(".bed$",dir(bed_dir))]
bed <- paste0(bed_dir,bed_file)
bam_dir <- c("03_bam_neg/","02_bam_pos/")
bam_strand <- c("_neg.bam","_pos.bam")
dir.create("07_featurecount_input_read_in_peak")
dir.create("08_bed_filtered")

#fracOverlapFeature <- 0.8
fracOverlap <- 0.2

for (i in bed) 
{
  print(i)
  prefix <- strsplit(basename(i),"_")[[1]][c(1,2)]
  bed_df <- as.data.frame(import(i))
  saf <- data.frame(bed_df$name,bed_df[c(1,2,3,5)])
  colnames(saf)<-c("GeneID","Chr","Start","End","Strand")
  bam <- paste0(bam_dir,prefix[1],"_input_",prefix[2],bam_strand)
  print(bam)
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
  write.csv(count$counts,paste0("07_featurecount_input_read_in_peak/",prefix[1],"_",prefix[2],".csv"))
  print("fragments are")
  print(colSums(count$stat[,-1]))
  peak_to_retain <- NA
  for (j in 1:length(rownames(saf)))
  {
    if (saf$Strand[j] == "-")
    {
      peak_to_retain[j] <- count$counts[j,1] >= 1
    }
    else if (saf$Strand[j] == "+")
    {
      peak_to_retain[j] <- count$counts[j,2] >= 1
    }
  }
  export(bed_df[peak_to_retain,],paste0("08_bed_filtered/",basename(i)))
  print(paste0(length(peak_to_retain)," peaks, keep ",sum(peak_to_retain)," peaks"))
}
