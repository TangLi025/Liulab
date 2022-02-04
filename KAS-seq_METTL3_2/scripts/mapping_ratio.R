group <- "KAS-seq_ALKBH5"

summary <- read.table(paste0("~/KAS-METTL/",group,"/05_bedtools/bedGraph/bam_summary.txt"))

mean <- mean(summary[,1])

ratio <- mean(summary[,1])/summary[,1]


write(ratio,paste0("~/KAS-METTL/",group,"/05_bedtools/bedGraph/bam_ratio.txt"),sep=" ",ncolumns=8)
