summary <- read.table("~/KAS-METTL/A375/05_bedtools/bedGraph/bam_summary.txt")

mean <- mean(summary[,1])

ratio <- mean(summary[,1])/summary[,1]


write(ratio,"~/KAS-METTL/A375/05_bedtools/bedGraph/bam_ratio.txt",sep=" ",ncolumns=8)
