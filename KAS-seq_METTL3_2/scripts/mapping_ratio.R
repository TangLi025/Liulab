summary <- read.table("~/KAS-METTL/METTL3_2/04_bam_rmdup/bam_summary.txt")

mean <- mean(summary[,1])

ratio <- mean(summary[,1])/summary[,1]


write(ratio,"~/KAS-METTL/METTL3_2/04_bam_rmdup/bam_ratio.txt",sep="\n")
