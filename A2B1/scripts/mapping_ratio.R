group <- "dedup"

summary <- read.table(paste0("~/LinLong/05_bedtools/",group,"/bedGraph/bam_summary.txt"))

mean <- mean(summary[,1])

ratio <- mean(summary[,1])/summary[,1]

ratio

write(ratio,paste0("~/LinLong/05_bedtools/",group,"/bedGraph/bam_ratio.txt"),sep=" ",ncolumns=8)
