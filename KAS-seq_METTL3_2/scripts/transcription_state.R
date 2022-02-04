rm(list=ls())
library(rtracklayer)

gene <- read.delim("~/reference/annotation/mm19/mm19_Refseq.bed",header = F)
gene <- gene[(gene[,3]-gene[,2]) >500,]

TSS_bed <- gene
TSS_bed[TSS_bed[,6]=="+",3] <- TSS_bed[TSS_bed[,6]=="+",2]+400
TSS_bed[TSS_bed[,6]=="+",2] <- TSS_bed[TSS_bed[,6]=="+",2]-200

TSS_bed[TSS_bed[,6]=="-",2] <- TSS_bed[TSS_bed[,6]=="-",3]-400
TSS_bed[TSS_bed[,6]=="-",3] <- TSS_bed[TSS_bed[,6]=="-",3]+200
TSS_bed[,2] <- as.integer(TSS_bed[,2])
TSS_bed[,3] <- as.integer(TSS_bed[,3])

write.table(TSS_bed[,1:6],file="~/reference/annotation/mm19/mm19_Refseq.TSS.bed",sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

GB_bed <- gene
GB_bed[GB_bed[,6]=="+",2] <- GB_bed[GB_bed[,6]=="+",2]+400

GB_bed[GB_bed[,6]=="-",3] <- GB_bed[GB_bed[,6]=="-",3]-400
GB_bed[,2] <- as.integer(GB_bed[,2])
GB_bed[,3] <- as.integer(GB_bed[,3])
write.table(GB_bed[,1:6],file="~/reference/annotation/mm19/mm19_Refseq.GB.bed",sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)



#### 使用multiBigWigSummary BED-file计算reads density

reads_density <- read.delim("~/KAS-seq_mouse/METTL3_3/07_deeptools/reads_density/reads_density_TSS.tab",check.names = FALSE)
summary(is.na(reads_density))
reads_density[is.na(reads_density)] <- 0
reads_density_mean <- mean(reads_density[,6])

summary(reads_density[,4] > reads_density_mean*20)

summary <- read.table("~/KAS-METTL/METTL3_3/05_bedtools/bedGraph/bam_summary.txt")

mean_reads <- mean(summary[,1])

TSS <- read.delim("~/KAS-METTL/METTL3_3/07_deeptools/computeMatrix/KO_input_rep1_TSS.tab", header=T, skip=2)
TSS <- TSS[,-13]
gene_body <- read.delim("~/KAS-METTL/METTL3_3/07_deeptools/computeMatrix/KO_input_rep1_body.tab", header=T, skip=2)
gene_body <- gene_body[,-21]

density_mean <- data.frame(apply(TSS,1,mean),apply(gene_body,1,mean))

TSS_mean <- mean(density_mean[,1])

GB_mean <- mean(density_mean[,2])

bw <- import("~/KAS-METTL/METTL3_3/05_bedtools/bigWig/KAS-seq_METTL3_3_KO_input_rep1_ext.bw")


mean_total <- sum(bw$score*bw@ranges@width)/2494787188/50

class1_gene <- density_mean[density_mean[,1]>=mean_total*20 & density_mean[,2]>=mean_total*10,]
class2_gene <- density_mean[density_mean[,1]>=mean_total*20 & density_mean[,2]<mean_total*10,]
class3_gene <- density_mean[density_mean[,1]<mean_total*20 & density_mean[,2]>=mean_total*10,]
class4_gene <- density_mean[density_mean[,1]<mean_total*20 & density_mean[,2]<mean_total*10,]








