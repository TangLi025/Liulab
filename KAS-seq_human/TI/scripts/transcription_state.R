rm(list=ls())
library(rtracklayer)

gene <- read.delim("~/reference/annotation/hg19/gencode.v19.annotation.protein_coding.chr.gene.bed",header = F)
gene <- gene[(gene[,3]-gene[,2]) >500,]

TSS_bed <- gene
TSS_bed[TSS_bed[,6]=="+",3] <- TSS_bed[TSS_bed[,6]=="+",2]+400
TSS_bed[TSS_bed[,6]=="+",2] <- TSS_bed[TSS_bed[,6]=="+",2]-200

TSS_bed[TSS_bed[,6]=="-",2] <- TSS_bed[TSS_bed[,6]=="-",3]-400
TSS_bed[TSS_bed[,6]=="-",3] <- TSS_bed[TSS_bed[,6]=="-",3]+200
TSS_bed[,2] <- as.integer(TSS_bed[,2])
TSS_bed[,3] <- as.integer(TSS_bed[,3])

write.table(TSS_bed[,1:6],file="~/reference/annotation/hg19/gencode.v19.annotation.protein_coding.chr.gene.TSS.bed",sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

GB_bed <- gene
GB_bed[GB_bed[,6]=="+",2] <- GB_bed[GB_bed[,6]=="+",2]+400

GB_bed[GB_bed[,6]=="-",3] <- GB_bed[GB_bed[,6]=="-",3]-400
GB_bed[,2] <- as.integer(GB_bed[,2])
GB_bed[,3] <- as.integer(GB_bed[,3])
write.table(GB_bed[,1:6],file="~/reference/annotation/hg19/gencode.v19.annotation.protein_coding.chr.gene.GB.bed",sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

TSS <- read.delim("~/KAS/TI_KAS-seq/07_deeptools/computeMatrix/DMSO_IP_merge_TSS.tab", header=T, skip=2)
TSS <- TSS[,-13]
gene_body_10000 <- read.delim("~/KAS/TI_KAS-seq/07_deeptools/computeMatrix/DMSO_IP_merge_body_10000.tab", header=T, skip=2)
gene_body_10000 <- gene_body_10000[,-201]
gene_body_1000 <- read.delim("~/KAS/TI_KAS-seq/07_deeptools/computeMatrix/DMSO_IP_merge_body_1000.tab", header=T, skip=2)
gene_body_1000 <- gene_body_1000[,-21]

density_mean <- data.frame(gene[,4],apply(TSS,1,mean),apply(gene_body_10000,1,mean),apply(gene_body_1000,1,mean))

TSS_mean <- mean(density_mean[,2])
GB_mean <- mean(density_mean[,3])

bw <- import("~/KAS/TI_KAS-seq/07_deeptools/bamCoverage/DMSO_IP_merge.bw")
mean_total <- 55262849*3/2827437033


class1_gene <- density_mean[density_mean[,2]>=mean_total*20 & density_mean[,3]>=mean_total*10,]
class2_gene <- density_mean[density_mean[,2]>=mean_total*20 & density_mean[,3]<mean_total*10,]
class3_gene <- density_mean[density_mean[,2]<mean_total*20 & density_mean[,3]>=mean_total*10,]
class4_gene <- density_mean[density_mean[,2]<mean_total*20 & density_mean[,3]<mean_total*10,]

class1_gene <- density_mean[density_mean[,2]>=TSS_mean*1 & density_mean[,3]>=GB_mean*1,]
class2_gene <- density_mean[density_mean[,2]>=TSS_mean*1 & density_mean[,3]<GB_mean*1,]
class3_gene <- density_mean[density_mean[,2]<TSS_mean*1 & density_mean[,3]>=GB_mean*1,]
class4_gene <- density_mean[density_mean[,2]<TSS_mean*1 & density_mean[,3]<GB_mean*1,]



