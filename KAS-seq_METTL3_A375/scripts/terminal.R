library(rtracklayer)

gene <- read.delim("~/reference/annotation/hg19/gencode.v19.annotation.protein_coding.chr.gene.bed",header = F)
gene_1 <- rbind(gene[1,],gene)
gene_0 <- rbind(gene,gene[1,])

gene[gene[,6]=="+",][gene[gene[,6]=="+",3]+10000 < gene[which(gene[,6]=="+")+1,2] | gene[which(gene[,6]=="+"),1]!=gene[which(gene[,6]=="+")+1,1],]
gene[gene[,6]=="-",][gene[gene[,6]=="-",2]-10000 > gene[which(gene[,6]=="-")-1,3] | gene[which(gene[,6]=="-"),1]!=gene[which(gene[,6]=="-")-1,1],]
gene_terminal <- rbind(gene[gene[,6]=="+",][gene[gene[,6]=="+",3]+10000 < gene[which(gene[,6]=="+")+1,2] | gene[which(gene[,6]=="+"),1]!=gene[which(gene[,6]=="+")+1,1],],gene[gene[,6]=="-",][gene[gene[,6]=="-",2]-10000 > gene[which(gene[,6]=="-")-1,3] | gene[which(gene[,6]=="-"),1]!=gene[which(gene[,6]=="-")-1,1],],gene[20332,])


TERMINAL_bed <- gene_terminal
TERMINAL_bed[which(TERMINAL_bed[,6]=="+"),2] <- TERMINAL_bed[which(TERMINAL_bed[,6]=="+"),3]
TERMINAL_bed[which(TERMINAL_bed[,6]=="+"),3] <- TERMINAL_bed[which(TERMINAL_bed[,6]=="+"),3]+10000

TERMINAL_bed[which(TERMINAL_bed[,6]=="-"),3] <- TERMINAL_bed[which(TERMINAL_bed[,6]=="-"),2]
TERMINAL_bed[which(TERMINAL_bed[,6]=="-"),2] <- TERMINAL_bed[which(TERMINAL_bed[,6]=="-"),2]-10000
TERMINAL_bed[which(TERMINAL_bed[,2]<0),2] <- 1
TERMINAL_bed[,2] <- as.integer(TERMINAL_bed[,2])
TERMINAL_bed[,3] <- as.integer(TERMINAL_bed[,3])
TERMINAL_bed <- na.omit(TERMINAL_bed)

write.table(TERMINAL_bed[,1:6],file="~/reference/annotation/hg19/gencode.v19.annotation.protein_coding.chr.gene.terminal.bed",sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

terminal <- read.delim("~/KAS/TI_KAS-seq/07_deeptools/computeMatrix/DMSO_IP_merge_terminal.tab", header=T, skip=2)
terminal <- terminal[,-501]

terminal_pos <- terminal >= 5

terminal_rank <- data.frame(TERMINAL_bed[,4],rowSums(terminal_pos))
terminal_rank <- terminal_rank[order(terminal_rank[,2],decreasing=TRUE),]

length <- dim(terminal_rank)[1]
long <- terminal_rank[1:round(length/3),]
medium <- terminal_rank[(round(length/3)+1):round(length*2/3),]
short <- terminal_rank[(round(length*2/3)+1):length,]

library(rtracklayer)
gtf=import('/disk1/home/user_09/reference/annotation/hg19/gencode.v19.annotation.protein_coding.chr.gtf')

long_gtf=gtf[gtf$gene_id %in% long[,1]]
export(long_gtf, "~/KAS/TI_KAS-seq/long.gtf")

medium_gtf=gtf[gtf$gene_id %in% medium[,1]]
export(medium_gtf, "~/KAS/TI_KAS-seq/medium.gtf")

short_gtf=gtf[gtf$gene_id %in% short[,1]]
export(short_gtf, "~/KAS/TI_KAS-seq/short.gtf")

