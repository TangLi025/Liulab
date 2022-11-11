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

GROUP="YTHDC1"

reads_density_TSS <- read.delim(paste0("~/KAS-seq_mouse/",GROUP,"/07_deeptools/reads_density/reads_density_TSS.tab"),check.names = FALSE)
summary(is.na(reads_density_TSS))
reads_density_TSS[is.na(reads_density_TSS)] <- 0
summary(is.na(reads_density_TSS))

reads_density_TSS_IP <- data.frame(CTRL_rep1=reads_density_TSS[,4],CTRL_rep2=reads_density_TSS[,5],KO_rep1=reads_density_TSS[,8],KO_rep2=reads_density_TSS[,9])

mean(reads_density_TSS_IP[,1])
mean(reads_density_TSS_IP[,2])
mean(reads_density_TSS_IP[,3])
mean(reads_density_TSS_IP[,4])

#reads_density_TSS_rep <- data.frame(CTRL_IP=(reads_density_TSS[,4]+reads_density_TSS[,5])/2,CTRL_input=(reads_density_TSS[,6]+reads_density_TSS[,7])/2)

reads_density_GB <- read.delim(paste0("~/KAS-seq_mouse/",GROUP,"/07_deeptools/reads_density/reads_density_GB.tab"),check.names = FALSE)
summary(is.na(reads_density_GB))
reads_density_GB[is.na(reads_density_GB)] <- 0
summary(is.na(reads_density_GB))

reads_density_GB_IP <- data.frame(CTRL_rep1=reads_density_GB[,4],CTRL_rep2=reads_density_GB[,5],KO_rep1=reads_density_GB[,8],KO_rep2=reads_density_GB[,9])

mean(reads_density_GB_IP[,1])
mean(reads_density_GB_IP[,2])
mean(reads_density_GB_IP[,3])
mean(reads_density_GB_IP[,4])

transcription_state <- matrix(nrow=4,ncol=4)
for (i in 1:4){
  transcription_state[1,i] <- as.numeric(summary(reads_density_TSS_IP[,i]>=5 & reads_density_GB_IP[,i]>=2.5)[3])
  
  transcription_state[2,i] <- as.numeric(summary(reads_density_TSS_IP[,i]>=5 & reads_density_GB_IP[,i]<2.5)[3])
  
  transcription_state[3,i] <- as.numeric(summary(reads_density_TSS_IP[,i]<5 & reads_density_GB_IP[,i]>=2.5)[3])
  transcription_state[4,i] <- as.numeric(summary(reads_density_TSS_IP[,i]<5 & reads_density_GB_IP[,i]<2.5)[3])
}

colnames(transcription_state) <- c("CTRL1","CTRL2","KO1","KO2")
rownames(transcription_state) <- c("class_1","class_2","class_3","class_4")

library(tidyverse)
library(reshape2)
transcription_state_t <- melt(transcription_state)
colnames(transcription_state_t) <- c("transcription_state","sample","gene_number")


library(ggplot2)

ggplot(data=transcription_state_t)+
  geom_col(mapping = aes(x=transcription_state,y=gene_number,fill=sample),position = 'dodge')+
  labs(title = GROUP)

