GROUP="METTL3_3"


CTRL_bed1 <- read.delim("/disk1/home/user_09/carRNA_science/m6A-seq_narrowPeak/Mettl3_Control_commen_peaks_annotation_homer.xls",header = T,)

m6A_plus_name <- gtf$gene_name[gtf$gene_name %in% m6A_gene_name]

library(rtracklayer)
gtf=import('/disk1/home/user_09/KAS/fpkm_m6A_Liu/fpkms_high.gene.gtf')

m6A_gene_name <- CTRL_bed1$Gene.Name

m6A_plus <- gtf[gtf$gene_name %in% m6A_gene_name,]
m6A_minus <- gtf[!(gtf$gene_name %in% m6A_gene_name),]

export(m6A_plus,)

refseq <- read.delim(file="~/reference/annotation/mm19/mm19_Refseq.bed", header = FALSE)

refseq_high <- refseq[refseq$V4 %in% gtf$gene_name,]

summary(refseq[,4] %in% m6A_gene_name)

m6A_plus <- refseq_high[refseq_high[,4] %in% m6A_gene_name,]
m6A_minus <- refseq_high[!(refseq_high[,4] %in% m6A_gene_name),]

write.table(m6A_plus,file="/disk1/home/user_09/carRNA_science/m6A-seq_narrowPeak/m6A_plus.bed",quote = FALSE,sep='\t',row.names = FALSE,col.names = FALSE)
write.table(m6A_minus,file="/disk1/home/user_09/carRNA_science/m6A-seq_narrowPeak/m6A_minus.bed",quote = FALSE,sep='\t',row.names = FALSE,col.names = FALSE)
