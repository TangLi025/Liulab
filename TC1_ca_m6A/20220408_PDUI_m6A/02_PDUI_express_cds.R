## 把featureCounts换成cds,避免3UTR长度的影响
library(Rsubread)

setwd("/disk/user_08/Data/TC1-planB/05_bam_change_index_hisat2/01_bam_sorted/")
path <- "/disk/user_08/Data/TC1-planB/05_bam_change_index_hisat2/01_bam_sorted/11_bam_merge/"
sample_name <- c("p0_input_rep1","p0_input_rep2","p5_input_rep1","p5_input_rep2",
                 "p10_input_rep1","p10_input_rep2","rp2_input_rep1","rp2_input_rep2",
                 "p0_ip_rep1","p0_ip_rep2","p5_ip_rep1","p5_ip_rep2",
                 "p10_ip_rep1","p10_ip_rep2","rp2_ip_rep1","rp2_ip_rep2")

caRNA_count_cds <- featureCounts(files=paste0(path,sample_name,".bam"),
                                   annot.ext = "/disk1/home/user_08/Data/annotation/mm39/gencode.vM28.annotation.gtf",
                                   isGTFAnnotationFile = TRUE, 
                                   GTF.featureType = "CDS", 
                                   GTF.attrType = "transcript_id",
                                   GTF.attrType.extra = "gene_id",
                                   isPairedEnd=TRUE,
                                   requireBothEndsMapped=TRUE,
                                   allowMultiOverlap = TRUE,
                                   fracOverlap = 0.5,
                                   countChimericFragments=FALSE,
                                   countMultiMappingReads = FALSE,
                                   strandSpecific = 2,checkFragLength = TRUE,
                                   minFragLength = 40, maxFragLength = 2000, 
                                   verbose = FALSE, nthreads=50)

caRNA_count_trans_df <- caRNA_count_cds$counts
colnames(caRNA_count_trans_df) <- sample_name

### caRNA expression level
caRNA_input_counts_p0p10 <- caRNA_count_trans_df[,c(1,2,5,6)]

library(DESeq2)
library(GEOquery)
library(GenomicFeatures)
library(BiocParallel)
library(GenomicAlignments)
library(rtracklayer)

# get annotation
anno <- rtracklayer::import("~/reference/annotation/mm19/gencode.vM28.annotation.gtf")
anno <- rtracklayer::as.data.frame(anno)

# get meta data
meta_data <- data.frame(rep = rep(c("rep1","rep2"),times=2),stage = rep(c("p0","p10"),each=2),row.names = colnames(caRNA_input_counts_p0p10),stringsAsFactors = TRUE)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = caRNA_input_counts_p0p10,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,name = "stage_p10_vs_p0")
summary(res)
res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange),]
rownames(res) <- substr(rownames(res),1,18)

## volcano
### PDUI p0p10
dapars_ca_p0p10 <- read.table("/disk/user_09/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p0p10.txt",
                              header = T,row.names = 1)

dapars_caRNA_p0p10_change <- dapars_ca_p0p10[which(dapars_ca_p0p10$filter != "NC"),]

dapars_caRNA_p0p10_change_hm <- as.matrix(dapars_caRNA_p0p10_change[,c(1,2)])
rownames(dapars_caRNA_p0p10_change_hm) <- substr(rownames(dapars_caRNA_p0p10_change_hm),1,18)


express_PDUI_change <- na.omit(res[rownames(dapars_caRNA_p0p10_change_hm),c("log2FoldChange","pvalue")])

as.data.frame(express_PDUI_change)
dapars_caRNA_p0p10_change_hm <- dapars_caRNA_p0p10_change_hm[rownames(express_PDUI_change),]
express_PDUI_change$PDUI <- log2((dapars_caRNA_p0p10_change_hm[,2]+0.0001)/(dapars_caRNA_p0p10_change_hm[,1]+0.0001))

express_PDUI_change$class <- ifelse(express_PDUI_change$PDUI>0,"lengthen","shorten")
express_PDUI_change <- as.data.frame(express_PDUI_change)

plot(x=express_PDUI_change$PDUI,y=express_PDUI_change$log2FoldChange)

library(ggpubr)
ggplot(data=express_PDUI_change)+
  geom_boxplot(mapping = aes(x=class,y=log2FoldChange,fill=class))


t.test(express_PDUI_change[which(express_PDUI_change$class=="lengthen"),"log2FoldChange"],express_PDUI_change[which(express_PDUI_change$class=="shorten"),"log2FoldChange"])


### mRNA_cds_expression

library(data.table)
mRNA_count_cds <- fread("/disk/user_09/user_08_TC1/10_mRNA/05_bam_hisat2/01_bam_sorted/11_bam_merge/05_gene_expression_tl/CDS_express_counts.tab")

mRNA_count_cds <- as.data.frame(mRNA_count_cds)
rownames(mRNA_count_cds) <- mRNA_count_cds$Geneid
mRNA_count_cds <- mRNA_count_cds[,-c(1:7)]
colnames(mRNA_count_cds) <- sample_name

mRNA_input_counts_p0p10 <- mRNA_count_cds[,c(1,2,5,6)]
# 
# # get annotation
# anno <- rtracklayer::import("~/reference/annotation/mm19/gencode.vM28.annotation.gtf")
# anno <- rtracklayer::as.data.frame(anno)
# 
# # get meta data
# meta_data <- data.frame(rep = rep(c("rep1","rep2"),times=2),stage = rep(c("p0","p10"),each=2),row.names = colnames(mRNA_input_counts_p0p10),stringsAsFactors = TRUE)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = mRNA_input_counts_p0p10,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,name = "stage_p10_vs_p0")
summary(res)
res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange),]
rownames(res) <- substr(rownames(res),1,18)


express_PDUI_change <- na.omit(res[rownames(dapars_caRNA_p0p10_change_hm),c("log2FoldChange","pvalue")])

as.data.frame(express_PDUI_change)
dapars_caRNA_p0p10_change_hm <- dapars_caRNA_p0p10_change_hm[rownames(express_PDUI_change),]
express_PDUI_change$PDUI <- log2((dapars_caRNA_p0p10_change_hm[,2]+0.0001)/(dapars_caRNA_p0p10_change_hm[,1]+0.0001))

express_PDUI_change$class <- ifelse(express_PDUI_change$PDUI>0,"lengthen","shorten")
express_PDUI_change <- as.data.frame(express_PDUI_change)

plot(x=express_PDUI_change$PDUI,y=express_PDUI_change$log2FoldChange)

library(ggpubr)
ggplot(data=express_PDUI_change)+
  geom_boxplot(mapping = aes(x=class,y=log2FoldChange,fill=class))


t.test(express_PDUI_change[which(express_PDUI_change$class=="lengthen"),"log2FoldChange"],express_PDUI_change[which(express_PDUI_change$class=="shorten"),"log2FoldChange"])




### mRNA_trans_expression

library(data.table)
mRNA_count_trans <- fread("/disk/user_09/user_08_TC1/10_mRNA/05_bam_hisat2/01_bam_sorted/11_bam_merge/05_gene_expression_tl/trans_express_counts.tab")

mRNA_count_trans <- as.data.frame(mRNA_count_trans)
rownames(mRNA_count_trans) <- mRNA_count_trans$Geneid
mRNA_count_trans <- mRNA_count_trans[,-c(1:7)]
colnames(mRNA_count_trans) <- sample_name

mRNA_input_counts_p0p10 <- mRNA_count_trans[,c(1,2,5,6)]
# 
# # get annotation
# anno <- rtracklayer::import("~/reference/annotation/mm19/gencode.vM28.annotation.gtf")
# anno <- rtracklayer::as.data.frame(anno)
# 
# # get meta data
# meta_data <- data.frame(rep = rep(c("rep1","rep2"),times=2),stage = rep(c("p0","p10"),each=2),row.names = colnames(mRNA_input_counts_p0p10),stringsAsFactors = TRUE)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = mRNA_input_counts_p0p10,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,name = "stage_p10_vs_p0")
summary(res)
res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange),]
rownames(res) <- substr(rownames(res),1,18)


express_PDUI_change <- na.omit(res[rownames(dapars_caRNA_p0p10_change_hm),c("log2FoldChange","pvalue")])

as.data.frame(express_PDUI_change)
dapars_caRNA_p0p10_change_hm <- dapars_caRNA_p0p10_change_hm[rownames(express_PDUI_change),]
express_PDUI_change$PDUI <- log2((dapars_caRNA_p0p10_change_hm[,2]+0.0001)/(dapars_caRNA_p0p10_change_hm[,1]+0.0001))

express_PDUI_change$class <- ifelse(express_PDUI_change$PDUI>0,"lengthen","shorten")
express_PDUI_change <- as.data.frame(express_PDUI_change)

plot(x=express_PDUI_change$PDUI,y=express_PDUI_change$log2FoldChange)

ggplot(data=express_PDUI_change)+
  geom_boxplot(mapping = aes(x=class,y=log2FoldChange,fill=class))


t.test(express_PDUI_change[which(express_PDUI_change$class=="lengthen"),"log2FoldChange"],express_PDUI_change[which(express_PDUI_change$class=="shorten"),"log2FoldChange"])

### mRNA_PDUI

dapars_m_p0p10 <- read.table("/disk/user_09/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_mRNA_p0p10.txt",
                              header = T,row.names = 1)

dapars_caRNA_p0p10_change <- dapars_ca_p0p10[which(dapars_ca_p0p10$filter != "NC"),]

dapars_caRNA_p0p10_change_hm <- as.matrix(dapars_caRNA_p0p10_change[,c(1,2)])
rownames(dapars_caRNA_p0p10_change_hm) <- substr(rownames(dapars_caRNA_p0p10_change_hm),1,18)


express_PDUI_change <- na.omit(res[rownames(dapars_caRNA_p0p10_change_hm),c("log2FoldChange","pvalue")])

as.data.frame(express_PDUI_change)
dapars_caRNA_p0p10_change_hm <- dapars_caRNA_p0p10_change_hm[rownames(express_PDUI_change),]
express_PDUI_change$PDUI <- log2((dapars_caRNA_p0p10_change_hm[,2]+0.0001)/(dapars_caRNA_p0p10_change_hm[,1]+0.0001))

express_PDUI_change$class <- ifelse(express_PDUI_change$PDUI>0,"lengthen","shorten")
express_PDUI_change <- as.data.frame(express_PDUI_change)


express_PDUI_change <- na.omit(res[rownames(dapars_caRNA_p0p10_change_hm),c("log2FoldChange","pvalue")])

as.data.frame(express_PDUI_change)
dapars_caRNA_p0p10_change_hm <- dapars_caRNA_p0p10_change_hm[rownames(express_PDUI_change),]
express_PDUI_change$PDUI <- log2((dapars_caRNA_p0p10_change_hm[,2]+0.0001)/(dapars_caRNA_p0p10_change_hm[,1]+0.0001))

express_PDUI_change$class <- ifelse(express_PDUI_change$PDUI>0,"lengthen","shorten")
express_PDUI_change <- as.data.frame(express_PDUI_change)

plot(x=express_PDUI_change$PDUI,y=express_PDUI_change$log2FoldChange)

ggplot(data=express_PDUI_change)+
  geom_boxplot(mapping = aes(x=class,y=log2FoldChange,fill=class))


t.test(express_PDUI_change[which(express_PDUI_change$class=="lengthen"),"log2FoldChange"],express_PDUI_change[which(express_PDUI_change$class=="shorten"),"log2FoldChange"])



