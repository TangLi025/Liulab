### 把featureCounts换成cds,避免3UTR长度的影响

#featureCounts对cds计数
library(Rsubread)
library(ggplot2)
library(ggpubr)

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


### p0p10
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
res_ca <- results(dds,name = "stage_p10_vs_p0")
summary(res_ca)
res_ca <- res_ca[!is.na(res_ca$padj) & !is.na(res_ca$log2FoldChange),]
rownames(res_ca) <- substr(rownames(res_ca),1,18)

### mRNA_cds_expression

library(data.table)
mRNA_count_cds <- fread("/disk/user_09/Data/user_08_TC1/10_mRNA/05_bam_hisat2/01_bam_sorted/11_bam_merge/05_gene_expression_tl/CDS_express_counts.tab")


mRNA_count_cds <- as.data.frame(mRNA_count_cds)
rownames(mRNA_count_cds) <- mRNA_count_cds$Geneid
mRNA_count_cds <- mRNA_count_cds[,-c(1:7)]
colnames(mRNA_count_cds) <- sample_name

mRNA_counts_cds_p0p10 <- mRNA_count_cds[,c(1,2,5,6)]

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = mRNA_counts_cds_p0p10,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_m <- results(dds,name = "stage_p10_vs_p0")
summary(res_m)
res_m <- res_m[!is.na(res_m$padj) & !is.na(res_m$log2FoldChange),]
rownames(res_m) <- substr(rownames(res_m),1,18)

### 整理表达量与PDUI数据

#### PDUI mRNA表达量

express_PDUI <- as.data.frame(na.omit(res_m[rownames(dapars_caRNA_p0p10),c("log2FoldChange","pvalue")]))

dapars_caRNA_p0p10_hm <- dapars_caRNA_p0p10[rownames(express_PDUI),]
express_PDUI$PDUI <- dapars_caRNA_p0p10_hm[,2]-dapars_caRNA_p0p10_hm[,1]

express_PDUI$class <- ifelse(dapars_caRNA_p0p10_hm$filter=="NC","NC",ifelse(dapars_caRNA_p0p10_hm$filter=="DOWN","lengthen","shorten"))

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_mDEG_p0p10_v2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$log2FoldChange)

ggplot(data=express_PDUI,aes(x=class,y=log2FoldChange,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p0p10")+
  theme_bw()
dev.off()

t.test(express_PDUI[which(express_PDUI$class=="lengthen"),"log2FoldChange"],express_PDUI[which(express_PDUI$class=="shorten"),"log2FoldChange"])


### PDUI caRNA p0p10
res_ca_hm <- na.omit(as.data.frame(res_ca)[rownames(express_PDUI),])
express_PDUI <- express_PDUI[rownames(res_ca_hm),]
express_PDUI$caRNA_lfc2 <- res_ca_hm$log2FoldChange

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_caDEG_p0p10.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$caRNA_lfc2)

ggplot(data=express_PDUI,aes(x=class,y=caRNA_lfc2,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p0p10")+
  theme_bw()
dev.off()

t.test(express_PDUI[which(express_PDUI$class=="lengthen"),"log2FoldChange"],express_PDUI[which(express_PDUI$class=="shorten"),"log2FoldChange"])

### PDUI caRNA-mRNA p0p10
express_PDUI$ca_m_lfc2 <- express_PDUI$caRNA_lfc2-express_PDUI$log2FoldChange

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca_mDEG_p0p10_v2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$ca_m_lfc2)

ggplot(data=express_PDUI,aes(x=class,y=ca_m_lfc2,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p0p10")+
  theme_bw()
dev.off()


### mRNA_PDUI
dapars_m_p0p10 <- read.table("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_mRNA_p0p10.txt",
                             header = T,row.names = 1)

rownames(dapars_m_p0p10) <- substr(rownames(dapars_m_p0p10),1,18)
dapars_m_p0p10 <- dapars_m_p0p10[rownames(express_PDUI),]


express_PDUI$mPDUI <- dapars_m_p0p10[,2]-dapars_m_p0p10[,1]

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca_mPDUI_p0p10_v2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$mPDUI)

ggplot(data=express_PDUI,aes(x=class,y=mPDUI,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p0p10")+
  theme_bw()
dev.off()


#################### p0p5 ####################

### p0p5
### caRNA expression level

caRNA_input_counts_p0p5 <- caRNA_count_trans_df[,c(1,2,3,4)]

# get meta data
meta_data <- data.frame(rep = rep(c("rep1","rep2"),times=2),stage = rep(c("p0","p5"),each=2),
                        row.names = colnames(caRNA_input_counts_p0p5),stringsAsFactors = TRUE)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = caRNA_input_counts_p0p5,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_ca <- results(dds,name = "stage_p5_vs_p0")
summary(res_ca)
res_ca <- res_ca[!is.na(res_ca$padj) & !is.na(res_ca$log2FoldChange),]
rownames(res_ca) <- substr(rownames(res_ca),1,18)

### mRNA_cds_expression

mRNA_counts_cds_p0p5 <- mRNA_count_cds[,c(1,2,3,4)]

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = mRNA_counts_cds_p0p5,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_m <- results(dds,name = "stage_p5_vs_p0")
summary(res_m)
res_m <- res_m[!is.na(res_m$padj) & !is.na(res_m$log2FoldChange),]
rownames(res_m) <- substr(rownames(res_m),1,18)

### 整理表达量与PDUI数据

#### PDUI mRNA表达量
dapars_caRNA_p0p5 <- as.data.frame(fread("~/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p0p5.txt"))
rownames(dapars_caRNA_p0p5) <- substr(dapars_caRNA_p0p5$Gene,1,18)
dapars_caRNA_p0p5 <- dapars_caRNA_p0p5[,-1]

express_PDUI <- as.data.frame(na.omit(res_m[rownames(dapars_caRNA_p0p5),c("log2FoldChange","pvalue")]))

dapars_caRNA_p0p5_hm <- dapars_caRNA_p0p5[rownames(express_PDUI),]
express_PDUI$PDUI <- dapars_caRNA_p0p5_hm[,2]-dapars_caRNA_p0p5_hm[,1]

express_PDUI$class <- ifelse(dapars_caRNA_p0p5_hm$filter=="NC","NC",ifelse(dapars_caRNA_p0p5_hm$filter=="DOWN","lengthen","shorten"))

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_mDEG_p0p5.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$log2FoldChange)

ggplot(data=express_PDUI,aes(x=class,y=log2FoldChange,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p0p5")+
  theme_bw()
dev.off()

t.test(express_PDUI[which(express_PDUI$class=="lengthen"),"log2FoldChange"],express_PDUI[which(express_PDUI$class=="shorten"),"log2FoldChange"])


### PDUI caRNA p0p5
res_ca_hm <- na.omit(as.data.frame(res_ca)[rownames(express_PDUI),])
express_PDUI <- express_PDUI[rownames(res_ca_hm),]
express_PDUI$caRNA_lfc2 <- res_ca_hm$log2FoldChange

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_caDEG_p0p5.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$caRNA_lfc2)

ggplot(data=express_PDUI,aes(x=class,y=caRNA_lfc2,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p0p5")+
  theme_bw()
dev.off()

t.test(express_PDUI[which(express_PDUI$class=="lengthen"),"log2FoldChange"],express_PDUI[which(express_PDUI$class=="shorten"),"log2FoldChange"])

### PDUI caRNA-mRNA p0p5
express_PDUI$ca_m_lfc2 <- express_PDUI$caRNA_lfc2-express_PDUI$log2FoldChange

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca_mDEG_p0p5.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$ca_m_lfc2)

ggplot(data=express_PDUI,aes(x=class,y=ca_m_lfc2,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p0p5")+
  theme_bw()
dev.off()


### mRNA_PDUI
dapars_m_p0p5 <- read.table("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_mRNA_p0p5.txt",
                             header = T,row.names = 1)

rownames(dapars_m_p0p5) <- substr(rownames(dapars_m_p0p5),1,18)
dapars_m_p0p5 <- dapars_m_p0p5[rownames(express_PDUI),]


express_PDUI$mPDUI <- dapars_m_p0p5[,2]-dapars_m_p0p5[,1]

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca_mPDUI_p0p5_v2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$mPDUI)

ggplot(data=express_PDUI,aes(x=class,y=mPDUI,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p0p5")+
  theme_bw()
dev.off()

#################### p0rp2 ####################

### p0rp2
### caRNA expression level

caRNA_input_counts_p0rp2 <- caRNA_count_trans_df[,c(1,2,7,8)]

# get meta data
meta_data <- data.frame(rep = rep(c("rep1","rep2"),times=2),stage = rep(c("p0","rp2"),each=2),
                        row.names = colnames(caRNA_input_counts_p0rp2),stringsAsFactors = TRUE)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = caRNA_input_counts_p0rp2,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_ca <- results(dds,name = "stage_rp2_vs_p0")
summary(res_ca)
res_ca <- res_ca[!is.na(res_ca$padj) & !is.na(res_ca$log2FoldChange),]
rownames(res_ca) <- substr(rownames(res_ca),1,18)

### mRNA_cds_expression

mRNA_counts_cds_p0rp2 <- mRNA_count_cds[,c(1,2,7,8)]

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = mRNA_counts_cds_p0rp2,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_m <- results(dds,name = "stage_rp2_vs_p0")
summary(res_m)
res_m <- res_m[!is.na(res_m$padj) & !is.na(res_m$log2FoldChange),]
rownames(res_m) <- substr(rownames(res_m),1,18)

### 整理表达量与PDUI数据

#### PDUI mRNA表达量
dapars_caRNA_p0rp2 <- as.data.frame(fread("~/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p0rp2.txt"))
rownames(dapars_caRNA_p0rp2) <- substr(dapars_caRNA_p0rp2$Gene,1,18)
dapars_caRNA_p0rp2 <- dapars_caRNA_p0rp2[,-1]

express_PDUI <- as.data.frame(na.omit(res_m[rownames(dapars_caRNA_p0rp2),c("log2FoldChange","pvalue")]))

dapars_caRNA_p0rp2_hm <- dapars_caRNA_p0rp2[rownames(express_PDUI),]
express_PDUI$PDUI <- dapars_caRNA_p0rp2_hm[,2]-dapars_caRNA_p0rp2_hm[,1]

express_PDUI$class <- ifelse(dapars_caRNA_p0rp2_hm$filter=="NC","NC",ifelse(dapars_caRNA_p0rp2_hm$filter=="DOWN","lengthen","shorten"))

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_mDEG_p0rp2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$log2FoldChange)

ggplot(data=express_PDUI,aes(x=class,y=log2FoldChange,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p0rp2")+
  theme_bw()
dev.off()

t.test(express_PDUI[which(express_PDUI$class=="lengthen"),"log2FoldChange"],express_PDUI[which(express_PDUI$class=="shorten"),"log2FoldChange"])


### PDUI caRNA p0rp2
res_ca_hm <- na.omit(as.data.frame(res_ca)[rownames(express_PDUI),])
express_PDUI <- express_PDUI[rownames(res_ca_hm),]
express_PDUI$caRNA_lfc2 <- res_ca_hm$log2FoldChange

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_caDEG_p0rp2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$caRNA_lfc2)

ggplot(data=express_PDUI,aes(x=class,y=caRNA_lfc2,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p0rp2")+
  theme_bw()
dev.off()

t.test(express_PDUI[which(express_PDUI$class=="lengthen"),"log2FoldChange"],express_PDUI[which(express_PDUI$class=="shorten"),"log2FoldChange"])

### PDUI caRNA-mRNA p0rp2
express_PDUI$ca_m_lfc2 <- express_PDUI$caRNA_lfc2-express_PDUI$log2FoldChange

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca_mDEG_p0rp2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$ca_m_lfc2)

ggplot(data=express_PDUI,aes(x=class,y=ca_m_lfc2,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p0rp2")+
  theme_bw()
dev.off()


### mRNA_PDUI
dapars_m_p0rp2 <- read.table("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_mRNA_p0rp2.txt",
                            header = T,row.names = 1)

rownames(dapars_m_p0rp2) <- substr(rownames(dapars_m_p0rp2),1,18)
dapars_m_p0rp2 <- dapars_m_p0rp2[rownames(express_PDUI),]


express_PDUI$mPDUI <- dapars_m_p0rp2[,2]-dapars_m_p0rp2[,1]

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca_mPDUI_p0rp2_v2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$mPDUI)

ggplot(data=express_PDUI,aes(x=class,y=mPDUI,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p0rp2")+
  theme_bw()
dev.off()

#################### p5rp2 ####################

### p5rp2
### caRNA expression level

caRNA_input_counts_p5rp2 <- caRNA_count_trans_df[,c(3,4,7,8)]

# get meta data
meta_data <- data.frame(rep = rep(c("rep1","rep2"),times=2),stage = rep(c("p5","rp2"),each=2),
                        row.names = colnames(caRNA_input_counts_p5rp2),stringsAsFactors = TRUE)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = caRNA_input_counts_p5rp2,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_ca <- results(dds,name = "stage_rp2_vs_p5")
summary(res_ca)
res_ca <- res_ca[!is.na(res_ca$padj) & !is.na(res_ca$log2FoldChange),]
rownames(res_ca) <- substr(rownames(res_ca),1,18)

### mRNA_cds_expression

mRNA_counts_cds_p5rp2 <- mRNA_count_cds[,c(3,4,7,8)]

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = mRNA_counts_cds_p5rp2,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_m <- results(dds,name = "stage_rp2_vs_p5")
summary(res_m)
res_m <- res_m[!is.na(res_m$padj) & !is.na(res_m$log2FoldChange),]
rownames(res_m) <- substr(rownames(res_m),1,18)

### 整理表达量与PDUI数据

#### PDUI mRNA表达量
dapars_caRNA_p5rp2 <- as.data.frame(fread("~/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p5rp2.txt"))
rownames(dapars_caRNA_p5rp2) <- substr(dapars_caRNA_p5rp2$Gene,1,18)
dapars_caRNA_p5rp2 <- dapars_caRNA_p5rp2[,-1]

express_PDUI <- as.data.frame(na.omit(res_m[rownames(dapars_caRNA_p5rp2),c("log2FoldChange","pvalue")]))

dapars_caRNA_p5rp2_hm <- dapars_caRNA_p5rp2[rownames(express_PDUI),]
express_PDUI$PDUI <- dapars_caRNA_p5rp2_hm[,2]-dapars_caRNA_p5rp2_hm[,1]

express_PDUI$class <- ifelse(dapars_caRNA_p5rp2_hm$filter=="NC","NC",ifelse(dapars_caRNA_p5rp2_hm$filter=="DOWN","lengthen","shorten"))

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_mDEG_p5rp2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$log2FoldChange)

ggplot(data=express_PDUI,aes(x=class,y=log2FoldChange,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p5rp2")+
  theme_bw()
dev.off()

t.test(express_PDUI[which(express_PDUI$class=="lengthen"),"log2FoldChange"],express_PDUI[which(express_PDUI$class=="shorten"),"log2FoldChange"])


### PDUI caRNA p5rp2
res_ca_hm <- na.omit(as.data.frame(res_ca)[rownames(express_PDUI),])
express_PDUI <- express_PDUI[rownames(res_ca_hm),]
express_PDUI$caRNA_lfc2 <- res_ca_hm$log2FoldChange

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_caDEG_p5rp2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$caRNA_lfc2)

ggplot(data=express_PDUI,aes(x=class,y=caRNA_lfc2,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p5rp2")+
  theme_bw()
dev.off()

t.test(express_PDUI[which(express_PDUI$class=="lengthen"),"log2FoldChange"],express_PDUI[which(express_PDUI$class=="shorten"),"log2FoldChange"])

### PDUI caRNA-mRNA p5rp2
express_PDUI$ca_m_lfc2 <- express_PDUI$caRNA_lfc2-express_PDUI$log2FoldChange

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca_mDEG_p5rp2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$ca_m_lfc2)

ggplot(data=express_PDUI,aes(x=class,y=ca_m_lfc2,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p5rp2")+
  theme_bw()
dev.off()


### mRNA_PDUI
dapars_m_p5rp2 <- read.table("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_mRNA_p5rp2.txt",
                             header = T,row.names = 1)

rownames(dapars_m_p5rp2) <- substr(rownames(dapars_m_p5rp2),1,18)
dapars_m_p5rp2 <- dapars_m_p5rp2[rownames(express_PDUI),]


express_PDUI$mPDUI <- dapars_m_p5rp2[,2]-dapars_m_p5rp2[,1]

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca_mPDUI_p5rp2_v2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$mPDUI)

ggplot(data=express_PDUI,aes(x=class,y=mPDUI,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p5rp2")+
  theme_bw()
dev.off()

#################### p10rp2 ####################

### p10rp2
### caRNA expression level

caRNA_input_counts_p10rp2 <- caRNA_count_trans_df[,c(5,6,7,8)]

# get meta data
meta_data <- data.frame(rep = rep(c("rep1","rep2"),times=2),stage = rep(c("p10","rp2"),each=2),
                        row.names = colnames(caRNA_input_counts_p10rp2),stringsAsFactors = TRUE)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = caRNA_input_counts_p10rp2,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_ca <- results(dds,name = "stage_rp2_vs_p10")
summary(res_ca)
res_ca <- res_ca[!is.na(res_ca$padj) & !is.na(res_ca$log2FoldChange),]
rownames(res_ca) <- substr(rownames(res_ca),1,18)

### mRNA_cds_expression

mRNA_counts_cds_p10rp2 <- mRNA_count_cds[,c(5,6,7,8)]

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = mRNA_counts_cds_p10rp2,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_m <- results(dds,name = "stage_rp2_vs_p10")
summary(res_m)
res_m <- res_m[!is.na(res_m$padj) & !is.na(res_m$log2FoldChange),]
rownames(res_m) <- substr(rownames(res_m),1,18)

### 整理表达量与PDUI数据

#### PDUI mRNA表达量
dapars_caRNA_p10rp2 <- as.data.frame(fread("~/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p10rp2.txt"))
rownames(dapars_caRNA_p10rp2) <- substr(dapars_caRNA_p10rp2$Gene,1,18)
dapars_caRNA_p10rp2 <- dapars_caRNA_p10rp2[,-1]

express_PDUI <- as.data.frame(na.omit(res_m[rownames(dapars_caRNA_p10rp2),c("log2FoldChange","pvalue")]))

dapars_caRNA_p10rp2_hm <- dapars_caRNA_p10rp2[rownames(express_PDUI),]
express_PDUI$PDUI <- dapars_caRNA_p10rp2_hm[,2]-dapars_caRNA_p10rp2_hm[,1]

express_PDUI$class <- ifelse(dapars_caRNA_p10rp2_hm$filter=="NC","NC",ifelse(dapars_caRNA_p10rp2_hm$filter=="DOWN","lengthen","shorten"))

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_mDEG_p10rp2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$log2FoldChange)

ggplot(data=express_PDUI,aes(x=class,y=log2FoldChange,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p10rp2")+
  theme_bw()
dev.off()

t.test(express_PDUI[which(express_PDUI$class=="lengthen"),"log2FoldChange"],express_PDUI[which(express_PDUI$class=="shorten"),"log2FoldChange"])


### PDUI caRNA p10rp2
res_ca_hm <- na.omit(as.data.frame(res_ca)[rownames(express_PDUI),])
express_PDUI <- express_PDUI[rownames(res_ca_hm),]
express_PDUI$caRNA_lfc2 <- res_ca_hm$log2FoldChange

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_caDEG_p10rp2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$caRNA_lfc2)

ggplot(data=express_PDUI,aes(x=class,y=caRNA_lfc2,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p10rp2")+
  theme_bw()
dev.off()

t.test(express_PDUI[which(express_PDUI$class=="lengthen"),"log2FoldChange"],express_PDUI[which(express_PDUI$class=="shorten"),"log2FoldChange"])

### PDUI caRNA-mRNA p10rp2
express_PDUI$ca_m_lfc2 <- express_PDUI$caRNA_lfc2-express_PDUI$log2FoldChange

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca_mDEG_p10rp2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$ca_m_lfc2)

ggplot(data=express_PDUI,aes(x=class,y=ca_m_lfc2,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p10rp2")+
  theme_bw()
dev.off()


### mRNA_PDUI
dapars_m_p10rp2 <- read.table("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_mRNA_p10rp2.txt",
                             header = T,row.names = 1)

rownames(dapars_m_p10rp2) <- substr(rownames(dapars_m_p10rp2),1,18)
dapars_m_p10rp2 <- dapars_m_p10rp2[rownames(express_PDUI),]


express_PDUI$mPDUI <- dapars_m_p10rp2[,2]-dapars_m_p10rp2[,1]

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca_mPDUI_p10rp2_v2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$mPDUI)

ggplot(data=express_PDUI,aes(x=class,y=mPDUI,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p10rp2")+
  theme_bw()
dev.off()

#################### p5p10 ####################

### p5p10
### caRNA expression level

caRNA_input_counts_p5p10 <- caRNA_count_trans_df[,c(3,4,5,6)]

# get meta data
meta_data <- data.frame(rep = rep(c("rep1","rep2"),times=2),stage = rep(c("p5","p10"),each=2),
                        row.names = colnames(caRNA_input_counts_p5p10),stringsAsFactors = TRUE)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = caRNA_input_counts_p5p10,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_ca <- results(dds,name = "stage_p10_vs_p5")
summary(res_ca)
res_ca <- res_ca[!is.na(res_ca$padj) & !is.na(res_ca$log2FoldChange),]
rownames(res_ca) <- substr(rownames(res_ca),1,18)

### mRNA_cds_expression

mRNA_counts_cds_p5p10 <- mRNA_count_cds[,c(3,4,5,6)]

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = mRNA_counts_cds_p5p10,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_m <- results(dds,name = "stage_p10_vs_p5")
summary(res_m)
res_m <- res_m[!is.na(res_m$padj) & !is.na(res_m$log2FoldChange),]
rownames(res_m) <- substr(rownames(res_m),1,18)

### 整理表达量与PDUI数据

#### PDUI mRNA表达量
dapars_caRNA_p5p10 <- as.data.frame(fread("~/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p5p10.txt"))
rownames(dapars_caRNA_p5p10) <- substr(dapars_caRNA_p5p10$Gene,1,18)
dapars_caRNA_p5p10 <- dapars_caRNA_p5p10[,-1]

express_PDUI <- as.data.frame(na.omit(res_m[rownames(dapars_caRNA_p5p10),c("log2FoldChange","pvalue")]))

dapars_caRNA_p5p10_hm <- dapars_caRNA_p5p10[rownames(express_PDUI),]
express_PDUI$PDUI <- dapars_caRNA_p5p10_hm[,2]-dapars_caRNA_p5p10_hm[,1]

express_PDUI$class <- ifelse(dapars_caRNA_p5p10_hm$filter=="NC","NC",ifelse(dapars_caRNA_p5p10_hm$filter=="DOWN","lengthen","shorten"))

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_mDEG_p5p10.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$log2FoldChange)

ggplot(data=express_PDUI,aes(x=class,y=log2FoldChange,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p5p10")+
  theme_bw()
dev.off()

t.test(express_PDUI[which(express_PDUI$class=="lengthen"),"log2FoldChange"],express_PDUI[which(express_PDUI$class=="shorten"),"log2FoldChange"])


### PDUI caRNA p5p10
res_ca_hm <- na.omit(as.data.frame(res_ca)[rownames(express_PDUI),])
express_PDUI <- express_PDUI[rownames(res_ca_hm),]
express_PDUI$caRNA_lfc2 <- res_ca_hm$log2FoldChange

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_caDEG_p5p10.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$caRNA_lfc2)

ggplot(data=express_PDUI,aes(x=class,y=caRNA_lfc2,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p5p10")+
  theme_bw()
dev.off()

t.test(express_PDUI[which(express_PDUI$class=="lengthen"),"log2FoldChange"],express_PDUI[which(express_PDUI$class=="shorten"),"log2FoldChange"])

### PDUI caRNA-mRNA p5p10
express_PDUI$ca_m_lfc2 <- express_PDUI$caRNA_lfc2-express_PDUI$log2FoldChange

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca_mDEG_p5p10.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$ca_m_lfc2)

ggplot(data=express_PDUI,aes(x=class,y=ca_m_lfc2,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p5p10")+
  theme_bw()
dev.off()


### mRNA_PDUI
dapars_m_p5p10 <- read.table("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_mRNA_p5p10.txt",
                              header = T,row.names = 1)

rownames(dapars_m_p5p10) <- substr(rownames(dapars_m_p5p10),1,18)
dapars_m_p5p10 <- dapars_m_p5p10[rownames(express_PDUI),]


express_PDUI$mPDUI <- dapars_m_p5p10[,2]-dapars_m_p5p10[,1]

pdf("~/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca_mPDUI_p5p10_v2.pdf",width = 3.5,height = 3.5)
plot(x=express_PDUI$PDUI,y=express_PDUI$mPDUI)

ggplot(data=express_PDUI,aes(x=class,y=mPDUI,fill=class))+
  geom_boxplot()+
  stat_compare_means(comparisons=list(c("NC","shorten"),c("lengthen","shorten"),c("NC","lengthen")),
                     label="p.signif",method="t.test")+
  ggtitle("p5p10")+
  theme_bw()
dev.off()