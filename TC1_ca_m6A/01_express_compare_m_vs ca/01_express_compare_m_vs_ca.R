## project: 01_express_compare_m_vs_ca
## Data: 20220413

library(Rsubread)
library(ggpubr)
library(ggplot2)
library(DESeq2)

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
                                 GTF.attrType.extra = c("gene_id","gene_name"),
                                 isPairedEnd=TRUE,
                                 requireBothEndsMapped=TRUE,
                                 allowMultiOverlap = TRUE,
                                 fracOverlap = 0.5,
                                 countChimericFragments=FALSE,
                                 countMultiMappingReads = FALSE,
                                 strandSpecific = 2,checkFragLength = TRUE,
                                 minFragLength = 40, maxFragLength = 2000, 
                                 verbose = FALSE, nthreads=50)

caRNA_count_cds_df <- caRNA_count_cds$counts
colnames(caRNA_count_cds_df) <- sample_name
caRNA_input_cds_counts_p0p10 <- caRNA_count_cds_df[,c(1,2,5,6)]

# get annotation
anno <- rtracklayer::import("~/reference/annotation/mm19/gencode.vM28.annotation.gtf")
anno <- rtracklayer::as.data.frame(anno)

# get meta data
meta_data <- data.frame(rep = rep(c("rep1","rep2"),times=2),stage = rep(c("p0","p10"),each=2),row.names = colnames(caRNA_input_cds_counts_p0p10),stringsAsFactors = TRUE)

# create DESeqDataSet

dds <- DESeqDataSetFromMatrix(countData = caRNA_input_cds_counts_p0p10,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_ca <- results(dds,name = "stage_p10_vs_p0")
summary(res_ca)
res_ca <- res_ca[!is.na(res_ca$padj) & !is.na(res_ca$log2FoldChange),]
rownames(res_ca) <- substr(rownames(res_ca),1,18)

### PDUI caRNA p0p10
dapars_ca_p0p10 <- read.table("/disk/user_09/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p0p10.txt",
                              header = T,row.names = 1)

dapars_ca_p0p10_hm <- dapars_ca_p0p10
rownames(dapars_ca_p0p10_hm) <- substr(rownames(dapars_ca_p0p10_hm),1,18)


express_PDUI_change <- na.omit(res_ca[rownames(dapars_ca_p0p10_hm),c("log2FoldChange","pvalue")])

#as.data.frame(express_PDUI_change)
dapars_ca_p0p10_hm <- dapars_ca_p0p10_hm[rownames(express_PDUI_change),]
express_PDUI_change$PDUI <- dapars_ca_p0p10_hm[,2]-dapars_ca_p0p10_hm[,1]

express_PDUI_change$class <- dapars_ca_p0p10_hm$filter
express_PDUI_change <- as.data.frame(express_PDUI_change)

dir.create("/disk/user_09/Data/03_TC1_caRNA/05_PDUI_expression")
pdf(file = "/disk/user_09/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_caDEG_p0p10.pdf")
plot(x=express_PDUI_change$PDUI,y=express_PDUI_change$log2FoldChange)


ggplot(data=express_PDUI_change)+
  geom_boxplot(mapping = aes(x=class,y=log2FoldChange,fill=class))
dev.off()


t.test(express_PDUI_change[which(express_PDUI_change$class=="NC"),"log2FoldChange"],express_PDUI_change[which(express_PDUI_change$class=="DOWN"),"log2FoldChange"])

t.test(express_PDUI_change[which(express_PDUI_change$class=="NC"),"log2FoldChange"],express_PDUI_change[which(express_PDUI_change$class=="UP"),"log2FoldChange"])


### mRNA p0p10
library(data.table)
mRNA_count_cds <- fread("/disk/user_09/user_08_TC1/10_mRNA/05_bam_hisat2/01_bam_sorted/11_bam_merge/05_gene_expression_tl/CDS_express_counts.tab")

mRNA_count_cds <- as.data.frame(mRNA_count_cds)
rownames(mRNA_count_cds) <- substr(mRNA_count_cds$Geneid,1,18)
mRNA_count_cds <- mRNA_count_cds[,-c(1:7)]
colnames(mRNA_count_cds) <- sample_name


mRNA_input_cds_counts_p0p10 <- mRNA_count_cds[,c(1,2,5,6)]

# get meta data
meta_data <- data.frame(rep = rep(c("rep1","rep2"),times=2),stage = rep(c("p0","p10"),each=2),row.names = colnames(mRNA_input_cds_counts_p0p10),stringsAsFactors = TRUE)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = mRNA_input_cds_counts_p0p10,colData = meta_data,design = ~ stage)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res_m <- results(dds,name = "stage_p10_vs_p0")
summary(res_m)
res_m <- res_m[!is.na(res_m$padj) & !is.na(res_m$log2FoldChange),]
rownames(res_m) <- substr(rownames(res_m),1,18)

### PDUI mRNA p0p10
dapars_ca_p0p10 <- read.table("/disk/user_09/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p0p10.txt",
                              header = T,row.names = 1)
rownames(dapars_ca_p0p10) <- substr(rownames(dapars_ca_p0p10),1,18)
dapars_ca_p0p10_hm <- dapars_ca_p0p10

express_PDUI_change <- na.omit(res_m[rownames(dapars_ca_p0p10_hm),c("log2FoldChange","pvalue")])

dapars_ca_p0p10_hm <- dapars_ca_p0p10_hm[rownames(express_PDUI_change),]
express_PDUI_change$PDUI <- dapars_ca_p0p10_hm[,2]-dapars_ca_p0p10_hm[,1]

express_PDUI_change$class <- dapars_ca_p0p10_hm$filter
express_PDUI_change <- as.data.frame(express_PDUI_change)

pdf(file = "/disk/user_09/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_mDEG_p0p10.pdf")
plot(x=express_PDUI_change$PDUI,y=express_PDUI_change$log2FoldChange)

ggplot(data=express_PDUI_change)+
  geom_boxplot(mapping = aes(x=class,y=log2FoldChange,fill=class))
dev.off()

t.test(express_PDUI_change[which(express_PDUI_change$class=="UP"),"log2FoldChange"],express_PDUI_change[which(express_PDUI_change$class=="DOWN"),"log2FoldChange"])

t.test(express_PDUI_change[which(express_PDUI_change$class=="NC"),"log2FoldChange"],express_PDUI_change[which(express_PDUI_change$class=="DOWN"),"log2FoldChange"])

t.test(express_PDUI_change[which(express_PDUI_change$class=="NC"),"log2FoldChange"],express_PDUI_change[which(express_PDUI_change$class=="UP"),"log2FoldChange"])

express_PDUI_change$ca_lfc <- res_ca[rownames(express_PDUI_change),"log2FoldChange"]
express_PDUI_change <- na.omit(express_PDUI_change)

express_PDUI_change$ca_m_lfc2 <- express_PDUI_change$ca_lfc - express_PDUI_change$log2FoldChange 

pdf(file = "/disk/user_09/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_ca-mDEG_p0p10.pdf")
plot(x=express_PDUI_change$PDUI,y=express_PDUI_change$ca_m_lfc2)

ggplot(data=express_PDUI_change)+
  geom_boxplot(mapping = aes(x=class,y=ca_m_lfc2,fill=class))
dev.off()

t.test(express_PDUI_change[which(express_PDUI_change$class=="UP"),"ca_m_lfc2"],express_PDUI_change[which(express_PDUI_change$class=="DOWN"),"ca_m_lfc2"])

t.test(express_PDUI_change[which(express_PDUI_change$class=="NC"),"ca_m_lfc2"],express_PDUI_change[which(express_PDUI_change$class=="DOWN"),"ca_m_lfc2"])

t.test(express_PDUI_change[which(express_PDUI_change$class=="NC"),"ca_m_lfc2"],express_PDUI_change[which(express_PDUI_change$class=="UP"),"ca_m_lfc2"])


### mito
p0p10_short_mito <- read.delim("/disk/user_09/Data/03_TC1_caRNA/05_PDUI_expression/mito_short_GO.txt",header=T)

p0p10_all_mito <- read.delim("/disk/user_09/Data/03_TC1_caRNA/05_PDUI_expression/mito_all_GO.txt",header = T)

PDUI_p0p10_mito <- dapars_ca_p0p10[rownames(p0p10_all_mito),]

express_mito_ca <- na.omit(res_ca[which(rownames(res_ca) %in% rownames(PDUI_p0p10_mito)),c("log2FoldChange","pvalue")])

PDUI_p0p10_mito <- PDUI_p0p10_mito[rownames(express_mito_ca),]
express_mito_ca$PDUI <- PDUI_p0p10_mito[,2]-PDUI_p0p10_mito[,1]

express_mito_ca$class <- PDUI_p0p10_mito$filter
express_mito_ca <- as.data.frame(express_mito_ca)


pdf(file = "/disk/user_09/Data/03_TC1_caRNA/05_PDUI_expression/caPDUI_caDEG_p0p10.pdf")
plot(x=express_mito_ca$PDUI,y=express_mito_ca$log2FoldChange)


ggplot(data=express_mito_ca)+
  geom_boxplot(mapping = aes(x=class,y=log2FoldChange,fill=class))
dev.off()


t.test(express_mito_ca[which(express_mito_ca$class=="NC"),"log2FoldChange"],express_mito_ca[which(express_mito_ca$class=="DOWN"),"log2FoldChange"])

t.test(express_mito_ca[which(express_mito_ca$class=="NC"),"log2FoldChange"],express_mito_ca[which(express_mito_ca$class=="UP"),"log2FoldChange"])

