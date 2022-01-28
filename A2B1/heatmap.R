rm(list=ls())

library(DESeq2)
library(GEOquery)
library(GenomicFeatures)
library(BiocParallel)
library(GenomicAlignments)
library(rtracklayer)

#generate bam list
sample_name <- c("Lysate_input_rep1","Lysate_input_rep2","Result_input_rep1","Result_input_rep2","Lysate_IP_rep1","Lysate_IP_rep2","Result_IP_rep1","Result_IP_rep2")
setwd("~/LinLong/04_bam_raw/")
dir <- getwd()
fls <- file.path(dir,paste0(sample_name,'.bam'))
file.exists(fls)

bamlist <- BamFileList(fls,yieldSize = 2000000)

#generate features
txdb <- makeTxDbFromGFF("~/reference/annotation/mm19/gencode.vM28.annotation.gtf")
exonByGene <- exonsBy(txdb,by='gene')

#generate SE object
register(MulticoreParam(workers = 50))
SE <- summarizeOverlaps(features = exonByGene,reads = bamlist,mode = 'Union',singleEnd = F,ignore.strand = T)
#add coldata
colnames(SE) <- sample_name

# get raw counts
countdata <- assay(SE)

# get annotation
anno <- rtracklayer::import("~/reference/annotation/mm19/gencode.vM28.annotation.gtf")
anno <- rtracklayer::as.data.frame(anno)

# get meta data
express_matrix <- countdata[,1:4]
meta_data <- data.frame(rep = c("rep1","rep2","rep1","rep2"),treatment = c("Lysate","Lysate","Result","Result"),row.names = sample_name[1:4],stringsAsFactors = TRUE)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = express_matrix,colData = meta_data,design = ~ treatment)

dds <- dds[rowSums(counts(dds))>1,]

# pca plot
library(ggplot2)
rld <- rlog(dds,blind = FALSE)

pcaData <- plotPCA(rld,intgroup = "treatment",returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData,"percentVar"))
percentVar

ggplot(pcaData,aes(x=PC1,y=PC2,color=treatment))+
  geom_point(size=3)+
  xlab(paste0("PC1:",percentVar[1],"% variance"))+
  ylab(paste0("PC2:",percentVar[2],"% variance"))+
  coord_fixed()+
  ggtitle("PCA plot")

# DEG
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,name = "treatment_Result_vs_Lysate")
summary(res)
res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange),]

res_select <- res[abs(res$log2FoldChange) > 0.5 & res$padj < 0.05,]
summary(res_select)

# convert to TPM
# modify annotation
anno <- anno[anno$type == 'gene',]
table(duplicated(anno$gene_id))
rownames(anno) <- as.character(anno$gene_id)
table(rownames(express_matrix) %in% rownames(anno))

tpm_matrix <- do.call(cbind,base::lapply(colnames(express_matrix),FUN = function(x){
  counts <- express_matrix[,x]
  effLen <- anno$width
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  return(exp(rate - denom + log(1e6)))
}))
colnames(tpm_matrix) <- colnames(express_matrix)
table(colSums(tpm_matrix))

# heatmap
library(ComplexHeatmap)
tpm_matrix <- tpm_matrix[rownames(res_select),]
tpm_matrix <- t(scale(t(tpm_matrix)))

# create annatation
cell_anno <- HeatmapAnnotation(treat = meta_data[colnames(tpm_matrix),"treatment"])

res_select$express <- NA
res_select[res_select$log2FoldChange >= 0,"express"] <- "UP"
res_select[res_select$log2FoldChange < 0,"express"] <- "DOWN"

gene_anno <- HeatmapAnnotation(express = res_select[rownames(tpm_matrix),"express"],which = 'row')
Heatmap(matrix = tpm_matrix,show_column_names =FALSE,show_row_names=FALSE,
        cluster_rows = TRUE,cluster_columns = TRUE,
        top_annotation = cell_anno, left_annotation = gene_anno,
        width = unit(2,'inches'))

# m6A level
# scale factor
total_reads <- colSums(countdata)

library(edgeR)
myCPM <- cpm(countdata)

compare <- data.frame(Lysate_rep1=log2(myCPM[,5]/myCPM[,1]*0.93),Lysate_rep2=log2(myCPM[,6]/myCPM[,2]*1.07),Result_rep1=log2(myCPM[,7]/myCPM[,3]*0.77),Result_rep2=log2(myCPM[,8]/myCPM[,4]*0.77))

select_m6A <- as.matrix(compare[rownames(compare) %in% rownames(tpm_matrix),])




select_m6A <- t(scale(t(select_m6A)))
summary(is.na(select_m6A))
summary(is.nan(select_m6A))
summary(is.infinite(select_m6A))

select_m6A <- select_m6A[which(rowSums(!is.na(select_m6A))==4),]
tpm_matrix <- tpm_matrix[which(rowSums(!is.na(select_m6A))==4),]

Heatmap(matrix = select_m6A,show_column_names =TRUE,show_row_names=FALSE,
        cluster_rows = TRUE,cluster_columns = TRUE,
        width = unit(2,'inches'))


compare_scale <- t(scale(t(compare)))
summary(is.na(compare_scale))
summary(is.nan(compare_scale))
summary(is.infinite(compare_scale))
compare_scale <- compare_scale[which(rowSums(!is.na(compare_scale))==4),]

Heatmap(matrix = compare_scale,show_column_names =TRUE,show_row_names=FALSE,
        cluster_rows = TRUE,cluster_columns = TRUE,
        width = unit(2,'inches'))


pdf("./20211230/16_diff_peaks/A2B1_m6A_RNA.pdf",width = 5, height = 8)
mat_m6A <- read.csv()
fivenum(tpm_matrix)
heatmap_ca <- Heatmap(tpm_matrix,show_row_dend = FALSE,  
                      #row_km=3,
                      #column_km=2, 
                      #cluster_rows = FALSE,
                      show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                      heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
fivenum(mat_rko[select,])
heatmap_total <- Heatmap(select_m6A,show_row_dend = FALSE,  
                         #row_km=3,
                         #column_km=2, 
                         #cluster_rows = FALSE,
                         show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                         heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
ht_list =heatmap_ca+heatmap_total
draw(ht_list,heatmap_legend_side = "top")
dev.off()