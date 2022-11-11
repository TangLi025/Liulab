rm(list=ls())

library(DESeq2)
library(GEOquery)
library(GenomicFeatures)
library(BiocParallel)
library(GenomicAlignments)
library(rtracklayer)

#generate bam list
sample_name <- c("Lysate.rep1.input","Lysate.rep2.input","Result.rep1.input","Result.rep2.input","Lysate.rep1.IP","Lysate.rep2.IP","Result.rep1.IP","Result.rep2.IP")
setwd("~/LinLong/04_bam_raw/")
dir <- getwd()
fls <- file.path(dir,paste0(sample_name,'.bam'))
file.exists(fls)

bamlist <- BamFileList(fls,yieldSize = 2000000)

library(ChIPpeakAnno)
#generate features
m6A_peak <- import("~/LinLong/10_bed_merge/raw/Lysate_Result_merge_peaks.bed")

#generate SE object
register(MulticoreParam(workers = 50))
SE <- summarizeOverlaps(features = m6A_peak,reads = bamlist,mode = 'Union',singleEnd = F,ignore.strand = T)
#add coldata
colnames(SE) <- sample_name

# get raw counts
countdata <- assay(SE)

write.table(countdata,"peak_counts.tab")
countdata <- read.table("peak_counts.tab")

peak_annotation <- read.delim("~/LinLong/10_bed_merge/raw/annotation_output/Lysate_Result_merge_peaks_annotation_homer_rearranged.xls")

# get annotation
anno <- rtracklayer::import("~/reference/annotation/mm19/gencode.vM28.annotation.gtf")
anno <- rtracklayer::as.data.frame(anno)

# get meta data
express_matrix <- countdata

nor_ratio <- c(1.950749499,1.063148729, 2.008697226,1.044198103,0.741532184,0.647182502,0.619590696,0.773117373)

express_matrix_nor <- data.frame(express_matrix[,5]*nor_ratio[5],express_matrix[,1]*nor_ratio[1],express_matrix[,6]*nor_ratio[6],express_matrix[,2]*nor_ratio[2],express_matrix[,7]*nor_ratio[7],express_matrix[,3]*nor_ratio[3],express_matrix[,8]*nor_ratio[8],express_matrix[,4]*nor_ratio[4])

express_matrix_log2FC <- data.frame(Lysate_rep1=log2(express_matrix_nor[,1]*0.93/express_matrix_nor[,2]),Lysate_rep2=log2(express_matrix_nor[,3]*1.07/express_matrix_nor[,4]),Result_rep1=log2(express_matrix_nor[,5]*0.78/express_matrix_nor[,6]),Result_rep2=log2(express_matrix_nor[,7]*0.78/express_matrix_nor[,8]))

express_matrix_log2FC <- as.matrix(express_matrix_log2FC)

peak_anno <- data.frame(peak_anno_merge=peak_annotation$Entrez.ID[which(rowSums(is.finite(express_matrix_log2FC))==4)])

express_matrix_nor <- express_matrix_nor[which(rowSums(is.finite(express_matrix_log2FC))==4),]

express_matrix_log2FC <- express_matrix_log2FC[which(rowSums(is.finite(express_matrix_log2FC))==4),]


peak_anno_pos <- peak_anno[which((express_matrix_log2FC[,3]-express_matrix_log2FC[,1]>0.5) & (express_matrix_log2FC[,4]-express_matrix_log2FC[,2]>0.5)),]
diff_peak_pos <- express_matrix_nor[which((express_matrix_log2FC[,3]-express_matrix_log2FC[,1]>0.5) & (express_matrix_log2FC[,4]-express_matrix_log2FC[,2]>0.5)),]

peak_anno_neg <- peak_anno[(express_matrix_log2FC[,3]-express_matrix_log2FC[,1]< -0.5) & (express_matrix_log2FC[,4]-express_matrix_log2FC[,2]< -0.5),]
diff_peak_neg <- express_matrix_nor[(express_matrix_log2FC[,3]-express_matrix_log2FC[,1]< -0.5) & (express_matrix_log2FC[,4]-express_matrix_log2FC[,2]< -0.5),]


diff_peak <- data.frame(peak_anno_merge=c(peak_anno_pos,peak_anno_neg))
diff_peak <- cbind(diff_peak,rbind(diff_peak_pos,diff_peak_neg))

agg_diff_peak <- aggregate(diff_peak[,2:9],by=list(peak_anno_merge=diff_peak$peak_anno_merge),FUN=sum)

agg_diff_peak_log2FC <- data.frame(Lysate_rep1=log2(agg_diff_peak[,2]*0.93/agg_diff_peak[,3]),Lysate_rep2=log2(agg_diff_peak[,4]*1.07/agg_diff_peak[,5]),Result_rep1=log2(agg_diff_peak[,6]*0.78/agg_diff_peak[,7]),Result_rep2=log2(agg_diff_peak[,8]*0.78/agg_diff_peak[,9]))

rownames(agg_diff_peak_log2FC) <- agg_diff_peak[,1]

diff_peak <- t(scale(t(agg_diff_peak_log2FC)))

Heatmap(matrix = diff_peak,show_column_names =TRUE,show_row_names=FALSE,
        cluster_rows = TRUE,cluster_columns = TRUE,
        width = unit(2,'inches'))


## express_gene
countdata_gene <- read.table("express_matrix.tab")


# get meta data
express_matrix_gene <- countdata_gene[,1:4]
meta_data <- data.frame(rep = c("rep1","rep2","rep1","rep2"),treatment = c("Lysate","Lysate","Result","Result"),row.names = sample_name[1:4],stringsAsFactors = TRUE)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = express_matrix_gene,colData = meta_data,design = ~ treatment)

dds <- dds[rowSums(counts(dds)>=10)>=4,]

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

res_select <- res[substr(rownames(res),1,18) %in% rownames(agg_diff_peak_log2FC),]
diff_peak <- diff_peak[rownames(agg_diff_peak_log2FC) %in% substr(rownames(res_select),1,18),]
summary(res_select)

# convert to TPM
# modify annotation
anno <- anno[anno$type == 'gene',]
table(duplicated(anno$gene_id))
rownames(anno) <- as.character(anno$gene_id)
table(rownames(express_matrix_gene) %in% rownames(anno))

tpm_matrix <- do.call(cbind,base::lapply(colnames(express_matrix_gene),FUN = function(x){
  counts <- express_matrix_gene[,x]
  effLen <- anno$width
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  return(exp(rate - denom + log(1e6)))
}))
colnames(tpm_matrix) <- colnames(express_matrix_gene)
rownames(tpm_matrix) <- rownames(express_matrix_gene)
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
Heatmap(matrix = tpm_matrix,show_column_names =TRUE,show_row_names=FALSE,
        cluster_rows = TRUE,cluster_columns = TRUE,
        top_annotation = cell_anno, left_annotation = gene_anno,
        width = unit(2,'inches'))



heatmap_ca <- Heatmap(diff_peak,show_row_dend = FALSE,  
                      #row_km=3,
                      #column_km=2, 
                      #cluster_rows = FALSE,
                      show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                      heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
heatmap_total <- Heatmap(tpm_matrix,show_row_dend = FALSE,  
                         #row_km=3,
                         #column_km=2, 
                         #cluster_rows = FALSE,
                         show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                         heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
ht_list =heatmap_ca+heatmap_total
draw(ht_list,heatmap_legend_side = "top")
