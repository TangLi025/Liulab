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

bed_df <- read.table("~/reference/annotation/mm19/repeats/mm19_repeats_family.bed")

saf <- data.frame(bed_df[,6],bed_df[,c(1,2,3,4)])
colnames(saf)<-c("GeneID","Chr","Start","End","Strand")

write.table(saf,file = "~/reference/annotation/mm19/repeats/mm19_repeats_family.saf",quote=FALSE,sep="\t",row.names = FALSE)



library(data.table)
counts <- fread(file="~/LinLong/13_repeats/repeats_count.tab",nThread = 40)

counts <- counts[,-c(2,3,4,5)]
rownames(counts) <- counts$Geneid
colnames(counts) <- c("Geneid","Length","Lysate_input_rep1","Lysate_input_rep2","Lysate_IP_rep1","Lysate_IP_rep2","Result_input_rep1","Result_input_rep2","Result_IP_rep1","Result_IP_rep2")

countdata <- as.data.frame(counts[,3:10])
rownames(countdata) <- counts$Geneid

# get meta data
express_matrix <- as.data.frame(counts[,c(3,4,7,8)])
rownames(express_matrix) <- counts$Geneid
meta_data <- data.frame(rep = c("rep1","rep2","rep1","rep2"),treatment = c("Lysate","Lysate","Result","Result"),row.names = c("Lysate_input_rep1","Lysate_input_rep2","Result_input_rep1","Result_input_rep2"),stringsAsFactors = TRUE)

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = express_matrix,colData = meta_data,design = ~ treatment)

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
res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange),]

res

res_select <- res[abs(res$log2FoldChange) > 0.3 & res$padj < 0.05,]
summary(res_select)

# convert to TPM
# modify annotation

tpm_matrix <- do.call(cbind,base::lapply(colnames(express_matrix),FUN = function(x){
  count <- express_matrix[,x]
  effLen <- counts$Length
  rate <- log(count) - log(effLen)
  denom <- log(sum(exp(rate)))
  return(exp(rate - denom + log(1e6)))
}))
colnames(tpm_matrix) <- colnames(express_matrix)
rownames(tpm_matrix) <- rownames(express_matrix)
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
Heatmap(matrix = tpm_matrix,show_column_names =TRUE,show_row_names=TRUE,
        cluster_rows = TRUE,cluster_columns = TRUE,
        top_annotation = cell_anno, left_annotation = gene_anno,
        width = unit(2,'inches'))

write.table(countdata,file="express_matrix.tab")
write.table(res,file="res.tab")
write.csv(res_select,file="res_select.csv")

# m6A level
# scale factor
total_reads <- colSums(countdata)

library(edgeR)
myCPM <- cpm(countdata)

### 20220213 先除后加
compare_spike_in <- data.frame(Lysate_rep1=myCPM[,3]/myCPM[,1]*0.94,Lysate_rep2=myCPM[,4]/myCPM[,2]*1.06,Result_rep1=myCPM[,7]/myCPM[,5]*0.69,Result_rep2=myCPM[,8]/myCPM[,6]*0.66)

compare_mean <- data.frame(Lysate=(compare_spike_in[,1]+compare_spike_in[,2])/2,Result=((compare_spike_in[,3]+compare_spike_in[,4])/2))
compare_log <- data.frame(family=counts$Geneid,log2_Result_Lysate=log2(compare_mean[,2]/compare_mean[,1]))
compare_log <- as.matrix(compare_log)

compare_log <- compare_log[which(is.finite(compare_log[,2] )),]
compare_log <- compare_log[order(compare_log[,2]),]

x1 = factor(compare_log$family, levels=compare_log$family)
ggplot(data=compare_log)+
  geom_col(mapping=aes(x=x1,y=log2_Result_Lysate,fill=log2_Result_Lysate))+
  theme(axis.text.x = element_text(size = 12,angle = 45,hjust = 1))
###

### 20220213 先加后除
compare_spike_in <- data.frame(Lysate=(myCPM[,3]+myCPM[,4])/(myCPM[,1]+myCPM[,2]),Result=(myCPM[,7]+myCPM[,8])*0.68/(myCPM[,5]+myCPM[,6]))

compare_log <- data.frame(family=counts$Geneid,log2_Result_Lysate=log2(compare_spike_in[,2]/compare_spike_in[,1]))

compare_log <- compare_log[which(is.finite(compare_log[,2] )),]
compare_log <- compare_log[order(compare_log[,2]),]

x1 = factor(compare_log$family, levels=compare_log$family)
ggplot(data=compare_log)+
  geom_col(mapping=aes(x=x1,y=log2_Result_Lysate,fill=log2_Result_Lysate))+
  theme(axis.text.x = element_text(size = 12,angle = 45,hjust = 1))

compare <- data.frame(Lysate_rep1=log2(myCPM[,3]/myCPM[,1]),Lysate_rep2=log2(myCPM[,4]/myCPM[,2]),Result_rep1=log2(myCPM[,7]/myCPM[,5]),Result_rep2=log2(myCPM[,8]/myCPM[,6]))

select_m6A <- as.matrix(compare_spike_in[rownames(compare_spike_in) %in% rownames(tpm_matrix),])




select_m6A <- t(scale(t(select_m6A)))
summary(is.na(select_m6A))
summary(is.nan(select_m6A))
summary(is.infinite(select_m6A))

select_m6A <- select_m6A[which(rowSums(!is.na(select_m6A))==4),]
tpm_matrix <- tpm_matrix[which(rowSums(!is.na(select_m6A))==4),]

Heatmap(matrix = select_m6A,show_column_names =TRUE,show_row_names=TRUE,
        cluster_rows = TRUE,cluster_columns = TRUE,
        width = unit(2,'inches'))


compare_scale <- t(scale(t(compare)))
summary(is.na(compare_scale))
summary(is.nan(compare_scale))
summary(is.infinite(compare_scale))
compare_scale <- compare_scale[which(rowSums(!is.na(compare_scale))==4),]

Heatmap(matrix = compare_scale,show_column_names =TRUE,show_row_names=TRUE,
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