library(Rsubread)

setwd("/disk/user_08/Data/TC1-planB/05_bam_change_index_hisat2/01_bam_sorted/")
path <- "/disk/user_08/Data/TC1-planB/05_bam_change_index_hisat2/01_bam_sorted/11_bam_merge/"
sample_name <- c("p0_input_rep1","p0_input_rep2","p5_input_rep1","p5_input_rep2",
                 "p10_input_rep1","p10_input_rep2","rp2_input_rep1","rp2_input_rep2",
                 "p0_ip_rep1","p0_ip_rep2","p5_ip_rep1","p5_ip_rep2",
                 "p10_ip_rep1","p10_ip_rep2","rp2_ip_rep1","rp2_ip_rep2")

caRNA_count_trans_df <- read.table("./11_bam_merge/05_gene_expression_tl/caRNA_trans_expression.tab")
colnames(caRNA_count_trans_df) <- sample_name

nor_ratio <- c(41608064,68341489,
               #               70516802,68367198,
               57358338,70205279,
               #               70496756,70107156,
               67542625,52345318,
               #               75806743,76023704,
               63139881,69162349)
#               78116866,66711470)
nor_ratio <- mean(nor_ratio)/nor_ratio

express_p0p10_nor <- caRNA_count_trans_df[,c(1,2,5,6,9,10,13,14)]
for ( i in 1:8){
  express_p0p10_nor[,i] <- express_p0p10_nor[,i]*nor_ratio[i]
}

express_p0p10_nor  <- express_p0p10_nor[which(apply(express_p0p10_nor>20,1,sum)>=4),]

scale_factor <- c(1,1.232340796)

m6A_level_trans_log2FC <- data.frame(trans_id =substr(rownames(express_p0p10_nor),1,18))

for (i in 1:2){
  m6A_level_trans_log2FC[,i+1]=log2((express_p0p10_nor[,2*i+3]+express_p0p10_nor[,2*i+4]+0.001)/(express_p0p10_nor[,2*i-1]+express_p0p10_nor[,2*i]+0.001))
}

colnames(m6A_level_trans_log2FC) <- c("trans_id","p0_m6A","p10_m6A")
rownames(m6A_level_trans_log2FC) <- m6A_level_trans_log2FC$trans_id

## PDUI

dapars2 <- read.table("~/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/DaPars2_caRNA_mRNA_mean_omitNA.txt",
                      header = T,row.names = 1)
dapars2_ca_scale <- t(scale(t(dapars2[,1:4])))
dapars2_m_scale <- t(scale(t(dapars2[,5:8])))
dapars2_scale <- na.omit(cbind(dapars2_ca_scale,dapars2_m_scale))
dapars2_ca_scale <- dapars2_scale[,1:4]
dapars2_m_scale <- dapars2_scale[,5:8]

library(ComplexHeatmap)
heatmap_PDUI_ca <- Heatmap(dapars2_ca_scale,show_row_dend = FALSE,  
                        #row_km=3,
                        #column_km=2, 
                        #cluster_rows = FALSE,
                        cluster_columns = T,
                        show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                        heatmap_legend_param=list(title = "",direction = "vertical", legend_width = unit(5, "cm")))

heatmap_PDUI_m <- Heatmap(dapars2_m_scale,show_row_dend = FALSE,  
        #row_km=3,
        #column_km=2, 
        #cluster_rows = FALSE,
        cluster_columns = T,
        show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
        heatmap_legend_param=list(title = "",direction = "vertical", legend_width = unit(5, "cm")))

ht_list =heatmap_PDUI_ca+heatmap_PDUI_m
draw(ht_list,heatmap_legend_side = "top")

ht_list =heatmap_PDUI_m+heatmap_PDUI_ca
draw(ht_list,heatmap_legend_side = "top")


### caRNA expression level

caRNA_express_counts <- read.table("~/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/11_bam_merge/05_gene_expression_tl/caRNA_trans_expression.tab")
colnames(caRNA_express_counts) <- sample_name
caRNA_input_counts_p0p10 <- caRNA_express_counts[,c(1,2,5,6)]

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

## volcano
### PDUI p0p10
dapars_ca_p0p10 <- read.table("/disk/user_09/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p0p10.txt",
                               header = T,row.names = 1)

dapars_caRNA_p0p10_change <- dapars_ca_p0p10[which(dapars_ca_p0p10$filter != "NC"),]

dapars_caRNA_p0p10_change_hm <- as.matrix(dapars_caRNA_p0p10_change[,c(1,2)])
rownames(dapars_caRNA_p0p10_change_hm) <- substr(rownames(dapars_caRNA_p0p10_change_hm),1,18)

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
  
# heatmap
library(ComplexHeatmap)

heatmap_PDUI <- Heatmap(t(scale(t(dapars_caRNA_p0p10_change_hm))),show_row_dend = FALSE,  
                        #row_km=3,
                        #column_km=2, 
                        #cluster_rows = FALSE,
                        show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                        heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
heatmap_m6A <- Heatmap(t(scale(t(m6A_dapars_change))),show_row_dend = FALSE,  
                       #row_km=3,
                       #column_km=2, 
                       #cluster_rows = FALSE,
                       show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                       heatmap_legend_param=list(title = "",direction = "vertical", legend_width = unit(5, "cm")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
ht_list =heatmap_m6A+heatmap_PDUI
draw(ht_list,heatmap_legend_side = "top")

# convert to TPM
# modify annotation
anno <- anno[anno$type == 'transcript',]
table(duplicated(anno$transcript_id))
rownames(anno) <- as.character(anno$transcript_id)
table(rownames(caRNA_input_counts) %in% rownames(anno))

tpm_matrix <- do.call(cbind,base::lapply(colnames(caRNA_input_counts),FUN = function(x){
  counts <- caRNA_input_counts[,x]
  effLen <- anno$width
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  return(exp(rate - denom + log(1e6)))
}))
colnames(tpm_matrix) <- colnames(caRNA_input_counts)
rownames(tpm_matrix) <- rownames(caRNA_input_counts)
table(colSums(tpm_matrix))

# heatmap
library(ComplexHeatmap)
tpm_matrix <- tpm_matrix[rownames(res_select),]
tpm_matrix <- t(scale(t(tpm_matrix)))

# create annatation
cell_anno <- HeatmapAnnotation(treat = meta_data[colnames(tpm_matrix),"stage"])

res_select$express <- NA
res_select[res_select$log2FoldChange >= 0,"express"] <- "UP"
res_select[res_select$log2FoldChange < 0,"express"] <- "DOWN"

gene_anno <- HeatmapAnnotation(express = res_select[rownames(tpm_matrix),"express"],which = 'row')
Heatmap(matrix = tpm_matrix,show_column_names =TRUE,show_row_names=FALSE,
        cluster_rows = TRUE,cluster_columns = TRUE,
        top_annotation = cell_anno, left_annotation = gene_anno,
        width = unit(2,'inches'))

write.table(countdata,file="express_matrix.tab")
write.table(res,file="res.tab")
write.csv(res_select,file="res_select.csv")


peak_anno_map <- data.frame(peak_Entrez.ID=peak_annotation$Entrez.ID[peak_annotation$Entrez.ID %in% substr(rownames(res_select),1,18)])

express_matrix_nor <- express_matrix_nor[peak_annotation$Entrez.ID %in% substr(rownames(res_select),1,18),]

peak <- as.data.frame(cbind(peak_anno_map,express_matrix_nor))

agg_peak <- aggregate(peak[,2:9],by=list(peak$peak_Entrez.ID),FUN=sum)

agg_peak_log2FC <- data.frame(Lysate_rep1=log2(agg_peak[,2]*0.93/agg_peak[,3]),Lysate_rep2=log2(agg_peak[,4]*1.07/agg_peak[,5]),Result_rep1=log2(agg_peak[,6]*0.78/agg_peak[,7]),Result_rep2=log2(agg_peak[,8]*0.78/agg_peak[,9]))

rownames(agg_peak_log2FC) <- agg_peak[,1]

diff_peak <- t(scale(t(agg_peak_log2FC)))

diff_peak <- diff_peak[rowSums(is.finite(diff_peak))==4,]
tpm_matrix <- tpm_matrix[substr(rownames(tpm_matrix),1,18) %in% rownames(diff_peak),]


Heatmap(matrix = diff_peak,show_column_names =TRUE,show_row_names=FALSE,
        cluster_rows = TRUE,cluster_columns = TRUE,
        width = unit(2,'inches'))

heatmap_ca <- Heatmap(tpm_matrix,show_row_dend = FALSE,  
                      #row_km=3,
                      #column_km=2, 
                      #cluster_rows = FALSE,
                      show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                      heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
heatmap_total <- Heatmap(diff_peak,show_row_dend = FALSE,  
                         #row_km=3,
                         #column_km=2, 
                         #cluster_rows = FALSE,
                         show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                         heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
ht_list =heatmap_ca+heatmap_total
draw(ht_list,heatmap_legend_side = "top")







mRNA_express_counts <- read.table("~/user_08_TC1/10_mRNA/05_bam_hisat2/01_bam_sorted/11_bam_merge/mRNA_gene_expression_matrix.txt")

dapars_caRNA_p0p10_change <- dapars_caRNA_p0p10[which(dapars_caRNA_p0p10$filter != "NC"),]

dapars_caRNA_p0p10_change_hm <- as.matrix(dapars_caRNA_p0p10_change[,c(2,3)])
rownames(dapars_caRNA_p0p10_change_hm) <- substr(dapars_caRNA_p0p10_change$Gene,1,18)

m6A_dapars_change <- m6A_level_trans_log2FC[,c(2:3)]
m6A_dapars_change <- m6A_dapars_change[rownames(dapars_caRNA_p0p10_change_hm),]

# heatmap
library(ComplexHeatmap)

heatmap_PDUI <- Heatmap(t(scale(t(dapars_caRNA_p0p10_change_hm))),show_row_dend = FALSE,  
                        #row_km=3,
                        #column_km=2, 
                        #cluster_rows = FALSE,
                        show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                        heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
heatmap_m6A <- Heatmap(t(scale(t(m6A_dapars_change))),show_row_dend = FALSE,  
                       #row_km=3,
                       #column_km=2, 
                       #cluster_rows = FALSE,
                       show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                       heatmap_legend_param=list(title = "",direction = "vertical", legend_width = unit(5, "cm")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.5), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
ht_list =heatmap_m6A+heatmap_PDUI
draw(ht_list,heatmap_legend_side = "top")


### 3UTR p0p10
caRNA_count_3UTR <- featureCounts(files=paste0(path,sample_name,".bam"),
                                  annot.ext = "/disk1/home/user_09/reference/annotation/mm19/dapars/mm39_gencode_dapars_longest_protein.saf",
                                  isGTFAnnotationFile = FALSE, 
                                  isPairedEnd=TRUE,
                                  requireBothEndsMapped=TRUE,
                                  allowMultiOverlap = TRUE,
                                  fracOverlap = 0.5,
                                  countChimericFragments=FALSE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific = 2,checkFragLength = TRUE,
                                  minFragLength = 40, maxFragLength = 2000, 
                                  verbose = FALSE, nthreads=40)

caRNA_count_3UTR_df <- caRNA_count_3UTR$counts
colnames(caRNA_count_3UTR_df) <- sample_name

apply(caRNA_count_3UTR_df,2,sum)

write.table(caRNA_count_3UTR_df,"./11_bam_merge/05_gene_expression_tl/caRNA_3UTR_expression.tab")

#nor_ratio <- c(41608064,68341489,70516802,68367198,57358338,70205279,70496756,70107156,
#               67542625,52345318,75806743,76023704,63139881,69162349,78116866,66711470)
#nor_ratio <- mean(nor_ratio)/nor_ratio

express_p0p10_3UTR_nor <- caRNA_count_3UTR_df[,c(1,2,5,6,9,10,13,14)]
for ( i in 1:8){
  express_p0p10_3UTR_nor[,i] <- express_p0p10_3UTR_nor[,i]*nor_ratio[i]
}

express_p0p10_3UTR_nor  <- express_p0p10_3UTR_nor[which(apply(express_p0p10_3UTR_nor>10,1,sum)>=4),]
# scale_factor <- c(0.995575897,1.004443762,
#                   #                  1.466860346,1.637934483,
#                   1.250982735,1.213976657)
# #                  0.954841898,0.98519686)
scale_factor <- c(1,1.232340796)

m6A_level_3UTR_log2FC <- data.frame(trans_id =substr(rownames(express_p0p10_3UTR_nor),1,18))
for (i in 1:2){
  m6A_level_3UTR_log2FC[,i+1]=log2((express_p0p10_3UTR_nor[,2*i+3]+express_p0p10_3UTR_nor[,2*i+4])*scale_factor[i]/(express_p0p10_3UTR_nor[,2*i-1]+express_p0p10_3UTR_nor[,2*i]))
}
# for (i in 1:2){
#   m6A_level_3UTR_log2FC[,i+1]=log2(((express_p0p10_3UTR_nor[,2*i+3]+express_p0p10_3UTR_nor[,2*i+4])+0.001)*scale_factor[i]/((express_p0p10_3UTR_nor[,2*i-1]+express_p0p10_3UTR_nor[,2*i]+0.001)))
# }
#colnames(m6A_level_3UTR_log2FC) <- c("trans_id","p0_rep1","p0_rep2","p10_rep1","p10_rep2")
colnames(m6A_level_3UTR_log2FC) <- c("trans_id","p0_m6A","p10_m6A")
rownames(m6A_level_3UTR_log2FC) <- m6A_level_3UTR_log2FC$trans_id

m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[,c(2:3)]
m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[which(((m6A_level_3UTR_log2FC$p0_m6A>0) + (m6A_level_3UTR_log2FC$p10_m6A>0))>=1),]


## PDUI
dapars_caRNA_p0p10_change <- dapars_caRNA_p0p10[which(dapars_caRNA_p0p10$filter != "NC"),]

dapars_caRNA_p0p10_change_hm <- as.matrix(dapars_caRNA_p0p10_change[,c(2,3)])
rownames(dapars_caRNA_p0p10_change_hm) <- substr(dapars_caRNA_p0p10_change$Gene,1,18)

m6A_dapars_3UTR_change <- m6A_level_3UTR_log2FC[rownames(dapars_caRNA_p0p10_change_hm),]

m6A_dapars_3UTR_change <- m6A_dapars_3UTR_change[which(apply(apply(m6A_dapars_3UTR_change,2,is.finite),1,sum)==2),]
dapars_caRNA_p0p10_change_3UTR_hm <- dapars_caRNA_p0p10_change_hm[rownames(m6A_dapars_3UTR_change),]

# heatmap
library(ComplexHeatmap)

heatmap_PDUI <- Heatmap(t(scale(t(dapars_caRNA_p0p10_change_3UTR_hm))),show_row_dend = FALSE,  
                        #row_km=3,
                        #column_km=2, 
                        #cluster_rows = FALSE,
                        show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                        heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")),
                        col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#B210FF", colorRampPalette(c("#B210FF", "white", "#EECE13"))(6), "#EECE13")))
#col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#81FFEF", colorRampPalette(c("#81FFEF", "white", "#F067B4"))(6), "#F067B4")))
#col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.5), 1), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(3), "#FF0000")))
heatmap_m6A <- Heatmap(t(scale(t(m6A_dapars_3UTR_change))),show_row_dend = FALSE,  
                       #row_km=3,
                       #column_km=2, 
                       #cluster_rows = FALSE,
                       show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                       heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")),
                       col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#B210FF", colorRampPalette(c("#B210FF", "white", "#EECE13"))(6), "#EECE13")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.2), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
ht_list =heatmap_m6A+heatmap_PDUI
draw(ht_list,heatmap_legend_side = "top")

ht_list =heatmap_PDUI+heatmap_m6A
draw(ht_list,heatmap_legend_side = "top")

### 3UTR p0p5

express_p0p5_3UTR_nor <- caRNA_count_3UTR_df[,c(1,2,3,4,9,10,11,12)]

nor_ratio_p0p5 <- c(41608064,68341489,
                    70516802,68367198,
                    #57358338,70205279,
                    #               70496756,70107156,
                    67542625,52345318,
                    75806743,76023704)
#63139881,69162349)
#               78116866,66711470)
nor_ratio_p0p5 <- mean(nor_ratio_p0p5)/nor_ratio_p0p5
for ( i in 1:8){
  express_p0p5_3UTR_nor[,i] <- express_p0p5_3UTR_nor[,i]*nor_ratio_p0p5[i]
}

express_p0p5_3UTR_nor  <- express_p0p5_3UTR_nor[which(apply(express_p0p5_3UTR_nor>10,1,sum)>=4),]
# scale_factor <- c(0.995575897,1.004443762,
#                   #                  1.466860346,1.637934483,
#                   1.250982735,1.213976657)
# #                  0.954841898,0.98519686)
scale_factor_p0p5 <- c(1,1.550039078)

m6A_level_3UTR_log2FC <- data.frame(trans_id =substr(rownames(express_p0p5_3UTR_nor),1,18))
# for (i in 1:2){
#   m6A_level_3UTR_log2FC[,i+1]=log2((express_p0p5_3UTR_nor[,2*i+3]+express_p0p5_3UTR_nor[,2*i+4])*scale_factor[i]/(express_p0p5_3UTR_nor[,2*i-1]+express_p0p5_3UTR_nor[,2*i]))
# }
for (i in 1:2){
  m6A_level_3UTR_log2FC[,i+1]=log2(((express_p0p5_3UTR_nor[,2*i+3]+express_p0p5_3UTR_nor[,2*i+4])+0.001)/((express_p0p5_3UTR_nor[,2*i-1]+express_p0p5_3UTR_nor[,2*i]+0.001)))
}
#colnames(m6A_level_3UTR_log2FC) <- c("trans_id","p0_rep1","p0_rep2","p10_rep1","p10_rep2")
colnames(m6A_level_3UTR_log2FC) <- c("trans_id","p0_m6A","p5_m6A")
rownames(m6A_level_3UTR_log2FC) <- m6A_level_3UTR_log2FC$trans_id

m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[,c(2:3)]
m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[which(((m6A_level_3UTR_log2FC$p0_m6A>0) + (m6A_level_3UTR_log2FC$p5_m6A>0))>=1),]


## PDUI
dapars_caRNA_p0p5_change <- dapars_caRNA_p0p5[which(dapars_caRNA_p0p5$filter != "NC"),]

dapars_caRNA_p0p5_change_hm <- as.matrix(dapars_caRNA_p0p5_change[,c(2,3)])
rownames(dapars_caRNA_p0p5_change_hm) <- substr(dapars_caRNA_p0p5_change$Gene,1,18)

m6A_dapars_3UTR_change <- m6A_level_3UTR_log2FC[rownames(dapars_caRNA_p0p5_change_hm),]

m6A_dapars_3UTR_change <- m6A_dapars_3UTR_change[which(apply(apply(m6A_dapars_3UTR_change,2,is.finite),1,sum)==2),]
dapars_caRNA_p0p5_change_3UTR_hm <- dapars_caRNA_p0p5_change_hm[rownames(m6A_dapars_3UTR_change),]

# heatmap
library(ComplexHeatmap)

heatmap_PDUI <- Heatmap(t(scale(t(dapars_caRNA_p0p5_change_3UTR_hm))),show_row_dend = FALSE,  
                        #row_km=3,
                        #column_km=2, 
                        #cluster_rows = FALSE,
                        show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                        heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")),
                        col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#B210FF", colorRampPalette(c("#B210FF", "white", "#EECE13"))(6), "#EECE13")))
#col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#81FFEF", colorRampPalette(c("#81FFEF", "white", "#F067B4"))(6), "#F067B4")))
#col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.5), 1), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(3), "#FF0000")))
heatmap_m6A <- Heatmap(t(scale(t(m6A_dapars_3UTR_change))),show_row_dend = FALSE,  
                       #row_km=3,
                       #column_km=2, 
                       #cluster_rows = FALSE,
                       show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                       heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")),
                       col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#B210FF", colorRampPalette(c("#B210FF", "white", "#EECE13"))(6), "#EECE13")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.2), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
ht_list =heatmap_m6A+heatmap_PDUI
draw(ht_list,heatmap_legend_side = "top")

ht_list =heatmap_PDUI+heatmap_m6A
draw(ht_list,heatmap_legend_side = "top")



### 3UTR p0rp2

express_p0rp2_3UTR_nor <- caRNA_count_3UTR_df[,c(1,2,7,8,9,10,15,16)]

nor_ratio_p0rp2 <- c(41608064,68341489,
                     #70516802,68367198,
                     #57358338,70205279,
                     70496756,70107156,
                     67542625,52345318,
                     #75806743,76023704)
                     #63139881,69162349
                     78116866,66711470)
nor_ratio_p0rp2 <- mean(nor_ratio_p0rp2)/nor_ratio_p0rp2
for ( i in 1:8){
  express_p0rp2_3UTR_nor[,i] <- express_p0rp2_3UTR_nor[,i]*nor_ratio_p0rp2[i]
}

express_p0rp2_3UTR_nor  <- express_p0rp2_3UTR_nor[which(apply(express_p0rp2_3UTR_nor>10,1,sum)>=4),]
# scale_factor <- c(0.995575897,1.004443762,
#                   #                  1.466860346,1.637934483,
#                   1.250982735,1.213976657)
# #                  0.954841898,0.98519686)
scale_factor_p0rp2 <- c(1,0.97)

m6A_level_3UTR_log2FC <- data.frame(trans_id =substr(rownames(express_p0rp2_3UTR_nor),1,18))
# for (i in 1:2){
#   m6A_level_3UTR_log2FC[,i+1]=log2((express_p0rp2_3UTR_nor[,2*i+3]+express_p0rp2_3UTR_nor[,2*i+4])*scale_factor[i]/(express_p0rp2_3UTR_nor[,2*i-1]+express_p0rp2_3UTR_nor[,2*i]))
# }
for (i in 1:2){
  m6A_level_3UTR_log2FC[,i+1]=log2(((express_p0rp2_3UTR_nor[,2*i+3]+express_p0rp2_3UTR_nor[,2*i+4])+0.001)/((express_p0rp2_3UTR_nor[,2*i-1]+express_p0rp2_3UTR_nor[,2*i]+0.001)))
}
#colnames(m6A_level_3UTR_log2FC) <- c("trans_id","p0_rep1","p0_rep2","p10_rep1","p10_rep2")
colnames(m6A_level_3UTR_log2FC) <- c("trans_id","p0_m6A","rp2_m6A")
rownames(m6A_level_3UTR_log2FC) <- m6A_level_3UTR_log2FC$trans_id

m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[,c(2:3)]
m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[which(((m6A_level_3UTR_log2FC$p0_m6A>0) + (m6A_level_3UTR_log2FC$rp2_m6A>0))>=1),]


## PDUI
dapars_caRNA_p0rp2_change <- dapars_caRNA_p0rp2[which(dapars_caRNA_p0rp2$filter != "NC"),]

dapars_caRNA_p0rp2_change_hm <- as.matrix(dapars_caRNA_p0rp2_change[,c(2,3)])
rownames(dapars_caRNA_p0rp2_change_hm) <- substr(dapars_caRNA_p0rp2_change$Gene,1,18)

m6A_dapars_3UTR_change <- m6A_level_3UTR_log2FC[rownames(dapars_caRNA_p0rp2_change_hm),]

m6A_dapars_3UTR_change <- m6A_dapars_3UTR_change[which(apply(apply(m6A_dapars_3UTR_change,2,is.finite),1,sum)==2),]
dapars_caRNA_p0rp2_change_3UTR_hm <- dapars_caRNA_p0rp2_change_hm[rownames(m6A_dapars_3UTR_change),]

# heatmap
library(ComplexHeatmap)

heatmap_PDUI <- Heatmap(t(scale(t(dapars_caRNA_p0rp2_change_3UTR_hm))),show_row_dend = FALSE,  
                        #row_km=3,
                        #column_km=2, 
                        #cluster_rows = FALSE,
                        show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                        heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")),
                        col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#B210FF", colorRampPalette(c("#B210FF", "white", "#EECE13"))(6), "#EECE13")))
#col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#81FFEF", colorRampPalette(c("#81FFEF", "white", "#F067B4"))(6), "#F067B4")))
#col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.5), 1), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(3), "#FF0000")))
heatmap_m6A <- Heatmap(t(scale(t(m6A_dapars_3UTR_change))),show_row_dend = FALSE,  
                       #row_km=3,
                       #column_km=2, 
                       #cluster_rows = FALSE,
                       show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                       heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")),
                       col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#B210FF", colorRampPalette(c("#B210FF", "white", "#EECE13"))(6), "#EECE13")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.2), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
ht_list =heatmap_m6A+heatmap_PDUI
draw(ht_list,heatmap_legend_side = "top")

ht_list =heatmap_PDUI+heatmap_m6A
draw(ht_list,heatmap_legend_side = "top")


### 3UTR p10rp2

express_p10rp2_3UTR_nor <- caRNA_count_3UTR_df[,c(5,6,7,8,13,14,15,16)]

nor_ratio_p10rp2 <- c(#41608064,68341489,
  #70516802,68367198,
  57358338,70205279,
  70496756,70107156,
  #67542625,52345318,
  #75806743,76023704)
  63139881,69162349,
  78116866,66711470)
nor_ratio_p10rp2 <- mean(nor_ratio_p10rp2)/nor_ratio_p10rp2
for ( i in 1:8){
  express_p10rp2_3UTR_nor[,i] <- express_p10rp2_3UTR_nor[,i]*nor_ratio_p10rp2[i]
}

express_p10rp2_3UTR_nor  <- express_p10rp2_3UTR_nor[which(apply(express_p10rp2_3UTR_nor>10,1,sum)>=4),]
# scale_factor <- c(0.995575897,1.004443762,
#                   #                  1.466860346,1.637934483,
#                   1.250982735,1.213976657)
# #                  0.954841898,0.98519686)
scale_factor_p10rp2 <- c(1.23,0.97)

m6A_level_3UTR_log2FC <- data.frame(trans_id =substr(rownames(express_p10rp2_3UTR_nor),1,18))
for (i in 1:2){
  m6A_level_3UTR_log2FC[,i+1]=log2(((express_p10rp2_3UTR_nor[,2*i+3]+express_p10rp2_3UTR_nor[,2*i+4])*scale_factor[i]+0.001)/(express_p10rp2_3UTR_nor[,2*i-1]+express_p10rp2_3UTR_nor[,2*i]+0.001))
}
# for (i in 1:2){
#   m6A_level_3UTR_log2FC[,i+1]=log2(((express_p10rp2_3UTR_nor[,2*i+3]+express_p10rp2_3UTR_nor[,2*i+4])+0.001)/((express_p10rp2_3UTR_nor[,2*i-1]+express_p10rp2_3UTR_nor[,2*i]+0.001)))
# }
#colnames(m6A_level_3UTR_log2FC) <- c("trans_id","p10_rep1","p10_rep2","p10_rep1","p10_rep2")
colnames(m6A_level_3UTR_log2FC) <- c("trans_id","p10_m6A","rp2_m6A")
rownames(m6A_level_3UTR_log2FC) <- m6A_level_3UTR_log2FC$trans_id

m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[,c(2:3)]
m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[which(((m6A_level_3UTR_log2FC$p10_m6A>0) + (m6A_level_3UTR_log2FC$rp2_m6A>0))>=1),]


## PDUI
dapars_caRNA_p10rp2_change <- dapars_caRNA_p10rp2[which(dapars_caRNA_p10rp2$filter != "NC"),]

dapars_caRNA_p10rp2_change_hm <- as.matrix(dapars_caRNA_p10rp2_change[,c(2,3)])
rownames(dapars_caRNA_p10rp2_change_hm) <- substr(dapars_caRNA_p10rp2_change$Gene,1,18)

m6A_dapars_3UTR_change <- m6A_level_3UTR_log2FC[rownames(dapars_caRNA_p10rp2_change_hm),]

m6A_dapars_3UTR_change <- m6A_dapars_3UTR_change[which(apply(apply(m6A_dapars_3UTR_change,2,is.finite),1,sum)==2),]
dapars_caRNA_p10rp2_change_3UTR_hm <- dapars_caRNA_p10rp2_change_hm[rownames(m6A_dapars_3UTR_change),]

# heatmap
library(ComplexHeatmap)

heatmap_PDUI <- Heatmap(t(scale(t(dapars_caRNA_p10rp2_change_3UTR_hm))),show_row_dend = FALSE,  
                        #row_km=3,
                        #column_km=2, 
                        #cluster_rows = FALSE,
                        show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                        heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")),
                        col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#B210FF", colorRampPalette(c("#B210FF", "white", "#EECE13"))(6), "#EECE13")))
#col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#81FFEF", colorRampPalette(c("#81FFEF", "white", "#F067B4"))(6), "#F067B4")))
#col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.5), 1), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(3), "#FF0000")))
heatmap_m6A <- Heatmap(t(scale(t(m6A_dapars_3UTR_change))),show_row_dend = FALSE,  
                       #row_km=3,
                       #column_km=2, 
                       #cluster_rows = FALSE,
                       show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                       heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")),
                       col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#B210FF", colorRampPalette(c("#B210FF", "white", "#EECE13"))(6), "#EECE13")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.2), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
ht_list =heatmap_m6A+heatmap_PDUI
draw(ht_list,heatmap_legend_side = "top")

ht_list =heatmap_PDUI+heatmap_m6A
draw(ht_list,heatmap_legend_side = "top")


### 5UTR p0p10
caRNA_count_5UTR <- featureCounts(files=paste0(path,sample_name,".bam"),
                                  annot.ext = "/disk1/home/user_09/reference/annotation/mm19/dapars/mm39_gencode_5UTR_longest_protein.saf",
                                  isGTFAnnotationFile = FALSE, 
                                  isPairedEnd=TRUE,
                                  requireBothEndsMapped=TRUE,
                                  allowMultiOverlap = TRUE,
                                  fracOverlap = 0.5,
                                  countChimericFragments=FALSE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific = 2,checkFragLength = TRUE,
                                  minFragLength = 40, maxFragLength = 2000, 
                                  verbose = FALSE, nthreads=40)

caRNA_count_5UTR_df <- caRNA_count_5UTR$counts
colnames(caRNA_count_5UTR_df) <- sample_name

write.table(caRNA_count_5UTR_df,"./11_bam_merge/05_gene_expression_tl/caRNA_5UTR_expression.tab")

nor_ratio <- c(41608064,68341489,
               #70516802,68367198,
               57358338,70205279,
               #70496756,70107156,
               67542625,52345318,
               #75806743,76023704,
               63139881,69162349)
#78116866,66711470)
nor_ratio <- mean(nor_ratio)/nor_ratio

express_p0p10_5UTR_nor <- caRNA_count_5UTR_df[,c(1,2,5,6,9,10,13,14)]
for ( i in 1:8){
  express_p0p10_5UTR_nor[,i] <- express_p0p10_5UTR_nor[,i]*nor_ratio[i]
}

express_p0p10_5UTR_nor  <- express_p0p10_5UTR_nor[which(apply(express_p0p10_5UTR_nor>10,1,sum)>=4),]
# scale_factor <- c(0.995575897,1.004443762,
#                   #                  1.466860346,1.637934483,
#                   1.250982735,1.213976657)
# #                  0.954841898,0.98519686)
scale_factor <- c(1,1.232340796)

m6A_level_5UTR_log2FC <- data.frame(trans_id =substr(rownames(express_p0p10_5UTR_nor),1,18))
# for (i in 1:2){
#   m6A_level_5UTR_log2FC[,i+1]=log2((express_p0p10_5UTR_nor[,2*i+3]+express_p0p10_5UTR_nor[,2*i+4])*scale_factor[i]/(express_p0p10_5UTR_nor[,2*i-1]+express_p0p10_5UTR_nor[,2*i]))
# }
for (i in 1:2){
  m6A_level_5UTR_log2FC[,i+1]=log2(((express_p0p10_5UTR_nor[,2*i+3]+express_p0p10_5UTR_nor[,2*i+4])*scale_factor[i]+0.001)/((express_p0p10_5UTR_nor[,2*i-1]+express_p0p10_5UTR_nor[,2*i]+0.001)))
}
#colnames(m6A_level_5UTR_log2FC) <- c("trans_id","p0_rep1","p0_rep2","p10_rep1","p10_rep2")
colnames(m6A_level_5UTR_log2FC) <- c("trans_id","p0_m6A","p10_m6A")
rownames(m6A_level_5UTR_log2FC) <- m6A_level_5UTR_log2FC$trans_id

m6A_level_5UTR_log2FC <- m6A_level_5UTR_log2FC[,c(2:3)]
m6A_level_5UTR_log2FC <- m6A_level_5UTR_log2FC[which(((m6A_level_5UTR_log2FC$p0_m6A>0) + (m6A_level_5UTR_log2FC$p10_m6A>0))>=1),]


## PDUI
dapars_caRNA_p0p10_change <- dapars_caRNA_p0p10[which(dapars_caRNA_p0p10$filter != "NC"),]

dapars_caRNA_p0p10_change_hm <- as.matrix(dapars_caRNA_p0p10_change[,c(2,3)])
rownames(dapars_caRNA_p0p10_change_hm) <- substr(dapars_caRNA_p0p10_change$Gene,1,18)

m6A_dapars_5UTR_change <- m6A_level_5UTR_log2FC[rownames(dapars_caRNA_p0p10_change_hm),]

m6A_dapars_5UTR_change <- m6A_dapars_5UTR_change[which(apply(apply(m6A_dapars_5UTR_change,2,is.finite),1,sum)==2),]
dapars_caRNA_p0p10_change_5UTR_hm <- dapars_caRNA_p0p10_change_hm[rownames(m6A_dapars_5UTR_change),]

# heatmap
library(ComplexHeatmap)

heatmap_PDUI <- Heatmap(t(scale(t(dapars_caRNA_p0p10_change_5UTR_hm))),show_row_dend = FALSE,  
                        #row_km=3,
                        #column_km=2, 
                        #cluster_rows = FALSE,
                        show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                        heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")),
                        col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#B210FF", colorRampPalette(c("#B210FF", "white", "#EECE13"))(6), "#EECE13")))
#col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#81FFEF", colorRampPalette(c("#81FFEF", "white", "#F067B4"))(6), "#F067B4")))
#col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.5), 1), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(3), "#FF0000")))
heatmap_m6A <- Heatmap(t(scale(t(m6A_dapars_5UTR_change))),show_row_dend = FALSE,  
                       #row_km=3,
                       #column_km=2, 
                       #cluster_rows = FALSE,
                       show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                       heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")),
                       col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#B210FF", colorRampPalette(c("#B210FF", "white", "#EECE13"))(6), "#EECE13")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.2), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
ht_list =heatmap_m6A+heatmap_PDUI
draw(ht_list,heatmap_legend_side = "top")

ht_list =heatmap_PDUI+heatmap_m6A
draw(ht_list,heatmap_legend_side = "top")



### CDS p0p10
caRNA_count_cds <- featureCounts(files=paste0(path,sample_name,".bam"),
                                 annot.ext = "/disk1/home/user_09/reference/annotation/mm19/dapars/mm39_gencode_cds_longest_protein.saf",
                                 isGTFAnnotationFile = FALSE, 
                                 isPairedEnd=TRUE,
                                 requireBothEndsMapped=TRUE,
                                 allowMultiOverlap = TRUE,
                                 fracOverlap = 0.5,
                                 countChimericFragments=FALSE,
                                 countMultiMappingReads = FALSE,
                                 strandSpecific = 2,checkFragLength = TRUE,
                                 minFragLength = 40, maxFragLength = 2000, 
                                 verbose = FALSE, nthreads=40)

caRNA_count_cds_df <- caRNA_count_cds$counts
colnames(caRNA_count_cds_df) <- sample_name

write.table(caRNA_count_cds_df,"./11_bam_merge/05_gene_expression_tl/caRNA_cds_expression.tab")

nor_ratio <- c(41608064,68341489,
               #70516802,68367198,
               57358338,70205279,
               #70496756,70107156,
               67542625,52345318,
               #75806743,76023704,
               63139881,69162349)
#78116866,66711470)
nor_ratio <- mean(nor_ratio)/nor_ratio

express_p0p10_cds_nor <- caRNA_count_cds_df[,c(1,2,5,6,9,10,13,14)]
for ( i in 1:8){
  express_p0p10_cds_nor[,i] <- express_p0p10_cds_nor[,i]*nor_ratio[i]
}

express_p0p10_cds_nor  <- express_p0p10_cds_nor[which(apply(express_p0p10_cds_nor>10,1,sum)>=4),]
# scale_factor <- c(0.995575897,1.004443762,
#                   #                  1.466860346,1.637934483,
#                   1.250982735,1.213976657)
# #                  0.954841898,0.98519686)
scale_factor <- c(1,1.232340796)

m6A_level_cds_log2FC <- data.frame(trans_id =substr(rownames(express_p0p10_cds_nor),1,18))
# for (i in 1:2){
#   m6A_level_cds_log2FC[,i+1]=log2((express_p0p10_cds_nor[,2*i+3]+express_p0p10_cds_nor[,2*i+4])*scale_factor[i]/(express_p0p10_cds_nor[,2*i-1]+express_p0p10_cds_nor[,2*i]))
# }
for (i in 1:2){
  m6A_level_cds_log2FC[,i+1]=log2(((express_p0p10_cds_nor[,2*i+3]+express_p0p10_cds_nor[,2*i+4])*scale_factor[i]+0.001)/((express_p0p10_cds_nor[,2*i-1]+express_p0p10_cds_nor[,2*i]+0.001)))
}
#colnames(m6A_level_cds_log2FC) <- c("trans_id","p0_rep1","p0_rep2","p10_rep1","p10_rep2")
colnames(m6A_level_cds_log2FC) <- c("trans_id","p0_m6A","p10_m6A")
rownames(m6A_level_cds_log2FC) <- m6A_level_cds_log2FC$trans_id

m6A_level_cds_log2FC <- m6A_level_cds_log2FC[,c(2:3)]
m6A_level_cds_log2FC <- m6A_level_cds_log2FC[which(((m6A_level_cds_log2FC$p0_m6A>0) + (m6A_level_cds_log2FC$p10_m6A>0))>=1),]


## PDUI
dapars_caRNA_p0p10_change <- dapars_caRNA_p0p10[which(dapars_caRNA_p0p10$filter != "NC"),]

dapars_caRNA_p0p10_change_hm <- as.matrix(dapars_caRNA_p0p10_change[,c(2,3)])
rownames(dapars_caRNA_p0p10_change_hm) <- substr(dapars_caRNA_p0p10_change$Gene,1,18)

m6A_dapars_cds_change <- m6A_level_cds_log2FC[rownames(dapars_caRNA_p0p10_change_hm),]

m6A_dapars_cds_change <- m6A_dapars_cds_change[which(apply(apply(m6A_dapars_cds_change,2,is.finite),1,sum)==2),]
dapars_caRNA_p0p10_change_cds_hm <- dapars_caRNA_p0p10_change_hm[rownames(m6A_dapars_cds_change),]

# heatmap
library(ComplexHeatmap)

heatmap_PDUI <- Heatmap(t(scale(t(dapars_caRNA_p0p10_change_cds_hm))),show_row_dend = FALSE,  
                        #row_km=3,
                        #column_km=2, 
                        #cluster_rows = FALSE,
                        show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                        heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")),
                        col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#B210FF", colorRampPalette(c("#B210FF", "white", "#EECE13"))(6), "#EECE13")))
#col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#81FFEF", colorRampPalette(c("#81FFEF", "white", "#F067B4"))(6), "#F067B4")))
#col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.5), 1), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(3), "#FF0000")))
heatmap_m6A <- Heatmap(t(scale(t(m6A_dapars_cds_change))),show_row_dend = FALSE,  
                       #row_km=3,
                       #column_km=2, 
                       #cluster_rows = FALSE,
                       show_row_names = FALSE, border=TRUE, #right_annotation = row_anno,
                       heatmap_legend_param=list(title = "",direction = "horizontal", legend_width = unit(5, "cm")),
                       col = circlize::colorRamp2(c(-1, seq(-0.5,0.5, by=0.2), 1), c("#B210FF", colorRampPalette(c("#B210FF", "white", "#EECE13"))(6), "#EECE13")))
#col = circlize::colorRamp2(c(-3, seq(-2.5,2.5, by=0.2), 3), c("#2166AC", colorRampPalette(c("#2166AC", "white", "#FF0000"))(11), "#FF0000")))
ht_list =heatmap_m6A+heatmap_PDUI
draw(ht_list,heatmap_legend_side = "top")

ht_list =heatmap_PDUI+heatmap_m6A
draw(ht_list,heatmap_legend_side = "top")



















#### peak

#generate features
m6A_peak_p0p10_saf <- read.table("./10_bed_merge/02_merge_peaks/p0p10.merge.saf")
colnames(m6A_peak_p0p10_saf) <- c("GeneID","Chr","Start","End","Strand")
m6A_peak_p0p10_bed <- cbind(m6A_peak_p0p10_saf,".")
m6A_peak_p0p10_bed <- m6A_peak_p0p10_bed[,c(2,3,4,1,6,5)]
write.table(m6A_peak_p0p10_bed,file ="./10_bed_merge/02_merge_peaks/p0p10.merge.bed",quote=F,sep="\t",row.names = F,col.names = F)

p0p10_peak_count <- featureCounts(files=paste0(path,c("p0_input_rep1","p0_input_rep2",
                                                      "p10_input_rep1","p10_input_rep2",
                                                      "p0_ip_rep1","p0_ip_rep2",
                                                      "p10_ip_rep1","p10_ip_rep2"),".bam"),
                                  annot.ext = m6A_peak_p0p10_saf,
                                  isGTFAnnotationFile = FALSE, 
                                  isPairedEnd=TRUE,
                                  requireBothEndsMapped=TRUE,
                                  allowMultiOverlap = TRUE,
                                  fracOverlap = 0.2,
                                  countChimericFragments=FALSE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific = 2,
                                  checkFragLength = TRUE,
                                  minFragLength = 40, 
                                  maxFragLength = 2000, 
                                  verbose = FALSE, nthreads=40)
p0p10_peak_count_df <- p0p10_peak_count$counts
write.table(p0p10_peak_count_df,"./10_bed_merge/02_merge_peaks/p0p10_peak_counts.tab")


peak_annotation <- read.delim("~/LinLong/10_bed_merge/raw/annotation_output/Lysate_Result_merge_peaks_annotation_homer_rearranged.xls")

# get annotation
#anno <- rtracklayer::import("~/reference/annotation/mm19/gencode.vM28.annotation.gtf")
#anno <- rtracklayer::as.data.frame(anno)


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



