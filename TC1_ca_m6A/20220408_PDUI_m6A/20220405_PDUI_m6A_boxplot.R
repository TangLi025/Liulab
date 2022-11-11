UTR3_bed <- import("/disk/user_08/Data/TC1-planB/07_APA/03_m6A_input/00_file/p0p10_3UTR_altered.bed")
UTR3_bed <- as.data.frame(UTR3_bed)
UTR3_saf <- UTR3_bed[,c(6,1,2,3,5)]
colnames(UTR3_saf) <- c("GeneID","Chr","Start","End","Strand")

caRNA_count_3UTR_alter <- featureCounts(files=paste0(path,sample_name,".bam"),
                                  annot.ext = UTR3_saf,
                                  isGTFAnnotationFile = FALSE, 
                                  isPairedEnd=TRUE,
                                  requireBothEndsMapped=TRUE,
                                  allowMultiOverlap = TRUE,
                                  fracOverlap = 1,
                                  fracOverlapFeature = 1,
                                  countChimericFragments=FALSE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific = 2,checkFragLength = TRUE,
                                  minFragLength = 40, maxFragLength = 2000, 
                                  verbose = FALSE, nthreads=40)

caRNA_count_3UTR_df <- caRNA_count_3UTR_alter$counts
colnames(caRNA_count_3UTR_df) <- sample_name

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
# for (i in 1:2){
#   m6A_level_3UTR_log2FC[,i+1]=log2((express_p0p10_3UTR_nor[,2*i+3]+express_p0p10_3UTR_nor[,2*i+4])*scale_factor[i]/(express_p0p10_3UTR_nor[,2*i-1]+express_p0p10_3UTR_nor[,2*i]))
# }
for (i in 1:2){
  m6A_level_3UTR_log2FC[,i+1]=log2(((express_p0p10_3UTR_nor[,2*i+3]+express_p0p10_3UTR_nor[,2*i+4])+0.001)/((express_p0p10_3UTR_nor[,2*i-1]+express_p0p10_3UTR_nor[,2*i]+0.001)))
}
colnames(m6A_level_3UTR_log2FC) <- c("trans_id","p0_m6A","p10_m6A")
rownames(m6A_level_3UTR_log2FC) <- m6A_level_3UTR_log2FC$trans_id

m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[,c(2:3)]
#m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[which(((m6A_level_3UTR_log2FC$p0_m6A>0) + (m6A_level_3UTR_log2FC$p10_m6A>0))>=1),]


# ## PDUI
# 
# 
# dapars_caRNA_p0p10_filter_hm <- dapars_caRNA_p0p10
# dapars_caRNA_p0p10_filter_hm$Gene <- substr(dapars_caRNA_p0p10_filter_hm$Gene,1,18)
# 
# dapars_caRNA_p0p10_filter_hm %>% mutate(m6A_p10_3UTR=m6A_level_3UTR_log2FC[dapars_caRNA_p0p10_filter_hm$Gene,"p10_m6A"],
#                                         m6A_p0_3UTR=m6A_level_3UTR_log2FC[dapars_caRNA_p0p10_filter_hm$Gene,"p0_m6A"]) -> dapars_caRNA_p0p10_filter_hm
# 
# dapars_caRNA_p0p10_filter_hm %>% mutate(log2_m6A_p10_p0_3UTR=(m6A_p10_3UTR+0.00001)/(m6A_p0_3UTR+0.00001)) -> dapars_caRNA_p0p10_filter_hm
#                            
# PDUI_comparisons <-list(c("UP","NC"),c("UP","DOWN"),c("NC","DOWN"))
# 
# ggplot(dapars_caRNA_p0p10_filter_hm, aes(x=filter, y=log2_m6A_p10_p0_3UTR, fill=filter)) + 
#   stat_boxplot(geom = "errorbar",width=0.15) + #添加虚线
#   geom_boxplot(outlier.shape = NA)+ 
#   stat_compare_means(comparisons = PDUI_comparisons,
#                      label = "p.signif",
#                      method = "t.test",
#                      label.y = c(3.1,3.2,3.7))


m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[which(((m6A_level_3UTR_log2FC$p0_m6A>0) + (m6A_level_3UTR_log2FC$p10_m6A>0))>=1),]


## PDUI
dapars_caRNA_p0p10 <- read.table("/disk/user_09/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p0p10.txt",
                                 header = T,row.names = 1)
dapars_caRNA_p0p10_change <- dapars_caRNA_p0p10[which(dapars_caRNA_p0p10$filter != "NC"),]

dapars_caRNA_p0p10_change_hm <- as.matrix(dapars_caRNA_p0p10_change[,c(1,2)])
rownames(dapars_caRNA_p0p10_change_hm) <- substr(rownames(dapars_caRNA_p0p10_change),1,18)

m6A_dapars_3UTR_change <- m6A_level_3UTR_log2FC[rownames(dapars_caRNA_p0p10_change_hm),]

m6A_dapars_3UTR_change <- m6A_dapars_3UTR_change[which(apply(apply(m6A_dapars_3UTR_change,2,is.finite),1,sum)==2),]
dapars_caRNA_p0p10_change_3UTR_hm <- dapars_caRNA_p0p10_change_hm[rownames(m6A_dapars_3UTR_change),]



m6A_dapars_3UTR_filter <- m6A_level_3UTR_log2FC[rownames(dapars_caRNA_p0p10_change_hm),]

m6A_dapars_3UTR_change <- m6A_dapars_3UTR_change[which(apply(apply(m6A_dapars_3UTR_change,2,is.finite),1,sum)==2),]
dapars_caRNA_p0p10_change_3UTR_hm <- as.data.frame(dapars_caRNA_p0p10_change_hm[rownames(m6A_dapars_3UTR_change),])

library(tidyverse)
compare_PDUI_m6A <- tibble(trans=rownames(dapars_caRNA_p0p10_change_3UTR_hm),
                           PDUI_p10_p0=ifelse((dapars_caRNA_p0p10_change_3UTR_hm$caRNA_p10+0.000001)/(dapars_caRNA_p0p10_change_3UTR_hm$caRNA_p0+0.000001)>1,"longer","shorter"),
                           m6A_p10_p0 = m6A_dapars_3UTR_change$p10_m6A-m6A_dapars_3UTR_change$p0_m6A)

ggplot(compare_PDUI_m6A, aes(x=PDUI_p10_p0, y=m6A_p10_p0, fill=PDUI_p10_p0)) + 
  stat_boxplot(geom = "errorbar",width=0.15) + #添加虚线
  geom_boxplot()+ 
  stat_compare_means(comparisons = list(c("longer","shorter")),
                     label = "p.signif",
                     method = "t.test")

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

UTR3_bed <- import("/disk/user_09/reference/annotation/mm19/dapars/mm39_gencode_dapars_longest_protein.bed")
UTR3_bed <- as.data.frame(UTR3_bed)
UTR3_saf <- UTR3_bed[,c(6,1,2,3,5)]
colnames(UTR3_saf) <- c("GeneID","Chr","Start","End","Strand")

path <- "/disk/user_08/Data/TC1-planB/05_bam_change_index_hisat2/01_bam_sorted/11_bam_merge/"
caRNA_count_3UTR <- featureCounts(files=paste0(path,sample_name,".bam"),
                                        annot.ext = UTR3_saf,
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
  m6A_level_3UTR_log2FC[,i+1]=log2(((express_p0p10_3UTR_nor[,2*i+3]+express_p0p10_3UTR_nor[,2*i+4])+0.001)/((express_p0p10_3UTR_nor[,2*i-1]+express_p0p10_3UTR_nor[,2*i]+0.001)))
}
colnames(m6A_level_3UTR_log2FC) <- c("trans_id","p0_m6A","p10_m6A")
rownames(m6A_level_3UTR_log2FC) <- m6A_level_3UTR_log2FC$trans_id

m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[,c(2:3)]
#m6A_level_3UTR_log2FC <- m6A_level_3UTR_log2FC[which(((m6A_level_3UTR_log2FC$p0_m6A>0) + (m6A_level_3UTR_log2FC$p10_m6A>0))>=1),]


# ## PDUI
# 
# 
# dapars_caRNA_p0p10_filter_hm <- dapars_caRNA_p0p10
# dapars_caRNA_p0p10_filter_hm$Gene <- substr(dapars_caRNA_p0p10_filter_hm$Gene,1,18)
# 
# dapars_caRNA_p0p10_filter_hm %>% mutate(m6A_p10_3UTR=m6A_level_3UTR_log2FC[dapars_caRNA_p0p10_filter_hm$Gene,"p10_m6A"],
#                                         m6A_p0_3UTR=m6A_level_3UTR_log2FC[dapars_caRNA_p0p10_filter_hm$Gene,"p0_m6A"]) -> dapars_caRNA_p0p10_filter_hm
# 
# dapars_caRNA_p0p10_filter_hm %>% mutate(log2_m6A_p10_p0_3UTR=(m6A_p10_3UTR+0.00001)/(m6A_p0_3UTR+0.00001)) -> dapars_caRNA_p0p10_filter_hm
#                            
# PDUI_comparisons <-list(c("UP","NC"),c("UP","DOWN"),c("NC","DOWN"))
# 
# ggplot(dapars_caRNA_p0p10_filter_hm, aes(x=filter, y=log2_m6A_p10_p0_3UTR, fill=filter)) + 
#   stat_boxplot(geom = "errorbar",width=0.15) + #添加虚线
#   geom_boxplot(outlier.shape = NA)+ 
#   stat_compare_means(comparisons = PDUI_comparisons,
#                      label = "p.signif",
#                      method = "t.test",
#                      label.y = c(3.1,3.2,3.7))

## PDUI
dapars_caRNA_p0p10_change <- dapars_caRNA_p0p10[which(dapars_caRNA_p0p10$filter != "NC"),]

dapars_caRNA_p0p10_change_hm <- as.matrix(dapars_caRNA_p0p10_change[,c(1,2)])
rownames(dapars_caRNA_p0p10_change_hm) <- substr(rownames(dapars_caRNA_p0p10_change),1,18)

m6A_dapars_3UTR_change <- m6A_level_3UTR_log2FC[rownames(dapars_caRNA_p0p10_change_hm),]

m6A_dapars_3UTR_change <- m6A_dapars_3UTR_change[which(apply(apply(m6A_dapars_3UTR_change,2,is.finite),1,sum)==2),]
dapars_caRNA_p0p10_change_3UTR_hm <- dapars_caRNA_p0p10_change_hm[rownames(m6A_dapars_3UTR_change),]



m6A_dapars_3UTR_filter <- m6A_level_3UTR_log2FC[rownames(dapars_caRNA_p0p10_change_hm),]

m6A_dapars_3UTR_change <- m6A_dapars_3UTR_change[which(apply(apply(m6A_dapars_3UTR_change,2,is.finite),1,sum)==2),]
dapars_caRNA_p0p10_change_3UTR_hm <- as.data.frame(dapars_caRNA_p0p10_change_hm[rownames(m6A_dapars_3UTR_change),])

compare_PDUI_m6A <- tibble(trans=rownames(dapars_caRNA_p0p10_change_3UTR_hm),
                           PDUI_p10_p0=ifelse((dapars_caRNA_p0p10_change_3UTR_hm$caRNA_p10+0.000001)/(dapars_caRNA_p0p10_change_3UTR_hm$caRNA_p0+0.000001)>1,"lengthen","shorten"),
                           m6A_p10_p0 = m6A_dapars_3UTR_change$p10_m6A-m6A_dapars_3UTR_change$p0_m6A)

ggplot(compare_PDUI_m6A, aes(x=PDUI_p10_p0, y=m6A_p10_p0, fill=PDUI_p10_p0)) + 
  stat_boxplot(geom = "errorbar",width=0.15) + #添加虚线
  geom_boxplot()+ 
  stat_compare_means(comparisons = list(c("lengthen","shorten")),
                     label = "p.signif",
                     method = "t.test")


library(ggplot2)
library(ggpubr)

  
