library(data.table)

fc_gene <- fread("/disk/user_08/Data/TC1-planB/10_mRNA/05_bam_hisat2/01_bam_sorted/11_bam_merge/mRNA_gene_expression_matrix.txt")

x_gene <- fc_gene[,c(8,9,12,13,16,17,20,21)]

counts <- x_gene
counts <- as.data.frame(counts)
rownames(counts) <- fc_gene$Geneid


# library_ratio <- c(85439160,141246126,141609654,137300068,115590222,141466700,142031524,141270996)
# library_ratio <- mean(library_ratio)/library_ratio
# 
# library(edgeR)
# 
# counts_library_nor <- counts
# for (i in 1:8){
#   counts_library_nor[,i] <- counts[,i]*library_ratio[i]
# }
# 
# counts_library_nor <- counts_library_nor[which(apply(counts_library_nor,1,sum)>10000),]
library(tidyverse)
library(edgeR)
cpm_counts <- cpm(counts)

# cpm_counts <- cpm_counts[which(apply(cpm_counts,1,sum)>4),]
library(reshape2)
library(ggplot2)

cpm_counts_toti <- cpm_counts[which(rownames(cpm_counts) %in% c("Cdkn1a","Btg1","Ctsb","Plk2","Ddit4l","Zscan4c","Zscan4f","Gm8300","Prnp","Zfp352","Bmpr2","Egfr","Nrp1","Mafk","Mmp19","Trp53inp2","Snai1","Btg2","Gm5662","Zscan4d","Gss","Rab43","Rab11b","Tgfbr3","Sp110","Ddr2")),]


counts_scale_toti <- as.data.frame(t(scale(t(cpm_counts_toti))))
#counts_scale$family_name <- rownames(counts_scale)

counts_scale_toti <- na.omit(counts_scale_toti)

#cpm_counts_plur <- cpm_counts[which(rownames(cpm_counts) %in% c("Oct4","Klf4","Nanog","Sox2","Sox2ot","Pou5f1")),]
#cpm_counts_plur <- cpm_counts[which(rownames(cpm_counts) %in% c("Zfp42","Tet1","Tcf15","Nr0b1","H2afz","Dppa5a")),]

#cpm_counts_plur <- cpm_counts[which(rownames(cpm_counts) %in% c("Idh2","Utf1","Tdgf1","Sox2","Prdm14","Upp1","Bhlhe40","Klf15","Pou5f1","Pou2f3","Grhl3")),]

#cpm_counts_plur <- cpm_counts[which(rownames(cpm_counts) %in% c("Oct4","Klf4","Nanog","Sox2","Sox2ot","Pou5f1","Zfp42","Tet1","Tcf15","Nr0b1","H2afz","Dppa5a","Idh2","Utf1","Tdgf1","Sox2","Prdm14","Upp1","Bhlhe40","Klf15","Pou5f1","Pou2f3","Grhl3")),]

cpm_counts_plur <- cpm_counts[which(rownames(cpm_counts) %in% c("Pou5f1","Tcf15","H2afz","Dppa5a","Utf1","Tdgf1","Sox2","Upp1","Bhlhe40","Klf15","Pou2f3","Grhl3")),]

counts_scale_plur <- as.data.frame(t(scale(t(cpm_counts_plur))))
#counts_scale$family_name <- rownames(counts_scale)

counts_scale_plur <- na.omit(counts_scale_plur)

#counts_m <- melt(counts_scale, id.vars=c("family_name"))

library(ComplexHeatmap)
meta_data <- data.frame(rep = rep(c("rep1","rep2"),4),treatment = rep(c("p0","p5","p10","rp2"),each=2),row.names = colnames(cpm_counts),stringsAsFactors = TRUE)
cell_anno <- HeatmapAnnotation(treat = meta_data[,"treatment"])

pdf("./04_gene_expression/01_pic/toti_pluri_genes_mRNA_cpm_dupaper.pdf")
Heatmap(matrix = counts_scale_toti,show_column_names =TRUE,show_row_names=TRUE,
        cluster_rows = TRUE,cluster_columns = FALSE,
        column_labels = c("p0_rep1","p0_rep2","p5_rep1","p5_rep2","p10_rep1","p10_rep2","rp2_rep1","rp2_rep2"),
        width = unit(2,'inches'))
Heatmap(matrix = counts_scale_plur,show_column_names =TRUE,show_row_names=TRUE,
        cluster_rows = TRUE,cluster_columns = FALSE,
        column_labels = c("p0_rep1","p0_rep2","p5_rep1","p5_rep2","p10_rep1","p10_rep2","rp2_rep1","rp2_rep2"),
        width = unit(2,'inches'))

dev.off()

#########################

# m6A genes

cpm_counts_m6A <- cpm_counts[which(rownames(cpm_counts) %in% c("Mettl3","Mettl14",
      "Mettl16","Wtap","Virma","Cbll1","Zc3h13","Rbm15","Rbm15b","Fto","Alkbh5",
      "Ythdf1","Ythdf2","Ythdf3","Ythdc1","Ythdc2","Hnrnpa2B1","Eif3A","Igf2bp1",
      "Igf2bp2","Igf2bp3","Fmr1","Hnrnpc","Rbmx","Elavl1","G3bp1","G3bp2","Nudt21","Dicer1")),]

counts_scale_m6A <- as.data.frame(t(scale(t(cpm_counts_m6A))))
#counts_scale$family_name <- rownames(counts_scale)

counts_scale_m6A <- na.omit(counts_scale_m6A)

Heatmap(matrix = counts_scale_m6A,show_column_names =TRUE,show_row_names=TRUE,
        cluster_rows = TRUE,cluster_columns = FALSE,
        column_labels = c("p0_rep1","p0_rep2","p5_rep1","p5_rep2","p10_rep1","p10_rep2","rp2_rep1","rp2_rep2"),
        width = unit(2,'inches'))

##### copper induced death

cpm_counts_copper <- cpm_counts[which(rownames(cpm_counts) %in% c("Fdx1","Lias","Lipt1","Dld","Dlat",
                                                                  "Pdha1","Pdhb",
                                                                  "Mtf1","Gls","Cdkn2a")),]

counts_scale_copper <- as.data.frame(t(scale(t(cpm_counts_copper))))
#counts_scale$family_name <- rownames(counts_scale)

counts_scale_copper <- na.omit(counts_scale_copper)

Heatmap(matrix = counts_scale_copper,show_column_names =TRUE,show_row_names=TRUE,
        cluster_rows = TRUE,cluster_columns = FALSE,
        column_labels = c("p0_rep1","p0_rep2","p5_rep1","p5_rep2","p10_rep1","p10_rep2","rp2_rep1","rp2_rep2"),
        width = unit(2,'inches'))
