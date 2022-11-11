library(data.table)

counts_repeats <- fread(file="/disk/user_08/Data/TC1-planB/05_bam_change_index_hisat2/01_bam_sorted/11_bam_merge/06_repeats/repeats_family_separate/repeats_ERVL_counts.tab",nThread = 40)
counts_repeats <- counts_repeats[,-c(2,3,4,5)]
rownames(counts_repeats) <- counts_repeats$Geneid
colnames(counts_repeats) <- c("Geneid","Length","p0_input_rep1","p0_input_rep2","p0_ip_rep1","p0_ip_rep2","p5_input_rep1","p5_input_rep2","p5_ip_rep1","p5_ip_rep2","p10_input_rep1","p10_input_rep2","p10_ip_rep1","p10_ip_rep2","rp2_input_rep1","rp2_input_rep2","rp2_ip_rep1","rp2_ip_rep2")


counts <- counts_repeats[,c(3,4,7,8,11,12,15,16)]
counts <- as.data.frame(counts)
rownames(counts) <- counts_repeats$Geneid


library_ratio <- c(41608064,68341489,70516802,68367198,57358338,70205279,70496756,70107156)
library_ratio <- mean(library_ratio)/library_ratio

library(edgeR)

counts_library_nor <- counts
for (i in 1:8){
  counts_library_nor[,i] <- counts[,i]*library_ratio[i]
}

counts_library_nor <- counts_library_nor[which(apply(counts_library_nor,1,sum)>10000),]

library(reshape2)
library(ggplot2)
counts_scale <- as.data.frame(t(scale(t(counts_library_nor))))
#counts_scale$family_name <- rownames(counts_scale)

counts_scale <- na.omit(counts_scale)

#counts_m <- melt(counts_scale, id.vars=c("family_name"))

library(ComplexHeatmap)
meta_data <- data.frame(rep = rep(c("rep1","rep2"),4),treatment = rep(c("p0","p5","p10","rp2"),each=2),row.names = colnames(counts_library_nor),stringsAsFactors = TRUE)
cell_anno <- HeatmapAnnotation(treat = meta_data[,"treatment"])


Heatmap(matrix = counts_scale,show_column_names =TRUE,show_row_names=TRUE,
        cluster_rows = TRUE,cluster_columns = FALSE,
        column_labels = c("p0_rep1","p0_rep2","p5_rep1","p5_rep2","p10_rep1","p10_rep2","rp2_rep1","rp2_rep2"),
        width = unit(2,'inches'))

counts_cpm <- cpm(counts)
counts_cpm <- counts_cpm[which(apply(counts_cpm,1,sum)>100),]


counts_cpm_scale <- as.data.frame(t(scale(t(counts_cpm))))

counts_cpm_scale <- na.omit(counts_cpm_scale)

Heatmap(matrix = counts_cpm_scale,show_column_names =TRUE,show_row_names=TRUE,
        cluster_rows = TRUE,cluster_columns = FALSE,
        column_labels = c("p0_rep1","p0_rep2","p5_rep1","p5_rep2","p10_rep1","p10_rep2","rp2_rep1","rp2_rep2"),
        width = unit(2,'inches'))