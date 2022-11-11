library(data.table)

#### p0
p0_peak_bed <- fread("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/10_bed_merge/00_common_peaks/p0_rep1_rep2_common_peaks.bed")

p0_peak_center <- p0_peak_bed

p0_peak_center$V2 <- round((p0_peak_bed$V2+p0_peak_bed$V3)/2)
p0_peak_center$V3 <- p0_peak_center$V2+1

write.table(p0_peak_center,"/disk/user_09/Data/03_TC1_caRNA/10_bed_merge/01_peak_center/p0_peak_center.bed",
            quote = F,sep = "\t",row.names = F,col.names = F)

#### p5
p5_peak_bed <- fread("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/10_bed_merge/00_common_peaks/p5_rep1_rep2_common_peaks.bed")

p5_peak_center <- p5_peak_bed

p5_peak_center$V2 <- round((p5_peak_bed$V2+p5_peak_bed$V3)/2)
p5_peak_center$V3 <- p5_peak_center$V2+1

write.table(p5_peak_center,"/disk/user_09/Data/03_TC1_caRNA/10_bed_merge/01_peak_center/p5_peak_center.bed",
            quote = F,sep = "\t",row.names = F,col.names = F)

#### p10
p10_peak_bed <- fread("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/10_bed_merge/00_common_peaks/p10_rep1_rep2_common_peaks.bed")

p10_peak_center <- p10_peak_bed

p10_peak_center$V2 <- round((p10_peak_bed$V2+p10_peak_bed$V3)/2)
p10_peak_center$V3 <- p10_peak_center$V2+1

write.table(p10_peak_center,"/disk/user_09/Data/03_TC1_caRNA/10_bed_merge/01_peak_center/p10_peak_center.bed",
            quote = F,sep = "\t",row.names = F,col.names = F)


#### rp2
rp2_peak_bed <- fread("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/10_bed_merge/00_common_peaks/rp2_rep1_rep2_common_peaks.bed")

rp2_peak_center <- rp2_peak_bed

rp2_peak_center$V2 <- round((rp2_peak_bed$V2+rp2_peak_bed$V3)/2)
rp2_peak_center$V3 <- rp2_peak_center$V2+1

write.table(rp2_peak_center,"/disk/user_09/Data/03_TC1_caRNA/10_bed_merge/01_peak_center/rp2_peak_center.bed",
            quote = F,sep = "\t",row.names = F,col.names = F)
