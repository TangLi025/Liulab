library(tidyverse)

PDUI_p0p10 <- read.table("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p0p10.txt",
                         header = T)
PDUI_p0p10 <- as_tibble(PDUI_p0p10)

p0p10_short <- PDUI_p0p10 %>% filter(filter=="UP")

APA_loci_p0p10_short <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p0p10_short$Gene))

p0p10_long <- PDUI_p0p10 %>% filter(filter=="DOWN")

APA_loci_p0p10_long <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p0p10_long$Gene))

p0p10_nc <- PDUI_p0p10 %>% filter(filter=="NC")

APA_loci_p0p10_nc <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p0p10_nc$Gene))

write.table(APA_loci_p0p10_short,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_region_p0p10_short.bed",quote = F,sep="\t",row.names = F,col.names = F)
write.table(APA_loci_p0p10_long,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_region_p0p10_long.bed",quote = F,sep="\t",row.names = F,col.names = F)
write.table(APA_loci_p0p10_nc,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_region_p0p10_nc.bed",quote = F,sep="\t",row.names = F,col.names = F)

## p0p5

PDUI_p0p5 <- read.table("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p0p5.txt",
                         header = T)
PDUI_p0p5 <- as_tibble(PDUI_p0p5)

p0p5_short <- PDUI_p0p5 %>% filter(filter=="UP")

APA_loci_p0p5_short <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p0p5_short$Gene))

p0p5_long <- PDUI_p0p5 %>% filter(filter=="DOWN")

APA_loci_p0p5_long <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p0p5_long$Gene))

p0p5_nc <- PDUI_p0p5 %>% filter(filter=="NC")

APA_loci_p0p5_nc <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p0p5_nc$Gene))

write.table(APA_loci_p0p5_short,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p0p5_short.bed",quote = F,sep="\t",row.names = F,col.names = F)
write.table(APA_loci_p0p5_long,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p0p5_long.bed",quote = F,sep="\t",row.names = F,col.names = F)
write.table(APA_loci_p0p5_nc,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p0p5_nc.bed",quote = F,sep="\t",row.names = F,col.names = F)

## p10rp2
PDUI_p10rp2 <- read.table("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p10rp2.txt",
                         header = T)
PDUI_p10rp2 <- as_tibble(PDUI_p10rp2)

p10rp2_short <- PDUI_p10rp2 %>% filter(filter=="UP")

APA_loci_p10rp2_short <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p10rp2_short$Gene))

p10rp2_long <- PDUI_p10rp2 %>% filter(filter=="DOWN")

APA_loci_p10rp2_long <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p10rp2_long$Gene))

p10rp2_nc <- PDUI_p10rp2 %>% filter(filter=="NC")

APA_loci_p10rp2_nc <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p10rp2_nc$Gene))

write.table(APA_loci_p10rp2_short,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p10rp2_short.bed",quote = F,sep="\t",row.names = F,col.names = F)
write.table(APA_loci_p10rp2_long,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p10rp2_long.bed",quote = F,sep="\t",row.names = F,col.names = F)
write.table(APA_loci_p10rp2_nc,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p10rp2_nc.bed",quote = F,sep="\t",row.names = F,col.names = F)

## p5rp2
PDUI_p5rp2 <- read.table("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p5rp2.txt",
                          header = T)
PDUI_p5rp2 <- as_tibble(PDUI_p5rp2)

p5rp2_short <- PDUI_p5rp2 %>% filter(filter=="UP")

APA_loci_p5rp2_short <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p5rp2_short$Gene))

p5rp2_long <- PDUI_p5rp2 %>% filter(filter=="DOWN")

APA_loci_p5rp2_long <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p5rp2_long$Gene))

p5rp2_nc <- PDUI_p5rp2 %>% filter(filter=="NC")

APA_loci_p5rp2_nc <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p5rp2_nc$Gene))

write.table(APA_loci_p5rp2_short,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p5rp2_short.bed",quote = F,sep="\t",row.names = F,col.names = F)
write.table(APA_loci_p5rp2_long,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p5rp2_long.bed",quote = F,sep="\t",row.names = F,col.names = F)
write.table(APA_loci_p5rp2_nc,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p5rp2_nc.bed",quote = F,sep="\t",row.names = F,col.names = F)

## p0rp2
PDUI_p0rp2 <- read.table("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p0rp2.txt",
                         header = T)
PDUI_p0rp2 <- as_tibble(PDUI_p0rp2)

p0rp2_short <- PDUI_p0rp2 %>% filter(filter=="UP")

APA_loci_p0rp2_short <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p0rp2_short$Gene))

p0rp2_long <- PDUI_p0rp2 %>% filter(filter=="DOWN")

APA_loci_p0rp2_long <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p0rp2_long$Gene))

p0rp2_nc <- PDUI_p0rp2 %>% filter(filter=="NC")

APA_loci_p0rp2_nc <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p0rp2_nc$Gene))

write.table(APA_loci_p0rp2_short,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p0rp2_short.bed",quote = F,sep="\t",row.names = F,col.names = F)
write.table(APA_loci_p0rp2_long,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p0rp2_long.bed",quote = F,sep="\t",row.names = F,col.names = F)
write.table(APA_loci_p0rp2_nc,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p0rp2_nc.bed",quote = F,sep="\t",row.names = F,col.names = F)

## p5p10
PDUI_p5p10 <- read.table("/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/group_table_0.2_0.58/DaPars2_caRNA_p5p10.txt",
                         header = T)
PDUI_p5p10 <- as_tibble(PDUI_p5p10)

p5p10_short <- PDUI_p5p10 %>% filter(filter=="UP")

APA_loci_p5p10_short <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p5p10_short$Gene))

p5p10_long <- PDUI_p5p10 %>% filter(filter=="DOWN")

APA_loci_p5p10_long <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p5p10_long$Gene))

p5p10_nc <- PDUI_p5p10 %>% filter(filter=="NC")

APA_loci_p5p10_nc <- APA_loci_bed %>% filter(name %in% gsub("\\..*","",p5p10_nc$Gene))

write.table(APA_loci_p5p10_short,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p5p10_short.bed",quote = F,sep="\t",row.names = F,col.names = F)
write.table(APA_loci_p5p10_long,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p5p10_long.bed",quote = F,sep="\t",row.names = F,col.names = F)
write.table(APA_loci_p5p10_nc,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/APA_loci_p5p10_nc.bed",quote = F,sep="\t",row.names = F,col.names = F)
