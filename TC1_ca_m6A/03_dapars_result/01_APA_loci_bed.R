dapars_all <- read.table("/disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_caRNA_new/DaPars2_caRNA_result.txt",
                         header = T)

dapars_all$chr <- gsub("\\:.*","",dapars_all$Loci)

dapars_all$strand <- substr(dapars_all$Gene,nchar(dapars_all$Gene),nchar(dapars_all$Gene))

for (i in 1:nrow(dapars_all)){
  dapars_all$start[i] <- as.numeric(strsplit(strsplit(dapars_all$Loci[i], split = ":", fixed = TRUE)[[1]][2],split = "-",fixed = TRUE)[[1]][1])
}

for (i in 1:nrow(dapars_all)){
  dapars_all$end[i] <- as.numeric(strsplit(strsplit(dapars_all$Loci[i], split = ":", fixed = TRUE)[[1]][2],split = "-",fixed = TRUE)[[1]][2])
}


library(tidyverse)
APA_loci_bed <- tibble(chr=dapars_all$chr,
                       start=dapars_all$Predicted_Proximal_APA,
                       end=dapars_all$Predicted_Proximal_APA+1,
                       name=gsub("\\..*","",dapars_all$Gene),
                       score=".",
                       strand=dapars_all$strand)

write.table(APA_loci_bed,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_caRNA_new/APA_loci.bed",quote = F,sep="\t",row.names = F,col.names = F)

APA_region_bed <- tibble(chr=dapars_all$chr,
                         start=ifelse(dapars_all$strand=="+",
                                      dapars_all$start,
                                      dapars_all$Predicted_Proximal_APA),
                       end=ifelse(dapars_all$strand=="+",
                                  dapars_all$Predicted_Proximal_APA,
                                  dapars_all$end),
                       name=gsub("\\..*","",dapars_all$Gene),
                       score=".",
                       strand=dapars_all$strand)

write.table(APA_region_bed,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_caRNA_new/pAPA_region.bed",quote = F,sep="\t",row.names = F,col.names = F)

dAPA_region_bed <- tibble(chr=dapars_all$chr,
                          start=ifelse(dapars_all$strand=="+",
                                       dapars_all$Predicted_Proximal_APA,
                                       dapars_all$start),
                         end=ifelse(dapars_all$strand=="+",
                                    dapars_all$end,
                                    dapars_all$Predicted_Proximal_APA),
                         name=gsub("\\..*","",dapars_all$Gene),
                         score=".",
                         strand=dapars_all$strand)

write.table(dAPA_region_bed,file="/disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_caRNA_new/dAPA_region.bed",quote = F,sep="\t",row.names = F,col.names = F)
