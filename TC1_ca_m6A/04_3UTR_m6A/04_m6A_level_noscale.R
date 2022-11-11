library(rtracklayer)
library(Rsubread)
library(data.table)
library(edgeR)
library(ggrepel)

library(ggthemes)
theme_set(theme_base(base_size = 12))

library_size <- c(
  41608064,68341489, 
  #70516802,68367198, 
  57358338,70205279, 
  #70496756,70107156, 
  67542625,52345318, 
  #75806743,76023704,
63139881,69162349) 
#78116866,66711470)

##################### p0p10 ####################

path="/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/11_bam_merge/"
sample_name <- c("p0_input_rep1","p0_input_rep2",
                 "p10_input_rep1","p10_input_rep2",
                 "p0_ip_rep1","p0_ip_rep2",
                 "p10_ip_rep1","p10_ip_rep2")
bam=paste0(path,sample_name,".bam")
saf="/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/10_bed_merge/02_merge_peaks/p0p10.merge.saf"

p0p10_peak <- featureCounts(files = bam, 
                           annot.ext = saf,
                           isGTFAnnotationFile = F,
                           minMQS = 20, 
                           strandSpecific = 2,
                           countMultiMappingReads = FALSE,
                           isPairedEnd = TRUE,
                           fracOverlap = 0.2,
                           maxFragLength = 2000,                         
                           nthreads = 40)

p0p10_peak_count <- p0p10_peak$counts
colnames(p0p10_peak_count) <- sample_name
colSums(p0p10_peak_count)

p0p10_peak_nor <- as.data.frame(p0p10_peak_count)
for (i in 1:8){
  p0p10_peak_nor[,i] <- p0p10_peak_count[,i]*mean(library_size)/library_size[i]
}

colSums(p0p10_peak_nor)


p0p10_peak_nor$p0_m6A <- log2((p0p10_peak_nor$p0_ip_rep1+p0p10_peak_nor$p0_ip_rep2+1)/(p0p10_peak_nor$p0_input_rep1+p0p10_peak_nor$p0_input_rep2+1))
p0p10_peak_nor$p10_m6A <- log2((p0p10_peak_nor$p10_ip_rep1+p0p10_peak_nor$p10_ip_rep2+1)/(p0p10_peak_nor$p10_input_rep1+p0p10_peak_nor$p10_input_rep2+1))

table(((p0p10_peak_nor$p0_m6A > 0) + (p0p10_peak_nor$p10_m6A > 0))>0)
p0p10_peak_nor <- p0p10_peak_nor[((p0p10_peak_nor$p0_m6A > 0) + (p0p10_peak_nor$p10_m6A > 0))>0,]

summary(p0p10_peak_nor$p0_m6A)
summary(p0p10_peak_nor$p10_m6A)

p0p10_peak_nor$lfc=p0p10_peak_nor$p10_m6A-p0p10_peak_nor$p0_m6A

p0p10_peak_nor$group1=ifelse(abs(p0p10_peak_nor$lfc) > 1,
                            ifelse(p0p10_peak_nor$lfc > 1,'Hyper','Hypo'),'Not')
table(p0p10_peak_nor$group1)
dir.create("/disk/user_09/Data/03_TC1_caRNA/13_m6A_level/01_noscale/")
Hyper=nrow(p0p10_peak_nor[p0p10_peak_nor$group1 =='Hyper',]) 
Hypo=nrow(p0p10_peak_nor[p0p10_peak_nor$group1 =='Hypo',])
Not=nrow(p0p10_peak_nor[p0p10_peak_nor$group1 =='Not',])
color <- c(Hyper = '#FF6666',Not = alpha("DimGray", 0.2),Hypo ='#0099CC')
ggplot(p0p10_peak_nor, aes(x = p10_m6A, y = p0_m6A,color=group1)) +
  geom_point(alpha=0.5,position = "jitter") + # color="#99CC33"
  scale_color_manual(values = color, breaks=c("Hyper", "Not", "Hypo"), 
                     labels=c(paste0("Hyper (", Hyper, ")"),
                              paste0("Not (", Not, ")"),
                              paste0("Hypo (",Hypo,")")))+
  labs(x=expression(m^6*A~level~'in'~p10), y=expression(m^6*A~level~'in'~p0))+
  theme(legend.position = c(0.19,0.81),
        legend.title = element_blank(),
        legend.text=element_text(size=7))+
  theme(axis.title.x = element_text(size =11))+
  theme(axis.title.y = element_text(size = 11))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) #调整与图片边缘的距离
ggsave('/disk/user_09/Data/03_TC1_caRNA/13_m6A_level/01_noscale/p0p10_volcano.pdf',width = 5,height = 4.5)


p0p10_peak_anno <- as.data.frame(fread("/disk/user_09/Data/03_TC1_caRNA/10_bed_merge/mm10/peak_annotation_homer/p0p10_annotation_homer_rearranged.xls"))
row.names(p0p10_peak_anno) <- p0p10_peak_anno$`PeakID (cmd=annotatePeaks.pl 10_bed_merge/mm10/p0p10_zero_deprived_rearranged.bed mm10 -cpu 10)`

p0p10_peak_nor_mm10 <- p0p10_peak_nor[rownames(p0p10_peak_nor) %in% p0p10_peak_anno$PeakID,]
p0p10_peak_nor_mm10$region <- gsub("\\ .*","",p0p10_peak_anno[rownames(p0p10_peak_nor_mm10),"Annotation"])


p0p10_peak_nor_mm10$region <- ifelse(p0p10_peak_nor_mm10$region=="3'","3'UTR",
                                     ifelse(p0p10_peak_nor_mm10$region=="5'","5'UTR",
                                            p0p10_peak_nor_mm10$region))

color <- c(`3'UTR` = 'red',other = alpha("Gray", 0.2),`5'UTR` ='#739211',exon='Green',intron='#0099CC')
ggplot(p0p10_peak_nor_mm10, aes(x = p10_m6A, y = p0_m6A,color=region)) +
  geom_point(alpha=0.3,position = "jitter") + # color="#99CC33"
  scale_color_manual(values = color, breaks=c("3'UTR", "other", "5'UTR","exon","intron"))+
  labs(x=expression(m^6*A~level~'in'~p10), y=expression(m^6*A~level~'in'~p0))+
  theme(legend.position = c(0.19,0.81),
        legend.title = element_blank(),
        legend.text=element_text(size=11))+
  theme(axis.title.x = element_text(size =11))+
  theme(axis.title.y = element_text(size = 11))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) #调整与图片边缘的距离
ggsave('/disk/user_09/Data/03_TC1_caRNA/13_m6A_level/01_noscale/p0p10_volcano_region.pdf',width = 5,height = 4.5)




volcanoM6APeak <- function(region,peak_nor_mm10anno,stage){
  peak_nor_region <- peak_nor_mm10anno[peak_nor_mm10anno$region==region,9:13]
  color <- c(Hyper = '#FF6666',Not = alpha("DimGray", 0.2),Hypo ='#0099CC')
  ggplot(peak_nor_region, aes(x = peak_nor_region[,paste0(stage[2],"_m6A")], 
                              y = peak_nor_region[,paste0(stage[1],"_m6A")],
                              color=group1)) +
    geom_point(alpha=0.3,position = "jitter") + # color="#99CC33"
    labs(x=paste0("m6A level in ",stage[2]), y=paste0("m6A level in ",stage[1]))+
    scale_color_manual(values = color, breaks=c("Hyper", "Not", "Hypo"), 
                       labels=c(paste0("Hyper (", nrow(peak_nor_region[peak_nor_region$group1=="Hyper",]), ")"),
                                paste0("Not (", nrow(peak_nor_region[peak_nor_region$group1=="Not",]), ")"),
                                paste0("Hypo (",nrow(peak_nor_region[peak_nor_region$group1=="Hypo",]),")")))+
    theme(legend.position = c(0.19,0.81),
          legend.title = element_blank(),
          legend.text=element_text(size=11))+
    theme(axis.title.x = element_text(size =11))+
    theme(axis.title.y = element_text(size = 11))+
    coord_cartesian(xlim=c(min(peak_nor_region[,1],
                               peak_nor_region[,2]),
                           max(peak_nor_region[,1],
                               peak_nor_region[,2])),
                    ylim=c(min(peak_nor_region[,1],
                               peak_nor_region[,2]),
                           max(peak_nor_region[,1],
                               peak_nor_region[,2])))+
    geom_abline()+
    ggtitle(paste0(stage[1],stage[2],"_",region))+
    theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) #调整与图片边缘的距离
  ggsave(paste0("/disk/user_09/Data/03_TC1_caRNA/13_m6A_level/01_noscale/",stage[1],stage[2],"/volcano_",region,".pdf"),width = 5,height = 4.5)
  
}


p0p10_peak_nor_mm10$region <- ifelse(p0p10_peak_nor_mm10$region=="3'","3'UTR",
                                     ifelse(p0p10_peak_nor_mm10$region=="5'","5'UTR",
                                            p0p10_peak_nor_mm10$region))

p10rp2_peak_nor_mm10$region <- ifelse(p10rp2_peak_nor_mm10$region=="3'","3'UTR",
                                      ifelse(p10rp2_peak_nor_mm10$region=="5'","5'UTR",
                                             p10rp2_peak_nor_mm10$region))

p10rp2_peak_nor_mm10$region <- ifelse(p10rp2_peak_nor_mm10$region=="3'","3'UTR",
                                     ifelse(p10rp2_peak_nor_mm10$region=="5'","5'UTR",
                                            p10rp2_peak_nor_mm10$region))

regions=c("3'UTR","5'UTR","exon","Intergenic","intron","promoter-TSS","TTS")

groups=list(c("p0","p5"),c("p0","p10"),c("p5","rp2"),c("p10","rp2"))

dir.create("/disk/user_09/Data/03_TC1_caRNA/13_m6A_level/01_noscale/p0p10")
#for (group in groups){
  for ( region in regions){
    volcanoM6APeak(region,get(paste0(paste0(group[1],group[2]),"_peak_nor_mm10")),group)
  }
#}
