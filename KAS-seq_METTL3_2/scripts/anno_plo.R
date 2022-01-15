GROUP="METTL14"

CTRL_bed1 <- read.delim(paste0("~/KAS-METTL/",GROUP,"/06_macs2/peak_annotation/broad/KAS-seq_",GROUP,"_CTRL_rep1_annotation_homer.xls"),header = T,)
CTRL_bed2 <- read.delim(paste0("~/KAS-METTL/",GROUP,"/06_macs2/peak_annotation/broad/KAS-seq_",GROUP,"_CTRL_rep2_annotation_homer.xls"),header = T)

KO_bed1 <- read.delim(paste0("~/KAS-METTL/",GROUP,"/06_macs2/peak_annotation/broad/KAS-seq_",GROUP,"_KO_rep1_annotation_homer.xls"),header = T)
KO_bed2 <- read.delim(paste0("~/KAS-METTL/",GROUP,"/06_macs2/peak_annotation/broad/KAS-seq_",GROUP,"_KO_rep2_annotation_homer.xls"),header = T)

sampleList <- list(CTRL_bed1,CTRL_bed2,KO_bed1,KO_bed2)

library(magrittr)
library(tidyverse)
library(reshape2)
peak=data.frame()
for(m6a in sampleList){
  peakanno=as.data.frame(m6a)
  anno=gsub("\\(.*).*","",peakanno$Annotation)
  number=as.data.frame(table(anno))
  number$ratio=number$Freq/length(anno)
  sample=strsplit(peakanno[1,1],split='_peak')[[1]][1]
  peak <- data.frame(sample=sample,Exon=number$Freq[1],ExonRatio=number$ratio[1],
                  Intergenic=number$Freq[2],IntergenicRatio=number$ratio[2],
                  Intron=number$Freq[3],IntronRatio=number$ratio[3],
                  Promoter=number$Freq[4],PromoterRatio=number$ratio[4],
                  TTS=number$Freq[5],TTSRatio=number$ratio[5]) %>% rbind(peak, .)
}

mydata=peak[c(1,3,5,7,9,11)]
plotdata=melt(mydata,id.vars=c("sample"),variable.name = "Genomic Regions", value.name = "m6A peaks fraction")
colnames(plotdata)=c("sample", "Genomic Regions","m6A peaks fraction")
plotdata$`Genomic Regions`=gsub('_Ratio|Ratio','',plotdata$`Genomic Regions`)
# plotdata$group=gsub('_.*','',plotdata$sample)
ggplot(plotdata, aes(x=`Genomic Regions`, y=`m6A peaks fraction`,fill=sample))+
  geom_bar(stat = "identity",position = 'dodge',color='white')+
  labs(x="mESC cell line", y=expression('KAS-seq'~'peak'~'fraction'))+
  scale_fill_manual("",values = c('#969696','#BDBDBD', # #B3B3B3
                                  "#E6550D","#FDAE6B"), # #FF6633
                    breaks=c('KAS-seq_METTL14_CTRL_rep1','KAS-seq_METTL14_CTRL_rep2','KAS-seq_METTL14_KO_rep1','KAS-seq_METTL14_KO_rep2'),
                    labels=c('WT-1','WT-2','METTL14_KO-1','METTL14_KO-2'))+
  theme(axis.text = element_text(size = 13), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
        axis.text.x  = element_text(hjust = 1, vjust = 1, angle = 45))+
  theme(legend.position = c(0.14, 0.72),legend.direction = "vertical",
        legend.text = element_text( size = 10))+ 
  theme(axis.title.x = element_text(size = 14,margin = margin(t=15)))+
  theme(axis.title.y = element_text(size = 14,margin = margin(r=11 )))+
  theme_bw()

