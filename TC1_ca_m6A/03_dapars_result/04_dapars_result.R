setwd("/disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_new/")
library(ggplot2)
library(tidyverse)

dapars2_caRNA_mRNA <- as_tibble(read.table("DaPars2_caRNA_mRNA_result.txt",header = TRUE))

colnames(dapars2_caRNA_mRNA) <- c("Gene","fit_value","Predicted_Proximal_APA","Loci","caRNA_p0_rep1_PDUI","caRNA_p0_rep2_PDUI","caRNA_p5_rep1_PDUI","caRNA_p5_rep2_PDUI","caRNA_p10_rep1_PDUI","caRNA_p10_rep2_PDUI","caRNA_rp2_rep1_PDUI","caRNA_rp2_rep2_PDUI","mRNA_p0_rep1_PDUI","mRNA_p0_rep2_PDUI","mRNA_p5_rep1_PDUI","mRNA_p5_rep2_PDUI","mRNA_p10_rep1_PDUI","mRNA_p10_rep2_PDUI","mRNA_rp2_rep1_PDUI","mRNA_rp2_rep2_PDUI")

#write.table(dapars2_caRNA_mRNA,"DaPars2_caRNA_mRNA.txt",quote = FALSE,sep="\t",row.names = FALSE)

dapars2_mean <- dapars2_caRNA_mRNA

dapars2_mean  <- tibble(Gene=dapars2_mean$Gene,caRNA_p0=apply(dapars2_mean[,c("caRNA_p0_rep1_PDUI","caRNA_p0_rep2_PDUI")],1,mean,na.rm=TRUE),
                        caRNA_p5=apply(dapars2_mean[,c("caRNA_p5_rep1_PDUI","caRNA_p5_rep2_PDUI")],1,mean,na.rm=TRUE),
                        caRNA_p10=apply(dapars2_mean[,c("caRNA_p10_rep1_PDUI","caRNA_p10_rep2_PDUI")],1,mean,na.rm=TRUE),
                        caRNA_rp2=apply(dapars2_mean[,c("caRNA_rp2_rep1_PDUI","caRNA_rp2_rep2_PDUI")],1,mean,na.rm=TRUE),
                        mRNA_p0=apply(dapars2_mean[,c("mRNA_p0_rep1_PDUI","mRNA_p0_rep2_PDUI")],1,mean,na.rm=TRUE),
                        mRNA_p5=apply(dapars2_mean[,c("mRNA_p5_rep1_PDUI","mRNA_p5_rep2_PDUI")],1,mean,na.rm=TRUE),
                        mRNA_p10=apply(dapars2_mean[,c("mRNA_p10_rep1_PDUI","mRNA_p10_rep2_PDUI")],1,mean,na.rm=TRUE),
                        mRNA_rp2=apply(dapars2_mean[,c("mRNA_rp2_rep1_PDUI","mRNA_rp2_rep2_PDUI")],1,mean,na.rm=TRUE)) 

dapars2_mean_omitNA <- na.omit(dapars2_mean)

write.table(dapars2_mean_omitNA,"DaPars2_caRNA_mRNA_mean_omitNA.txt",quote = FALSE,sep="\t",row.names = FALSE)

#volcano_color <- c(UP = alpha("#C01623", 0.7),NC = alpha("DimGray", 0.2),DOWN = alpha("#4431A5", 0.7))

## caRNA_p0p10
dapars_caRNA_p0p10 <- select(dapars2_mean_omitNA,c(1,2,4))
dapars_caRNA_p0p10$filter <- "NC"
dapars_caRNA_p0p10$filter[which(abs(dapars_caRNA_p0p10$caRNA_p0-dapars_caRNA_p0p10$caRNA_p10)>0.2 & log2(dapars_caRNA_p0p10$caRNA_p0/dapars_caRNA_p0p10$caRNA_p10)>0.58)] <- "UP"

dapars_caRNA_p0p10$filter[which(abs(dapars_caRNA_p0p10$caRNA_p0-dapars_caRNA_p0p10$caRNA_p10)>0.2 & log2(dapars_caRNA_p0p10$caRNA_p10/dapars_caRNA_p0p10$caRNA_p0)>0.58)] <- "DOWN"

dir.create("group_table_0.2_0.58")
write.table(dapars_caRNA_p0p10,"group_table_0.2_0.58/DaPars2_caRNA_p0p10.txt",quote = FALSE,sep="\t",row.names = FALSE)
table(dapars_caRNA_p0p10$filter)

ggplot(dapars_caRNA_p0p10)+
  geom_point(aes(x=caRNA_p10,y=caRNA_p0,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR p0_longer (", nrow(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR p10_longer (",nrow(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in caRNA_p0") + xlab("Mean PDUIs of genes in caRNA_p10")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_caRNA_p0p10_volcano.pdf",width = 5.5,height = 6)

## caRNA_p0p5
dapars_caRNA_p0p5 <- select(dapars2_mean_omitNA,c(1,2,3))
dapars_caRNA_p0p5$filter <- "NC"
dapars_caRNA_p0p5$filter[which(abs(dapars_caRNA_p0p5$caRNA_p0-dapars_caRNA_p0p5$caRNA_p5)>0.2 & log2(dapars_caRNA_p0p5$caRNA_p0/dapars_caRNA_p0p5$caRNA_p5)>0.58)] <- "UP"

dapars_caRNA_p0p5$filter[which(abs(dapars_caRNA_p0p5$caRNA_p0-dapars_caRNA_p0p5$caRNA_p5)>0.2 & log2(dapars_caRNA_p0p5$caRNA_p5/dapars_caRNA_p0p5$caRNA_p0)>0.58)] <- "DOWN"

table(dapars_caRNA_p0p5$filter)
write.table(dapars_caRNA_p0p5,"group_table_0.2_0.58/DaPars2_caRNA_p0p5.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_caRNA_p0p5)+
  geom_point(aes(x=caRNA_p5,y=caRNA_p0,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR p0_longer (", nrow(dapars_caRNA_p0p5[dapars_caRNA_p0p5$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR p5_longer (",nrow(dapars_caRNA_p0p5[dapars_caRNA_p0p5$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in caRNA_p0") + xlab("Mean PDUIs of genes in caRNA_p5")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_caRNA_p0p5_volcano.pdf",width = 5.5,height = 6)

## caRNA_p0rp2
dapars_caRNA_p0rp2 <- select(dapars2_mean_omitNA,c(1,2,5))
dapars_caRNA_p0rp2$filter <- "NC"
dapars_caRNA_p0rp2$filter[which(abs(dapars_caRNA_p0rp2$caRNA_p0-dapars_caRNA_p0rp2$caRNA_rp2)>0.2 & log2(dapars_caRNA_p0rp2$caRNA_p0/dapars_caRNA_p0rp2$caRNA_rp2)>0.58)] <- "UP"

dapars_caRNA_p0rp2$filter[which(abs(dapars_caRNA_p0rp2$caRNA_p0-dapars_caRNA_p0rp2$caRNA_rp2)>0.2 & log2(dapars_caRNA_p0rp2$caRNA_rp2/dapars_caRNA_p0rp2$caRNA_p0)>0.58)] <- "DOWN"

table(dapars_caRNA_p0rp2$filter)

write.table(dapars_caRNA_p0rp2,"group_table_0.2_0.58/DaPars2_caRNA_p0rp2.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_caRNA_p0rp2)+
  geom_point(aes(x=caRNA_rp2,y=caRNA_p0,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR p0_longer (", nrow(dapars_caRNA_p0rp2[dapars_caRNA_p0rp2$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR rp2_longer (",nrow(dapars_caRNA_p0rp2[dapars_caRNA_p0rp2$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in caRNA_p0") + xlab("Mean PDUIs of genes in caRNA_rp2")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_caRNA_p0rp2_volcano.pdf",width = 5.5,height = 6)

## caRNA_p5rp2
dapars_caRNA_p5rp2 <- select(dapars2_mean_omitNA,c(1,3,5))
dapars_caRNA_p5rp2$filter <- "NC"
dapars_caRNA_p5rp2$filter[which(abs(dapars_caRNA_p5rp2$caRNA_p5-dapars_caRNA_p5rp2$caRNA_rp2)>0.2 & log2(dapars_caRNA_p5rp2$caRNA_p5/dapars_caRNA_p5rp2$caRNA_rp2)>0.58)] <- "UP"

dapars_caRNA_p5rp2$filter[which(abs(dapars_caRNA_p5rp2$caRNA_p5-dapars_caRNA_p5rp2$caRNA_rp2)>0.2 & log2(dapars_caRNA_p5rp2$caRNA_rp2/dapars_caRNA_p5rp2$caRNA_p5)>0.58)] <- "DOWN"

table(dapars_caRNA_p5rp2$filter)

write.table(dapars_caRNA_p5rp2,"group_table_0.2_0.58/DaPars2_caRNA_p5rp2.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_caRNA_p5rp2)+
  geom_point(aes(x=caRNA_rp2,y=caRNA_p5,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR p5_longer (", nrow(dapars_caRNA_p5rp2[dapars_caRNA_p5rp2$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR rp2_longer (",nrow(dapars_caRNA_p5rp2[dapars_caRNA_p5rp2$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in caRNA_p5") + xlab("Mean PDUIs of genes in caRNA_rp2")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_caRNA_p5rp2_volcano.pdf",width = 5.5,height = 6)

## caRNA_p10rp2
dapars_caRNA_p10rp2 <- select(dapars2_mean_omitNA,c(1,4,5))
dapars_caRNA_p10rp2$filter <- "NC"
dapars_caRNA_p10rp2$filter[which(abs(dapars_caRNA_p10rp2$caRNA_p10-dapars_caRNA_p10rp2$caRNA_rp2)>0.2 & log2(dapars_caRNA_p10rp2$caRNA_p10/dapars_caRNA_p10rp2$caRNA_rp2)>0.58)] <- "UP"

dapars_caRNA_p10rp2$filter[which(abs(dapars_caRNA_p10rp2$caRNA_p10-dapars_caRNA_p10rp2$caRNA_rp2)>0.2 & log2(dapars_caRNA_p10rp2$caRNA_rp2/dapars_caRNA_p10rp2$caRNA_p10)>0.58)] <- "DOWN"

table(dapars_caRNA_p10rp2$filter)

write.table(dapars_caRNA_p10rp2,"group_table_0.2_0.58/DaPars2_caRNA_p10rp2.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_caRNA_p10rp2)+
  geom_point(aes(x=caRNA_rp2,y=caRNA_p10,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR p10_longer (", nrow(dapars_caRNA_p10rp2[dapars_caRNA_p10rp2$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR rp2_longer (",nrow(dapars_caRNA_p10rp2[dapars_caRNA_p10rp2$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in caRNA_p10") + xlab("Mean PDUIs of genes in caRNA_rp2")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_caRNA_p10rp2_volcano.pdf",width = 5.5,height = 6)

## caRNA_p5p10
dapars_caRNA_p5p10 <- select(dapars2_mean_omitNA,c(1,3,4))
dapars_caRNA_p5p10$filter <- "NC"
dapars_caRNA_p5p10$filter[which(abs(dapars_caRNA_p5p10$caRNA_p5-dapars_caRNA_p5p10$caRNA_p10)>0.2 & log2(dapars_caRNA_p5p10$caRNA_p5/dapars_caRNA_p5p10$caRNA_p10)>0.58)] <- "UP"

dapars_caRNA_p5p10$filter[which(abs(dapars_caRNA_p5p10$caRNA_p5-dapars_caRNA_p5p10$caRNA_p10)>0.2 & log2(dapars_caRNA_p5p10$caRNA_p10/dapars_caRNA_p5p10$caRNA_p5)>0.58)] <- "DOWN"

table(dapars_caRNA_p5p10$filter)

write.table(dapars_caRNA_p5p10,"group_table_0.2_0.58/DaPars2_caRNA_p5p10.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_caRNA_p5p10)+
  geom_point(aes(x=caRNA_p10,y=caRNA_p5,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR p5_longer (", nrow(dapars_caRNA_p5p10[dapars_caRNA_p5p10$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR p10_longer (",nrow(dapars_caRNA_p5p10[dapars_caRNA_p5p10$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in caRNA_p5") + xlab("Mean PDUIs of genes in caRNA_p10")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_caRNA_p5p10_volcano.pdf",width = 5.5,height = 6)

library(VennDiagram)
setwd("/disk/user_09/Data/03_TC1_caRNA/12_dapars/dapars2_new/")

dir.create("venn_diagram_0.2_0.58")

venn.diagram(x=list(p0p5_shorter=dapars_caRNA_p0p5$Gene[which(dapars_caRNA_p0p5$filter=="UP")],
                    p5rp2_longer=dapars_caRNA_p5rp2$Gene[which(dapars_caRNA_p5rp2$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p0p5rp2_caRNA_shorter_longer.tiff")

venn.diagram(x=list(p0p5_longer=dapars_caRNA_p0p5$Gene[which(dapars_caRNA_p0p5$filter=="DOWN")],
                    p5rp2_shorter=dapars_caRNA_p5rp2$Gene[which(dapars_caRNA_p5rp2$filter=="UP")]),
             "venn_diagram_0.2_0.58/p0p5rp2_caRNA_longer_shorter.tiff")

venn.diagram(x=list(p0p10_shorter=dapars_caRNA_p0p10$Gene[which(dapars_caRNA_p0p10$filter=="UP")],
                    p10rp2_longer=dapars_caRNA_p10rp2$Gene[which(dapars_caRNA_p10rp2$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p0p10rp2_caRNA_shorter_longer.tiff")

venn.diagram(x=list(p0p10_longer=dapars_caRNA_p0p10$Gene[which(dapars_caRNA_p0p10$filter=="DOWN")],
                    p10rp2_shorter=dapars_caRNA_p10rp2$Gene[which(dapars_caRNA_p10rp2$filter=="UP")]),
             "venn_diagram_0.2_0.58/p0p10rp2_caRNA_longer_shorter.tiff")

venn.diagram(x=list(p0p10_UP=dapars_caRNA_p0p10$Gene[which(dapars_caRNA_p0p10$filter=="UP")],
                    p0p10_DOWN=dapars_caRNA_p0p10$Gene[which(dapars_caRNA_p0p10$filter=="DOWN")],
                    p10rp2_UP=dapars_caRNA_p10rp2$Gene[which(dapars_caRNA_p10rp2$filter=="UP")],
                    p10rp2_DOWN=dapars_caRNA_p10rp2$Gene[which(dapars_caRNA_p10rp2$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p0p10rp2_caRNA.tiff")

venn.diagram(x=list(p0p5_shorter=dapars_caRNA_p0p5$Gene[which(dapars_caRNA_p0p5$filter=="UP")],
                    p0p10_shorter=dapars_caRNA_p0p10$Gene[which(dapars_caRNA_p0p10$filter=="UP")]),
             "venn_diagram_0.2_0.58/p0p5p10_caRNA_shorter.tiff")

venn.diagram(x=list(p0p5_longer=dapars_caRNA_p0p5$Gene[which(dapars_caRNA_p0p5$filter=="DOWN")],
                    p0p10_longer=dapars_caRNA_p0p10$Gene[which(dapars_caRNA_p0p10$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p0p5p10_caRNA_longer.tiff")

##############################################

## mRNA_p0p10
dapars_mRNA_p0p10 <- select(dapars2_mean_omitNA,c(1,6,8))
dapars_mRNA_p0p10$filter <- "NC"
dapars_mRNA_p0p10$filter[which(abs(dapars_mRNA_p0p10$mRNA_p0-dapars_mRNA_p0p10$mRNA_p10)>0.2 & log2(dapars_mRNA_p0p10$mRNA_p0/dapars_mRNA_p0p10$mRNA_p10)>0.58)] <- "UP"

dapars_mRNA_p0p10$filter[which(abs(dapars_mRNA_p0p10$mRNA_p0-dapars_mRNA_p0p10$mRNA_p10)>0.2 & log2(dapars_mRNA_p0p10$mRNA_p10/dapars_mRNA_p0p10$mRNA_p0)>0.58)] <- "DOWN"

table(dapars_mRNA_p0p10$filter)

write.table(dapars_mRNA_p0p10,"group_table_0.2_0.58/DaPars2_mRNA_p0p10.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_mRNA_p0p10)+
  geom_point(aes(x=mRNA_p10,y=mRNA_p0,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR p0_longer (", nrow(dapars_mRNA_p0p10[dapars_mRNA_p0p10$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR p10_longer (",nrow(dapars_mRNA_p0p10[dapars_mRNA_p0p10$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in mRNA_p0") + xlab("Mean PDUIs of genes in mRNA_p10")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_mRNA_p0p10_volcano.pdf",width = 5.5,height = 6)

## mRNA_p0p5
dapars_mRNA_p0p5 <- select(dapars2_mean_omitNA,c(1,6,7))
dapars_mRNA_p0p5$filter <- "NC"
dapars_mRNA_p0p5$filter[which(abs(dapars_mRNA_p0p5$mRNA_p0-dapars_mRNA_p0p5$mRNA_p5)>0.2 & log2(dapars_mRNA_p0p5$mRNA_p0/dapars_mRNA_p0p5$mRNA_p5)>0.58)] <- "UP"

dapars_mRNA_p0p5$filter[which(abs(dapars_mRNA_p0p5$mRNA_p0-dapars_mRNA_p0p5$mRNA_p5)>0.2 & log2(dapars_mRNA_p0p5$mRNA_p5/dapars_mRNA_p0p5$mRNA_p0)>0.58)] <- "DOWN"

table(dapars_mRNA_p0p5$filter)
write.table(dapars_mRNA_p0p5,"group_table_0.2_0.58/DaPars2_mRNA_p0p5.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_mRNA_p0p5)+
  geom_point(aes(x=mRNA_p5,y=mRNA_p0,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR p0_longer (", nrow(dapars_mRNA_p0p5[dapars_mRNA_p0p5$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR p5_longer (",nrow(dapars_mRNA_p0p5[dapars_mRNA_p0p5$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in mRNA_p0") + xlab("Mean PDUIs of genes in mRNA_p5")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_mRNA_p0p5_volcano.pdf",width = 5.5,height = 6)

## mRNA_p0rp2
dapars_mRNA_p0rp2 <- select(dapars2_mean_omitNA,c(1,6,9))
dapars_mRNA_p0rp2$filter <- "NC"
dapars_mRNA_p0rp2$filter[which(abs(dapars_mRNA_p0rp2$mRNA_p0-dapars_mRNA_p0rp2$mRNA_rp2)>0.2 & log2(dapars_mRNA_p0rp2$mRNA_p0/dapars_mRNA_p0rp2$mRNA_rp2)>0.58)] <- "UP"

dapars_mRNA_p0rp2$filter[which(abs(dapars_mRNA_p0rp2$mRNA_p0-dapars_mRNA_p0rp2$mRNA_rp2)>0.2 & log2(dapars_mRNA_p0rp2$mRNA_rp2/dapars_mRNA_p0rp2$mRNA_p0)>0.58)] <- "DOWN"

table(dapars_mRNA_p0rp2$filter)
write.table(dapars_mRNA_p0rp2,"group_table_0.2_0.58/DaPars2_mRNA_p0rp2.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_mRNA_p0rp2)+
  geom_point(aes(x=mRNA_rp2,y=mRNA_p0,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR p0_longer (", nrow(dapars_mRNA_p0rp2[dapars_mRNA_p0rp2$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR rp2_longer (",nrow(dapars_mRNA_p0rp2[dapars_mRNA_p0rp2$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in mRNA_p0") + xlab("Mean PDUIs of genes in mRNA_rp2")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_mRNA_p0rp2_volcano.pdf",width = 5.5,height = 6)

## mRNA_p5rp2
dapars_mRNA_p5rp2 <- select(dapars2_mean_omitNA,c(1,7,9))
dapars_mRNA_p5rp2$filter <- "NC"
dapars_mRNA_p5rp2$filter[which(abs(dapars_mRNA_p5rp2$mRNA_p5-dapars_mRNA_p5rp2$mRNA_rp2)>0.2 & log2(dapars_mRNA_p5rp2$mRNA_p5/dapars_mRNA_p5rp2$mRNA_rp2)>0.58)] <- "UP"

dapars_mRNA_p5rp2$filter[which(abs(dapars_mRNA_p5rp2$mRNA_p5-dapars_mRNA_p5rp2$mRNA_rp2)>0.2 & log2(dapars_mRNA_p5rp2$mRNA_rp2/dapars_mRNA_p5rp2$mRNA_p5)>0.58)] <- "DOWN"

table(dapars_mRNA_p5rp2$filter)
write.table(dapars_mRNA_p5rp2,"group_table_0.2_0.58/DaPars2_mRNA_p5rp2.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_mRNA_p5rp2)+
  geom_point(aes(x=mRNA_rp2,y=mRNA_p5,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR p5_longer (", nrow(dapars_mRNA_p5rp2[dapars_mRNA_p5rp2$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR rp2_longer (",nrow(dapars_mRNA_p5rp2[dapars_mRNA_p5rp2$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in mRNA_p5") + xlab("Mean PDUIs of genes in mRNA_rp2")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_mRNA_p5rp2_volcano.pdf",width = 5.5,height = 6)

## mRNA_p10rp2
dapars_mRNA_p10rp2 <- select(dapars2_mean_omitNA,c(1,8,9))
dapars_mRNA_p10rp2$filter <- "NC"
dapars_mRNA_p10rp2$filter[which(abs(dapars_mRNA_p10rp2$mRNA_p10-dapars_mRNA_p10rp2$mRNA_rp2)>0.2 & log2(dapars_mRNA_p10rp2$mRNA_p10/dapars_mRNA_p10rp2$mRNA_rp2)>0.58)] <- "UP"

dapars_mRNA_p10rp2$filter[which(abs(dapars_mRNA_p10rp2$mRNA_p10-dapars_mRNA_p10rp2$mRNA_rp2)>0.2 & log2(dapars_mRNA_p10rp2$mRNA_rp2/dapars_mRNA_p10rp2$mRNA_p10)>0.58)] <- "DOWN"

table(dapars_mRNA_p10rp2$filter)
write.table(dapars_mRNA_p10rp2,"group_table_0.2_0.58/DaPars2_mRNA_p10rp2.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_mRNA_p10rp2)+
  geom_point(aes(x=mRNA_rp2,y=mRNA_p10,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR p10_longer (", nrow(dapars_mRNA_p10rp2[dapars_mRNA_p10rp2$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR rp2_longer (",nrow(dapars_mRNA_p10rp2[dapars_mRNA_p10rp2$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in mRNA_p10") + xlab("Mean PDUIs of genes in mRNA_rp2")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_mRNA_p10rp2_volcano.pdf",width = 5.5,height = 6)

## mRNA_p5p10
dapars_mRNA_p5p10 <- select(dapars2_mean_omitNA,c(1,7,8))
dapars_mRNA_p5p10$filter <- "NC"
dapars_mRNA_p5p10$filter[which(abs(dapars_mRNA_p5p10$mRNA_p5-dapars_mRNA_p5p10$mRNA_p10)>0.2 & log2(dapars_mRNA_p5p10$mRNA_p5/dapars_mRNA_p5p10$mRNA_p10)>0.58)] <- "UP"

dapars_mRNA_p5p10$filter[which(abs(dapars_mRNA_p5p10$mRNA_p5-dapars_mRNA_p5p10$mRNA_p10)>0.2 & log2(dapars_mRNA_p5p10$mRNA_p10/dapars_mRNA_p5p10$mRNA_p5)>0.58)] <- "DOWN"

table(dapars_mRNA_p5p10$filter)
write.table(dapars_mRNA_p5p10,"group_table_0.2_0.58/DaPars2_mRNA_p5p10.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_mRNA_p5p10)+
  geom_point(aes(x=mRNA_p10,y=mRNA_p5,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR p5_longer (", nrow(dapars_mRNA_p5p10[dapars_mRNA_p5p10$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR p10_longer (",nrow(dapars_mRNA_p5p10[dapars_mRNA_p5p10$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in mRNA_p5") + xlab("Mean PDUIs of genes in mRNA_p10")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_mRNA_p5p10_volcano.pdf",width = 5.5,height = 6)

venn.diagram(x=list(p0p5_shorter=dapars_mRNA_p0p5$Gene[which(dapars_mRNA_p0p5$filter=="UP")],
                    p5rp2_longer=dapars_mRNA_p5rp2$Gene[which(dapars_mRNA_p5rp2$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p0p5rp2_mRNA_shorter_longer.tiff")

venn.diagram(x=list(p0p5_longer=dapars_mRNA_p0p5$Gene[which(dapars_mRNA_p0p5$filter=="DOWN")],
                    p5rp2_shorter=dapars_mRNA_p5rp2$Gene[which(dapars_mRNA_p5rp2$filter=="UP")]),
             "venn_diagram_0.2_0.58/p0p5rp2_mRNA_longer_shorter.tiff")

venn.diagram(x=list(p0p10_shorter=dapars_mRNA_p0p10$Gene[which(dapars_mRNA_p0p10$filter=="UP")],
                    p10rp2_longer=dapars_mRNA_p10rp2$Gene[which(dapars_mRNA_p10rp2$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p0p10rp2_mRNA_shorter_longer.tiff")

venn.diagram(x=list(p0p10_longer=dapars_mRNA_p0p10$Gene[which(dapars_mRNA_p0p10$filter=="DOWN")],
                    p10rp2_shorter=dapars_mRNA_p10rp2$Gene[which(dapars_mRNA_p10rp2$filter=="UP")]),
             "venn_diagram_0.2_0.58/p0p10rp2_mRNA_longer_shorter.tiff")

################################################

venn.diagram(x=list(p0p5_mRNA_shorter=dapars_mRNA_p0p5$Gene[which(dapars_mRNA_p0p5$filter=="UP")],
                    p0p5_mRNA__longer=dapars_mRNA_p0p5$Gene[which(dapars_mRNA_p0p5$filter=="DOWN")],
                    p0p5_caRNA_shorter=dapars_caRNA_p0p5$Gene[which(dapars_caRNA_p0p5$filter=="UP")],
                    p0p5_caRNA_longer=dapars_caRNA_p0p5$Gene[which(dapars_caRNA_p0p5$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p0p5_mRNA_caRNA.tiff")

venn.diagram(x=list(p0p10_mRNA_shorter=dapars_mRNA_p0p10$Gene[which(dapars_mRNA_p0p10$filter=="UP")],
                    p0p10_mRNA__longer=dapars_mRNA_p0p10$Gene[which(dapars_mRNA_p0p10$filter=="DOWN")],
                    p0p10_caRNA_shorter=dapars_caRNA_p0p10$Gene[which(dapars_caRNA_p0p10$filter=="UP")],
                    p0p10_caRNA_longer=dapars_caRNA_p0p10$Gene[which(dapars_caRNA_p0p10$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p0p10_mRNA_caRNA.tiff")

venn.diagram(x=list(p0rp2_mRNA_shorter=dapars_mRNA_p0rp2$Gene[which(dapars_mRNA_p0rp2$filter=="UP")],
                    p0rp2_mRNA__longer=dapars_mRNA_p0rp2$Gene[which(dapars_mRNA_p0rp2$filter=="DOWN")],
                    p0rp2_caRNA_shorter=dapars_caRNA_p0rp2$Gene[which(dapars_caRNA_p0rp2$filter=="UP")],
                    p0rp2_caRNA_longer=dapars_caRNA_p0rp2$Gene[which(dapars_caRNA_p0rp2$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p0rp2_mRNA_caRNA.tiff")

venn.diagram(x=list(p5p10_mRNA_shorter=dapars_mRNA_p5p10$Gene[which(dapars_mRNA_p5p10$filter=="UP")],
                    p5p10_mRNA__longer=dapars_mRNA_p5p10$Gene[which(dapars_mRNA_p5p10$filter=="DOWN")],
                    p5p10_caRNA_shorter=dapars_caRNA_p5p10$Gene[which(dapars_caRNA_p5p10$filter=="UP")],
                    p5p10_caRNA_longer=dapars_caRNA_p5p10$Gene[which(dapars_caRNA_p5p10$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p5p10_mRNA_caRNA.tiff")

venn.diagram(x=list(p5rp2_mRNA_shorter=dapars_mRNA_p5rp2$Gene[which(dapars_mRNA_p5rp2$filter=="UP")],
                    p5rp2_mRNA__longer=dapars_mRNA_p5rp2$Gene[which(dapars_mRNA_p5rp2$filter=="DOWN")],
                    p5rp2_caRNA_shorter=dapars_caRNA_p5rp2$Gene[which(dapars_caRNA_p5rp2$filter=="UP")],
                    p5rp2_caRNA_longer=dapars_caRNA_p5rp2$Gene[which(dapars_caRNA_p5rp2$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p5rp2_mRNA_caRNA.tiff")

venn.diagram(x=list(p10rp2_mRNA_shorter=dapars_mRNA_p10rp2$Gene[which(dapars_mRNA_p10rp2$filter=="UP")],
                    p10rp2_mRNA__longer=dapars_mRNA_p10rp2$Gene[which(dapars_mRNA_p10rp2$filter=="DOWN")],
                    p10rp2_caRNA_shorter=dapars_caRNA_p10rp2$Gene[which(dapars_caRNA_p10rp2$filter=="UP")],
                    p10rp2_caRNA_longer=dapars_caRNA_p10rp2$Gene[which(dapars_caRNA_p10rp2$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p10rp2_mRNA_caRNA.tiff")


venn.diagram(x=list(p0p5_longer=dapars_mRNA_p0p5$Gene[which(dapars_mRNA_p0p5$filter=="DOWN")],
                    p5rp2_shorter=dapars_mRNA_p5rp2$Gene[which(dapars_mRNA_p5rp2$filter=="UP")]),
             "venn_diagram_0.2_0.58/p0p5rp2_mRNA_longer_shorter.tiff")

venn.diagram(x=list(p0p10_shorter=dapars_mRNA_p0p10$Gene[which(dapars_mRNA_p0p10$filter=="UP")],
                    p10rp2_longer=dapars_mRNA_p10rp2$Gene[which(dapars_mRNA_p10rp2$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p0p10rp2_mRNA_shorter_longer.tiff")

venn.diagram(x=list(p0p10_longer=dapars_mRNA_p0p10$Gene[which(dapars_mRNA_p0p10$filter=="DOWN")],
                    p10rp2_shorter=dapars_mRNA_p10rp2$Gene[which(dapars_mRNA_p10rp2$filter=="UP")]),
             "venn_diagram_0.2_0.58/p0p10rp2_mRNA_longer_shorter.tiff")


venn.diagram(x=list(p0p10_caRNA_shorter=dapars_caRNA_p0p10$Gene[which(dapars_caRNA_p0p10$filter=="UP")],
                    p10rp2_caRNA_longer=dapars_caRNA_p10rp2$Gene[which(dapars_caRNA_p10rp2$filter=="DOWN")],
                    p0p5_caRNA_shorter=dapars_caRNA_p0p5$Gene[which(dapars_caRNA_p0p5$filter=="UP")],
                    p5rp2_caRNA_longer=dapars_caRNA_p5rp2$Gene[which(dapars_caRNA_p5rp2$filter=="DOWN")]),
             "venn_diagram_0.2_0.58/p0p5p10rp2_caRNA.tiff")


###############################################

## p0_caRNA_mRNA
dapars_caRNAmRNA <- select(dapars2_mean_omitNA,c(1,2,6))
dapars_caRNAmRNA$filter <- "NC"
dapars_caRNAmRNA$filter[which(abs(dapars_caRNAmRNA$caRNA_p0-dapars_caRNAmRNA$mRNA_p0)>0.2 & log2(dapars_caRNAmRNA$caRNA_p0/dapars_caRNAmRNA$mRNA_p0)>0.58)] <- "UP"

dapars_caRNAmRNA$filter[which(abs(dapars_caRNAmRNA$caRNA_p0-dapars_caRNAmRNA$mRNA_p0)>0.2 & log2(dapars_caRNAmRNA$mRNA_p0/dapars_caRNAmRNA$caRNA_p0)>0.58)] <- "DOWN"

table(dapars_caRNAmRNA$filter)
write.table(dapars_caRNAmRNA,"group_table_0.2_0.58/DaPars2_p0_caRNA_mRNA.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_caRNAmRNA)+
  geom_point(aes(x=mRNA_p0,y=caRNA_p0,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR caRNA_longer (", nrow(dapars_caRNAmRNA[dapars_caRNAmRNA$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR mRNA_longer (",nrow(dapars_caRNAmRNA[dapars_caRNAmRNA$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in caRNA_p0") + xlab("Mean PDUIs of genes in mRNA_p0")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_p0_caRNA_mRNA_volcano.pdf",width = 5.5,height = 6)

## p5_caRNA_mRNA
dapars_caRNAmRNA <- select(dapars2_mean_omitNA,c(1,3,7))
dapars_caRNAmRNA$filter <- "NC"
dapars_caRNAmRNA$filter[which(abs(dapars_caRNAmRNA$caRNA_p5-dapars_caRNAmRNA$mRNA_p5)>0.2 & log2(dapars_caRNAmRNA$caRNA_p5/dapars_caRNAmRNA$mRNA_p5)>0.58)] <- "UP"

dapars_caRNAmRNA$filter[which(abs(dapars_caRNAmRNA$caRNA_p5-dapars_caRNAmRNA$mRNA_p5)>0.2 & log2(dapars_caRNAmRNA$mRNA_p5/dapars_caRNAmRNA$caRNA_p5)>0.58)] <- "DOWN"

table(dapars_caRNAmRNA$filter)
write.table(dapars_caRNAmRNA,"group_table_0.2_0.58/DaPars2_p5_caRNA_mRNA.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_caRNAmRNA)+
  geom_point(aes(x=mRNA_p5,y=caRNA_p5,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR caRNA_longer (", nrow(dapars_caRNAmRNA[dapars_caRNAmRNA$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR mRNA_longer (",nrow(dapars_caRNAmRNA[dapars_caRNAmRNA$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in caRNA_p5") + xlab("Mean PDUIs of genes in mRNA_p5")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_p5_caRNA_mRNA_volcano.pdf",width = 5.5,height = 6)

## p10_caRNA_mRNA
dapars_caRNAmRNA <- select(dapars2_mean_omitNA,c(1,4,8))
dapars_caRNAmRNA$filter <- "NC"
dapars_caRNAmRNA$filter[which(abs(dapars_caRNAmRNA$caRNA_p10-dapars_caRNAmRNA$mRNA_p10)>0.2 & log2(dapars_caRNAmRNA$caRNA_p10/dapars_caRNAmRNA$mRNA_p10)>0.58)] <- "UP"

dapars_caRNAmRNA$filter[which(abs(dapars_caRNAmRNA$caRNA_p10-dapars_caRNAmRNA$mRNA_p10)>0.2 & log2(dapars_caRNAmRNA$mRNA_p10/dapars_caRNAmRNA$caRNA_p10)>0.58)] <- "DOWN"

table(dapars_caRNAmRNA$filter)
write.table(dapars_caRNAmRNA,"group_table_0.2_0.58/DaPars2_p10_caRNA_mRNA.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_caRNAmRNA)+
  geom_point(aes(x=mRNA_p10,y=caRNA_p10,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR caRNA_longer (", nrow(dapars_caRNAmRNA[dapars_caRNAmRNA$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR mRNA_longer (",nrow(dapars_caRNAmRNA[dapars_caRNAmRNA$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in caRNA_p10") + xlab("Mean PDUIs of genes in mRNA_p10")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_p10_caRNA_mRNA_volcano.pdf",width = 5.5,height = 6)


## rp2_caRNA_mRNA
dapars_caRNAmRNA <- select(dapars2_mean_omitNA,c(1,5,9))
dapars_caRNAmRNA$filter <- "NC"
dapars_caRNAmRNA$filter[which(abs(dapars_caRNAmRNA$caRNA_rp2-dapars_caRNAmRNA$mRNA_rp2)>0.2 & log2(dapars_caRNAmRNA$caRNA_rp2/dapars_caRNAmRNA$mRNA_rp2)>0.58)] <- "UP"

dapars_caRNAmRNA$filter[which(abs(dapars_caRNAmRNA$caRNA_rp2-dapars_caRNAmRNA$mRNA_rp2)>0.2 & log2(dapars_caRNAmRNA$mRNA_rp2/dapars_caRNAmRNA$caRNA_rp2)>0.58)] <- "DOWN"

table(dapars_caRNAmRNA$filter)
write.table(dapars_caRNAmRNA,"group_table_0.2_0.58/DaPars2_rp2_caRNA_mRNA.txt",quote = FALSE,sep="\t",row.names = FALSE)

ggplot(dapars_caRNAmRNA)+
  geom_point(aes(x=mRNA_rp2,y=caRNA_rp2,color=filter),size=2,position="jitter",alpha=0.5)+
  scale_color_manual(values = volcano_color, breaks=c("UP", "NC", "DOWN"),labels=c(paste0("3'UTR caRNA_longer (", nrow(dapars_caRNAmRNA[dapars_caRNAmRNA$filter == "UP",]), ")"), "3'UTR non-significant", paste0("3'UTR mRNA_longer (",nrow(dapars_caRNAmRNA[dapars_caRNAmRNA$filter == "DOWN",]),")")))+
  ylab("Mean PDUIs of genes in caRNA_rp2") + xlab("Mean PDUIs of genes in mRNA_rp2")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        #legend.position = c(0.85,0.75),
        legend.background = element_blank(),
        panel.grid =element_blank(),
        panel.background = element_rect(fill = "white",colour="black",size=2),
        legend.key = element_blank(),
        legend.text = element_text(size = 15,  face = 'bold'),
        legend.direction= "vertical")+
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.3), "cm"))+ #调整与图片边缘的距离
  theme(axis.title.x = element_text(size = 18,margin = margin(t=8)))+
  theme(axis.title.y = element_text(size = 18,margin = margin(r=5 )))
ggsave("group_table_0.2_0.58/dapars_rp2_caRNA_mRNA_volcano.pdf",width = 5.5,height = 6)




library(VennDiagram)
setwd("/disk1/home/user_08/Data/TC1-planB/05_bam_change_index_hisat2/01_bam_sorted/12_dapars2/DaPars2_caRNA_mRNA_protein_coding/")

dir.create("venn_diagram_0.2_0.58")
venn.diagram(x=list(p0p5_UP=dapars_p0p5$Gene[which(dapars_p0p5$filter=="UP")],p0p5_DOWN=dapars_p0p5$Gene[which(dapars_p0p5$filter=="DOWN")],p5rp2_UP=dapars_p5rp2$Gene[which(dapars_p5rp2$filter=="UP")],p5rp2_DOWN=dapars_p5rp2$Gene[which(dapars_p5rp2$filter=="DOWN")]),"venn_diagram_0.2_0.58/p0p5rp2_mRNA.tiff")

ggplot(dapars_p0p10)+
  geom_histogram(aes(X01_wiggle.TC1_RNAseq_p0_PDUI),color="red",fill="red",binwidth = 0.01,alpha=0.3)+
  geom_histogram(aes(X01_wiggle.TC1_RNAseq_p6_PDUI),color="blue",fill="blue",binwidth = 0.01,alpha=0.3)

ggplot(dapars_p10rp2)+
  geom_histogram(aes(X01_wiggle.TC1_RNAseq_p6_PDUI),color="red",fill="red",binwidth = 0.01,alpha=0.3)+
  geom_histogram(aes(X01_wiggle.TC1_RNAseq_rp2_PDUI),color="blue",fill="blue",binwidth = 0.01,alpha=0.3)
