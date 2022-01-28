rm(list=ls())

setwd("~/KAS-METTL/METTL3_2/08_FPKM/m6A")
bam=c("~/carRNA_science/m6A-seq/04_bam_raw/CTRL_input_rep1.bam","~/carRNA_science/m6A-seq/04_bam_raw/CTRL_input_rep2.bam")

library(Rsubread)
RNA <- featureCounts(files = bam,
                     annot.ext = '~/reference/annotation/mm19/gencode.vM28.annotation.protein_coding.chr.gtf',
                     isGTFAnnotationFile = TRUE,
                     GTF.featureType = "gene",
                     GTF.attrType = "gene_name",
                     GTF.attrType.extra = "gene_type",
                     useMetaFeatures = FALSE,
                     allowMultiOverlap = FALSE,
                     countMultiMappingReads = FALSE,
                     fracOverlap = 0.5,
                     isPairedEnd=TRUE,
                     verbose = TRUE, 
                     nthreads=50
)

countdata <- RNA[["counts"]]

write.table(countdata,"countdat_LiuJun.txt",quote = FALSE)

count2fpkm <- function(counts,effLen){
  N <- sum(counts)
  exp(log(counts)+log(1e9)-log(effLen)-log(N))
}
fpkms <- as.data.frame(apply(countdata,2,count2fpkm,effLen=RNA[["annotation"]][["Length"]]))

fpkms$mean <- (fpkms$CTRL_input_rep1.bam + fpkms$CTRL_input_rep2.bam)/2

fpkms_silent <- fpkms[fpkms$mean <= 0.1,]

fpkms_other <- fpkms[fpkms$mean > 0.1,]

fpkms_order <- fpkms_other[order(fpkms_other$mean,decreasing = TRUE),]

length <- dim(fpkms_order)[1]
fpkms_high <- fpkms_order[1:2000,]
fpkms_medium <- fpkms_order[(length/2-999):(length/2+1000),]
fpkms_low <- fpkms_order[(length-1999):length,]

random <- sample(1:dim(fpkms_silent)[1],2000)
fpkms_random <- fpkms_silent[random,]

library(rtracklayer)
gtf=import('/disk1/home/user_09/reference/annotation/mm19/gencode.vM28.annotation.protein_coding.chr.gtf')

fpkms_high_gtf=gtf[gtf$gene_name %in% rownames(fpkms_high)]
export(fpkms_high_gtf, "fpkms_high.gtf")

fpkms_medium_gtf=gtf[gtf$gene_name %in% rownames(fpkms_medium)]
export(fpkms_medium_gtf, "fpkms_medium.gtf")
fpkms_low_gtf=gtf[gtf$gene_name %in% rownames(fpkms_low)]
export(fpkms_low_gtf, "fpkms_low.gtf")
fpkms_silent_gtf=gtf[gtf$gene_name %in% rownames(fpkms_random)]
export(fpkms_silent_gtf, "fpkms_silent.gtf")

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggplot2)
data <- data.frame(name=c(rep("high",2000),rep("medium",2000),rep("low",2000),rep("silent",2000)),value=c(fpkms_high[,3],fpkms_medium[,3],fpkms_low[,3],fpkms_random[,3]))

ggplot(data, aes(x=name, y=value, fill=name)) +
geom_boxplot() +
  theme_bw()+
scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
theme(
  legend.position="none",
  plot.title = element_text(size=11)
) +
coord_cartesian(ylim=c(0,20))+
ggtitle("Basic boxplot") +
xlab("")



