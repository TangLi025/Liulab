rm(list=ls())
library(ChIPpeakAnno)

MODE <- "broad"
GROUP <- "KAS-seq_ALKBH5"
TYPE <- "broadPeak"

CTRL_bed1 <- "/disk1/home/user_08/custom/HepG2-caJ_rep1_peaks.bed"
CTRL_bed2 <- "/disk1/home/user_08/custom/HepG2-caJ_rep2_peaks.bed"

KO_bed1 <- "/disk1/home/user_08/custom/HepG2-toN_rep1_peaks.bed"
KO_bed2 <- "/disk1/home/user_08/custom/HepG2-toN_rep2_peaks.bed"

CTRL_1 <- ChIPpeakAnno::toGRanges(CTRL_bed1, format="BED", header=FALSE)
CTRL_2 <- ChIPpeakAnno::toGRanges(CTRL_bed2, format="BED", header=FALSE)
KO_1 <- ChIPpeakAnno::toGRanges(KO_bed1, format="BED", header=FALSE)
KO_2 <- ChIPpeakAnno::toGRanges(KO_bed2, format="BED", header=FALSE)


## must keep the class exactly same as gr1$score, i.e., numeric.
CTRL_1$score <- as.numeric(CTRL_1$score) 
CTRL_2$score <- as.numeric(CTRL_2$score) 
KO_1$score <- as.numeric(KO_1$score) 
KO_2$score <- as.numeric(KO_2$score) 

ol_CTRL <- findOverlapsOfPeaks(CTRL_1, CTRL_2)
## add metadata (mean of score) to the overlapping peaks
ol_CTRL <- addMetadata(ol_CTRL, colNames="score", FUN=base::mean) 
ol_CTRL$peaklist[["CTRL_1///CTRL_2"]]
makeVennDiagram(ol_CTRL, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 


ol_KO <- findOverlapsOfPeaks(KO_1, KO_2)
ol_KO <- addMetadata(ol_KO, colNames="score", FUN=base::mean) 
makeVennDiagram(ol_KO, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 



CTRL <- ol_CTRL$peaklist[["CTRL_1///CTRL_2"]]
CTRL_export <- as.data.frame(CTRL)

write.table(CTRL_export[,c(1,2,3,6,5,5)],paste0("~/KAS-METTL/",GROUP,"/06_macs2/broad/",GROUP,"_CTRL_common_nomodel_peaks.broadPeak"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

KO <- ol_KO$peaklist[["KO_1///KO_2"]]
KO_export <- as.data.frame(KO)
write.table(KO_export[,c(1,2,3,6,5,5)],paste0("~/KAS-METTL/",GROUP,"/06_macs2/broad/",GROUP,"_KO_common_nomodel_peaks.broadPeak"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

ol <- findOverlapsOfPeaks(CTRL, KO)
CTRL_unique <- ol$peaklist[["CTRL"]]
KO_unique <- ol$peaklist[["KO"]]


makeVennDiagram(ol, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 

library(GenomicFeatures)
txdb <- makeTxDbFromGFF('~/reference/annotation/hg38/gencode.v39.annotation.gtf')
annoData <- toGRanges(txdb, format='gene')
annoData[1:2]

overlaps_CTRL <- ol_CTRL$peaklist[["CTRL_1///CTRL_2"]]
binOverFeature(overlaps_CTRL, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(CTRL)")

overlaps_KO <- ol_KO$peaklist[["KO_1///KO_2"]]
binOverFeature(overlaps_KO, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(KO)")

overlaps <- ol$peaklist[["CTRL///KO"]]
binOverFeature(overlaps, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(KO)")

## check the genomic element distribution of the duplicates
## the genomic element distribution will indicates the 
## the correlation between duplicates.
peaks_CTRL <- GRangesList(rep1=CTRL_1,
                     rep2=CTRL_2)
genomicElementDistribution(peaks_CTRL, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=2000),
                           geneDownstream=c(upstream=0, downstream=2000))

peaks_KO <- GRangesList(rep1=KO_1,
                          rep2=KO_2)
genomicElementDistribution(peaks_KO, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=2000),
                           geneDownstream=c(upstream=0, downstream=2000))

peaks <- GRangesList(rep1=CTRL,
                        rep2=KO)
genomicElementDistribution(peaks, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=2000),
                           geneDownstream=c(upstream=0, downstream=2000))

## check the genomic element distribution for the overlaps
## the genomic element distribution will indicates the 
## the best methods for annotation.
## The percentages in the legend show the percentage of peaks in 
## each category.
out_CTRL <- genomicElementDistribution(overlaps_CTRL, 
                                  TxDb = txdb,
                                  promoterRegion=c(upstream=1000, downstream=100),
                                  geneDownstream=c(upstream=0, downstream=2000),
                                  promoterLevel=list(
                                    # from 5' -> 3', fixed precedence 3' -> 5'
                                    breaks = c( -1000, -500, 0, 100),
                                    labels = c( "upstream 0.5-1Kb", 
                                               "upstream <500b", "TSS - 100b"),
                                    colors = c("#FFE5CC", 
                                               "#FFAD65", "#FF8E32")))

out_KO <- genomicElementDistribution(overlaps_KO, 
                                  TxDb = txdb,
                                  promoterRegion=c(upstream=1000, downstream=100),
                                  geneDownstream=c(upstream=0, downstream=2000),
                                  promoterLevel=list(
                                    # from 5' -> 3', fixed precedence 3' -> 5'
                                    breaks = c( -1000, -500, 0, 100),
                                    labels = c( "upstream 0.5-1Kb", 
                                               "upstream <500b", "TSS - 100b"),
                                    colors = c("#FFE5CC",
                                               "#FFAD65", "#FF8E32")))
out_CTRL_unique <- genomicElementDistribution(CTRL_unique, 
                                  TxDb = txdb,
                                  promoterRegion=c(upstream=1000, downstream=100),
                                  geneDownstream=c(upstream=0, downstream=2000),
                                  promoterLevel=list(
                                    # from 5' -> 3', fixed precedence 3' -> 5'
                                    breaks = c(-1000, -500, 0, 100),
                                    labels = c( "upstream 0.5-1Kb", 
                                               "upstream <500b", "TSS - 100b"),
                                    colors = c("#FFE5CC", 
                                               "#FFAD65", "#FF8E32")))

out_KO_unique <- genomicElementDistribution(KO_unique, 
                                              TxDb = txdb,
                                              promoterRegion=c(upstream=1000, downstream=100),
                                              geneDownstream=c(upstream=0, downstream=2000),
                                              promoterLevel=list(
                                                # from 5' -> 3', fixed precedence 3' -> 5'
                                                breaks = c(-1000, -500, 0, 500),
                                                labels = c("upstream 0.5-1Kb", 
                                                           "upstream <500b", "TSS - 500b"),
                                                colors = c("#FFE5CC",
                                                           "#FFAD65", "#FF8E32")))

out <- genomicElementDistribution(overlaps, 
                                              TxDb = txdb,
                                              promoterRegion=c(upstream=2000, downstream=500),
                                              geneDownstream=c(upstream=0, downstream=2000),
                                              promoterLevel=list(
                                                # from 5' -> 3', fixed precedence 3' -> 5'
                                                breaks = c(-2000, -1000, -500, 0, 500),
                                                labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                                           "upstream <500b", "TSS - 500b"),
                                                colors = c("#FFE5CC", "#FFCA99", 
                                                           "#FFAD65", "#FF8E32")))
