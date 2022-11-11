rm(list=ls())
library(ChIPpeakAnno)

P0_bed1 <- "/disk1/home/user_09/TC1_m6A/08_bed_filtered/dedup/TC1_P0_rep1_peaks.bed"
P0_bed2 <- "/disk1/home/user_09/TC1_m6A/08_bed_filtered/dedup/TC1_P0_rep2_peaks.bed"

rP2_bed1 <- "/disk1/home/user_09/TC1_m6A/08_bed_filtered/dedup/TC1_rP2_rep1_peaks.bed"
rP2_bed2 <- "/disk1/home/user_09/TC1_m6A/08_bed_filtered/dedup/TC1_rP2_rep2_peaks.bed"

P0_1 <- ChIPpeakAnno::toGRanges(P0_bed1, format="BED", header=FALSE)
P0_2 <- ChIPpeakAnno::toGRanges(P0_bed2, format="BED", header=FALSE)
rP2_1 <- ChIPpeakAnno::toGRanges(rP2_bed1, format="BED", header=FALSE)
rP2_2 <- ChIPpeakAnno::toGRanges(rP2_bed2, format="BED", header=FALSE)


## must keep the class exactly same as gr1$score, i.e., numeric.
P0_1$score <- as.numeric(P0_1$score) 
P0_2$score <- as.numeric(P0_2$score) 
rP2_1$score <- as.numeric(rP2_1$score) 
rP2_2$score <- as.numeric(rP2_2$score) 

ol_P0 <- findOverlapsOfPeaks(P0_1, P0_2)
## add metadata (mean of score) to the overlapping peaks
ol_P0 <- addMetadata(ol_P0, colNames="score", FUN=base::mean) 
ol_P0$peaklist[["P0_1///P0_2"]]
makeVennDiagram(ol_P0, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 


ol_rP2 <- findOverlapsOfPeaks(rP2_1, rP2_2)
ol_rP2 <- addMetadata(ol_rP2, colNames="score", FUN=base::mean) 
makeVennDiagram(ol_rP2, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 



P0 <- ol_P0$peaklist[["P0_1///P0_2"]]
P0_export <- as.data.frame(P0)

write.table(P0_export[,c(1,2,3,6,5,5)],"~/LinLong/08_bed_filtered/dedup/P0_common_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

rP2 <- ol_rP2$peaklist[["rP2_1///rP2_2"]]
rP2_export <- as.data.frame(rP2)
write.table(rP2_export[,c(1,2,3,6,5,5)],"~/LinLong/08_bed_filtered/dedup/rP2_common_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

ol <- findOverlapsOfPeaks(P0, rP2)
P0_unique <- ol$peaklist[["P0"]]


makeVennDiagram(ol, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 

library(GenomicFeatures)
txdb <- makeTxDbFromGFF('~/reference/annotation/mm19/gencode.vM28.annotation.protein_coding.chr.gtf')
annoData <- toGRanges(txdb, format='gene')
annoData[1:2]

overlaps_P0 <- ol_P0$peaklist[["P0_1///P0_2"]]
binOverFeature(overlaps_P0, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(P0)")

overlaps_rP2 <- ol_rP2$peaklist[["rP2_1///rP2_2"]]
binOverFeature(overlaps_rP2, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(rP2)")

overlaps <- ol$peaklist[["P0///rP2"]]
binOverFeature(overlaps, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(rP2)")

## check the genomic element distribution of the duplicates
## the genomic element distribution will indicates the 
## the correlation between duplicates.
peaks_P0 <- GRangesList(rep1=P0_1,
                            rep2=P0_2)
genomicElementDistribution(peaks_P0, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=2000),
                           geneDownstream=c(upstream=0, downstream=2000))

peaks_rP2 <- GRangesList(rep1=rP2_1,
                            rep2=rP2_2)
genomicElementDistribution(peaks_rP2, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=2000),
                           geneDownstream=c(upstream=0, downstream=2000))

peaks <- GRangesList(rep1=P0,
                     rep2=rP2)
genomicElementDistribution(peaks, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=2000),
                           geneDownstream=c(upstream=0, downstream=2000))


## check the genomic element distribution for the overlaps
## the genomic element distribution will indicates the 
## the best methods for annotation.
## The percentages in the legend show the percentage of peaks in 
## each category.
out_P0 <- genomicElementDistribution(overlaps_P0, 
                                         TxDb = txdb,
                                         promoterRegion=c(upstream=2000, downstream=500),
                                         geneDownstream=c(upstream=0, downstream=5000),
                                         promoterLevel=list(
                                           # from 5' -> 3', fixed precedence 3' -> 5'
                                           breaks = c(-2000, -1000, -500, 0, 500),
                                           labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                                      "upstream <500b", "TSS - 500b"),
                                           colors = c("#FFE5CC", "#FFCA99", 
                                                      "#FFAD65", "#FF8E32")))

out_rP2 <- genomicElementDistribution(overlaps_rP2, 
                                         TxDb = txdb,
                                         promoterRegion=c(upstream=2000, downstream=500),
                                         geneDownstream=c(upstream=0, downstream=5000),
                                         promoterLevel=list(
                                           # from 5' -> 3', fixed precedence 3' -> 5'
                                           breaks = c(-2000, -1000, -500, 0, 500),
                                           labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                                      "upstream <500b", "TSS - 500b"),
                                           colors = c("#FFE5CC", "#FFCA99", 
                                                      "#FFAD65", "#FF8E32")))
out_P0_unique <- genomicElementDistribution(P0_unique, 
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
