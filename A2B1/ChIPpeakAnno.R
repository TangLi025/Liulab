rm(list=ls())
library(ChIPpeakAnno)

Lysate_bed1 <- "/disk1/home/user_09/LinLong/08_bed_filtered/dedup/Lysate_rep1_peaks.bed"
Lysate_bed2 <- "/disk1/home/user_09/LinLong/08_bed_filtered/dedup/Lysate_rep2_peaks.bed"

Result_bed1 <- "/disk1/home/user_09/LinLong/08_bed_filtered/dedup/Result_rep1_peaks.bed"
Result_bed2 <- "/disk1/home/user_09/LinLong/08_bed_filtered/dedup/Result_rep2_peaks.bed"

Lysate_1 <- ChIPpeakAnno::toGRanges(Lysate_bed1, format="BED", header=FALSE)
Lysate_2 <- ChIPpeakAnno::toGRanges(Lysate_bed2, format="BED", header=FALSE)
Result_1 <- ChIPpeakAnno::toGRanges(Result_bed1, format="BED", header=FALSE)
Result_2 <- ChIPpeakAnno::toGRanges(Result_bed2, format="BED", header=FALSE)


## must keep the class exactly same as gr1$score, i.e., numeric.
Lysate_1$score <- as.numeric(Lysate_1$score) 
Lysate_2$score <- as.numeric(Lysate_2$score) 
Result_1$score <- as.numeric(Result_1$score) 
Result_2$score <- as.numeric(Result_2$score) 

ol_Lysate <- findOverlapsOfPeaks(Lysate_1, Lysate_2)
## add metadata (mean of score) to the overlapping peaks
ol_Lysate <- addMetadata(ol_Lysate, colNames="score", FUN=base::mean) 
ol_Lysate$peaklist[["Lysate_1///Lysate_2"]]
makeVennDiagram(ol_Lysate, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 


ol_Result <- findOverlapsOfPeaks(Result_1, Result_2)
ol_Result <- addMetadata(ol_Result, colNames="score", FUN=base::mean) 
makeVennDiagram(ol_Result, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 



Lysate <- ol_Lysate$peaklist[["Lysate_1///Lysate_2"]]
Lysate_export <- as.data.frame(Lysate)

write.table(Lysate_export[,c(1,2,3,6,5,5)],"~/LinLong/08_bed_filtered/dedup/Lysate_common_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

Result <- ol_Result$peaklist[["Result_1///Result_2"]]
Result_export <- as.data.frame(Result)
write.table(Result_export[,c(1,2,3,6,5,5)],"~/LinLong/08_bed_filtered/dedup/Result_common_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

ol <- findOverlapsOfPeaks(Lysate, Result)
Lysate_unique <- ol$peaklist[["Lysate"]]


makeVennDiagram(ol, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 

library(GenomicFeatures)
txdb <- makeTxDbFromGFF('~/reference/annotation/mm19/gencode.vM28.annotation.protein_coding.chr.gtf')
annoData <- toGRanges(txdb, format='gene')
annoData[1:2]

overlaps_Lysate <- ol_Lysate$peaklist[["Lysate_1///Lysate_2"]]
binOverFeature(overlaps_Lysate, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(Lysate)")

overlaps_Result <- ol_Result$peaklist[["Result_1///Result_2"]]
binOverFeature(overlaps_Result, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(Result)")

overlaps <- ol$peaklist[["Lysate///Result"]]
binOverFeature(overlaps, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(Result)")

## check the genomic element distribution of the duplicates
## the genomic element distribution will indicates the 
## the correlation between duplicates.
peaks_Lysate <- GRangesList(rep1=Lysate_1,
                     rep2=Lysate_2)
genomicElementDistribution(peaks_Lysate, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=2000),
                           geneDownstream=c(upstream=0, downstream=2000))

peaks_Result <- GRangesList(rep1=Result_1,
                          rep2=Result_2)
genomicElementDistribution(peaks_Result, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=2000),
                           geneDownstream=c(upstream=0, downstream=2000))

peaks <- GRangesList(rep1=Lysate,
                        rep2=Result)
genomicElementDistribution(peaks, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=2000),
                           geneDownstream=c(upstream=0, downstream=2000))


## check the genomic element distribution for the overlaps
## the genomic element distribution will indicates the 
## the best methods for annotation.
## The percentages in the legend show the percentage of peaks in 
## each category.
out_Lysate <- genomicElementDistribution(overlaps_Lysate, 
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

out_Result <- genomicElementDistribution(overlaps_Result, 
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
out_Lysate_unique <- genomicElementDistribution(Lysate_unique, 
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
