rm(list=ls())
library(ChIPpeakAnno)

caN_bed1 <- "/disk1/home/user_08/custom/HepG2-caN_rep1_peaks.bed"
caN_bed2 <- "/disk1/home/user_08/custom/HepG2-caN_rep2_peaks.bed"

toN_bed1 <- "/disk1/home/user_08/custom/HepG2-toN_rep1_peaks.bed"
toN_bed2 <- "/disk1/home/user_08/custom/HepG2-toN_rep2_peaks.bed"

caN_1 <- ChIPpeakAnno::toGRanges(caN_bed1, format="BED", header=FALSE)
caN_2 <- ChIPpeakAnno::toGRanges(caN_bed2, format="BED", header=FALSE)
toN_1 <- ChIPpeakAnno::toGRanges(toN_bed1, format="BED", header=FALSE)
toN_2 <- ChIPpeakAnno::toGRanges(toN_bed2, format="BED", header=FALSE)


## must keep the class exactly same as gr1$score, i.e., numeric.
caN_1$score <- as.numeric(caN_1$score) 
caN_2$score <- as.numeric(caN_2$score) 
toN_1$score <- as.numeric(toN_1$score) 
toN_2$score <- as.numeric(toN_2$score) 

ol_caN <- findOverlapsOfPeaks(caN_1, caN_2)
## add metadata (mean of score) to the overlapping peaks
ol_caN <- addMetadata(ol_caN, colNames="score", FUN=base::mean) 
ol_caN$peaklist[["caN_1///caN_2"]]
makeVennDiagram(ol_caN, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 


ol_toN <- findOverlapsOfPeaks(toN_1, toN_2)
ol_toN <- addMetadata(ol_toN, colNames="score", FUN=base::mean) 
makeVennDiagram(ol_toN, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 

caN <- ol_caN$peaklist[["caN_1///caN_2"]]
caN_export <- as.data.frame(caN)
caN_export["width"] <- "."
caN_export$peakNames <- rownames(caN_export)


write.table(caN_export[,c(1,2,3,6,4,5)],"~/cato/peaklist/caN_common_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

toN <- ol_toN$peaklist[["toN_1///toN_2"]]
toN_export <- as.data.frame(toN)
toN_export["width"] <- "."
caN_export$peakNames <- rownames(caN_export)
write.table(toN_export[,c(1,2,3,6,4,5)],"~/cato/peaklist/toN_common_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

ol <- findOverlapsOfPeaks(caN, toN)

catoN_caN_unique <- ol$peaklist[["caN"]]
catoN_caN_unique_export <- as.data.frame(catoN_caN_unique)
catoN_caN_unique_export["width"] <- "."
write.table(catoN_caN_unique_export[,c(1,2,3,6,4,5)],"~/cato/peaklist/catoN_caN_unique_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

catoN_toN_unique <- ol$peaklist[["toN"]]
catoN_toN_unique_export <- as.data.frame(catoN_toN_unique)
catoN_toN_unique_export["width"] <- "."
write.table(catoN_toN_unique_export[,c(1,2,3,6,4,5)],"~/cato/peaklist/catoN_toN_unique_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

catoN_common <- ol$peaklist[["caN///toN"]]
catoN_common_export <- as.data.frame(catoN_common)
catoN_common_export["width"] <- "."
write.table(catoN_common_export[,c(1,2,3,6,4,5)],"~/cato/peaklist/catoN_common_peaks.bed",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")

makeVennDiagram(ol, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 

library(GenomicFeatures)
txdb <- makeTxDbFromGFF('~/reference/annotation/hg38/gencode.v39.annotation.gtf')
annoData <- toGRanges(txdb, format='gene')
annoData[1:2]

overlaps_caN <- ol_caN$peaklist[["caN_1///caN_2"]]
binOverFeature(overlaps_caN, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(caN)")

overlaps_toN <- ol_toN$peaklist[["toN_1///toN_2"]]
binOverFeature(overlaps_toN, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(toN)")

overlaps <- ol$peaklist[["caN///toN"]]
binOverFeature(overlaps, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(toN)")

## check the genomic element distribution of the duplicates
## the genomic element distribution will indicates the 
## the correlation between duplicates.
peaks_caN <- GRangesList(rep1=caN_1,
                     rep2=caN_2)
genomicElementDistribution(peaks_caN, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=2000),
                           geneDownstream=c(upstream=0, downstream=2000))

peaks_toN <- GRangesList(rep1=toN_1,
                          rep2=toN_2)
genomicElementDistribution(peaks_toN, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=2000),
                           geneDownstream=c(upstream=0, downstream=2000))

peaks <- GRangesList(rep1=caN,
                        rep2=toN)
genomicElementDistribution(peaks, 
                           TxDb = txdb,
                           promoterRegion=c(upstream=2000, downstream=2000),
                           geneDownstream=c(upstream=0, downstream=2000))

## check the genomic element distribution for the overlaps
## the genomic element distribution will indicates the 
## the best methods for annotation.
## The percentages in the legend show the percentage of peaks in 
## each category.
out_caN <- genomicElementDistribution(overlaps_caN, 
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

out_toN <- genomicElementDistribution(overlaps_toN, 
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
out_caN_unique <- genomicElementDistribution(caN_unique, 
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

out_toN_unique <- genomicElementDistribution(toN_unique, 
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
