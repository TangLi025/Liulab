rm(list=ls())
library(ChIPpeakAnno)

CTRL_bed1 <- paste0("~/carRNA_science/m6A-seq_narrowPeak/GSM3912804_Mettl3_KO2_r1_peaks.narrowPeak")
CTRL_bed2 <- paste0("~/carRNA_science/m6A-seq_narrowPeak/GSM3912805_Mettl3_KO2_r2_peaks.narrowPeak")

CTRL_1 <- ChIPpeakAnno::toGRanges(CTRL_bed1, format="BED", header=FALSE)
CTRL_2 <- ChIPpeakAnno::toGRanges(CTRL_bed2, format="BED", header=FALSE)

## must keep the class exactly same as gr1$score, i.e., numeric.
CTRL_1$score <- as.numeric(CTRL_1$score) 
CTRL_2$score <- as.numeric(CTRL_2$score) 

ol_CTRL <- findOverlapsOfPeaks(CTRL_1, CTRL_2)
## add metadata (mean of score) to the overlapping peaks
ol_CTRL <- addMetadata(ol_CTRL, colNames="score", FUN=mean) 
ol_CTRL$peaklist[["CTRL_1///CTRL_2"]]
makeVennDiagram(ol_CTRL, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 

CTRL <- ol_CTRL$peaklist[["CTRL_1///CTRL_2"]]
CTRL_export <- as.data.frame(CTRL)

write.table(CTRL_export[,c(1,2,3)],"~/carRNA_science/m6A-seq_narrowPeak/Mettl3_KO2_commen_peaks.narrowPeak",quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")


KO <- ol_KO$peaklist[["KO_1///KO_2"]]

ol <- findOverlapsOfPeaks(CTRL, KO)


makeVennDiagram(ol, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 

library(GenomicFeatures)
txdb <- makeTxDbFromGFF('~/reference/annotation/mm10/Mus_musculus.GRCm38.102.chr.gtf')
annoData <- toGRanges(txdb, format='gene')
annoData[1:2]

overlaps_CTRL <- ol_CTRL$peaklist[["CTRL_1///CTRL_2"]]
binOverFeature(overlaps_CTRL, annotationData=annoData,
               radius=5000, nbins=50, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS(CTRL)")

## check the genomic element distribution of the duplicates
## the genomic element distribution will indicates the 
## the correlation between duplicates.
peaks_CTRL <- GRangesList(rep1=CTRL_1,
                     rep2=CTRL_2)
genomicElementDistribution(peaks_CTRL, 
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
                                  promoterRegion=c(upstream=2000, downstream=500),
                                  geneDownstream=c(upstream=0, downstream=5000),
                                  promoterLevel=list(
                                    # from 5' -> 3', fixed precedence 3' -> 5'
                                    breaks = c(-2000, -1000, -500, 0, 500),
                                    labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                               "upstream <500b", "TSS - 500b"),
                                    colors = c("#FFE5CC", "#FFCA99", 
                                               "#FFAD65", "#FF8E32")))




overlaps.anno <- annotatePeakInBatch(overlaps, 
                                     AnnotationData=annoData, 
                                     output="nearestBiDirectionalPromoters",
                                     bindingRegion=c(-2000, 500))
library(org.Hs.eg.db)
overlaps.anno <- addGeneIDs(overlaps.anno,
                            "org.Hs.eg.db",
                            IDs2Add = "entrez_id")
head(overlaps.anno)










