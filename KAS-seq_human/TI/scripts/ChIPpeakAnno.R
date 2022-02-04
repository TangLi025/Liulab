library(ChIPpeakAnno)

bed1 <- "~/KAS/TI_KAS-seq/06_macs2/bed_merge/DMSO_common_peaks.broadPeak"
bed2 <- "~/KAS/TI_KAS-seq/06_macs2/bed_merge/DRB_common_peaks.broadPeak"
bed3 <- "~/KAS/TI_KAS-seq/06_macs2/bed_merge/TRIP_common_peaks.broadPeak"

Native <- ChIPpeakAnno::toGRanges(bed1, format="BED", header=FALSE)
DRB <- ChIPpeakAnno::toGRanges(bed2, format="BED", header=FALSE)
Triptolide <- ChIPpeakAnno::toGRanges(bed3, format="BED", header=FALSE)



## must keep the class exactly same as gr1$score, i.e., numeric.
Native$score <- as.numeric(Native$score) 
DRB$score <- as.numeric(DRB$score) 
Triptolide$score <- as.numeric(Triptolide$score) 
ol <- findOverlapsOfPeaks(Native, DRB,Triptolide)
## add metadata (mean of score) to the overlapping peaks
ol <- addMetadata(ol, colNames="score", FUN=mean) 
ol$peaklist[["rep1///rep2"]][1:2]

makeVennDiagram(ol, fill=c("#a1d8b1", "#edfcc2","#f88aaf"), # circle fill color
                col=c("#D55E00", "#0072B2","#f88aaf"), #circle border color
                cat.col=c("#D55E00", "#0072B2","#f88aaf"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 

library(EnsDb.Hsapiens.v75)
annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")
annoData[1:2]