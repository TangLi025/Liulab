library(ChIPpeakAnno)

CTRL_bed1 <- "~/KAS-METTL/METTL3_2/06_macs2/broad/KAS-seq_METTL3_2_CTRL_rep1_peaks.broadPeak"
CTRL_bed2 <- "~/KAS-METTL/METTL3_2/06_macs2/broad/KAS-seq_METTL3_2_CTRL_rep2_peaks.broadPeak"

KO_bed1 <- "~/KAS-METTL/METTL3_2/06_macs2/broad/KAS-seq_METTL3_2_KO_rep1_peaks.broadPeak"
KO_bed2 <- "~/KAS-METTL/METTL3_2/06_macs2/broad/KAS-seq_METTL3_2_KO_rep2_peaks.broadPeak"

CTRL_1 <- ChIPpeakAnno::toGRanges(CTRL_bed1, format="BED", header=FALSE)
CTRL_2 <- ChIPpeakAnno::toGRanges(CTRL_bed2, format="BED", header=FALSE)
KO_1 <- ChIPpeakAnno::toGRanges(KO_bed1, format="BED", header=FALSE)
KO_2 <- ChIPpeakAnno::toGRanges(KO_bed2, format="BED", header=FALSE)


## must keep the class exactly same as gr1$score, i.e., numeric.
CTRL$score <- as.numeric(CTRL$score) 
KO$score <- as.numeric(KO$score) 
ol_CTRL <- findOverlapsOfPeaks(CTRL_1, CTRL_2)


ol_KO <- findOverlapsOfPeaks(KO_1, KO_2)
## add metadata (mean of score) to the overlapping peaks
ol_CTRL <- addMetadata(ol_CTRL, colNames="score", FUN=mean) 
ol_CTRL$peaklist[["CTRL_1///CTRL_2"]]
CTRL <- ol_CTRL$peaklist[["CTRL_1///CTRL_2"]]
KO <- ol_KO$peaklist[["KO_1///KO_2"]]

ol <- findOverlapsOfPeaks(CTRL, KO)


makeVennDiagram(ol, fill=c("#a1d8b1", "#edfcc2"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2"),# label color, keep same as circle border color
                cex=3,
                cat.cex=2) 



library(EnsDb.Hsapiens.v75)
annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")
annoData[1:2]