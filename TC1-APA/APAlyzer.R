rm(list=ls())

#BiocManager::install('RJWANGbioinfo/APAlyzer')

library(APAlyzer)

suppressMessages(library("Rsamtools"))
flsall = c("/disk1/home/user_09/TC1_APA/04_bam_raw/TC1-AP_P0_rep1.bam","/disk1/home/user_09/TC1_APA/04_bam_raw/TC1-AP_P0_rep2.bam","/disk1/home/user_09/TC1_APA/04_bam_raw/TC1-AP_rP2n_rep1.bam","/disk1/home/user_09/TC1_APA/04_bam_raw/TC1-AP_rP2n_rep2.bam")
names(flsall) <- c("P0_rep1","P0_rep2","P10_rep1","P10_rep2")
flsall

library(repmis)

load("/disk1/home/user_09/reference/APAlyZer/mm10_REF.RData")

head(refUTRraw,2)
head(dfIPA,2)
head(dfLE,2)


PASREF=REF4PAS(refUTRraw,dfIPA,dfLE)
UTRdbraw=PASREF$UTRdbraw
dfIPAraw=PASREF$dfIPA
dfLEraw=PASREF$dfLE   

head(UTRdbraw,2)

DFUTRraw=PASEXP_3UTR(UTRdbraw, flsall, Strandtype="forward")

IPA_OUTraw=PASEXP_IPA(dfIPAraw, dfLEraw, flsall, Strandtype="forward", nts=40,minMQS = 10,SeqType = "SingleEnd")
head(IPA_OUTraw,2)

# Build the sample table with replicates
sampleTable1 = data.frame(samplename = c(names(flsall)),
                          condition = c(rep("P0",2),rep("P10",2)))
sampleTable1

# Analysis 3'UTR APA between KD and NT group using multi-repilicate design
test_3UTRmuti=APAdiff(sampleTable1,
                      DFUTRraw, 
                      conKET='P0',
                      trtKEY='P10',
                      PAS='3UTR',
                      CUTreads=1,
                      p_adjust_methods="fdr",
                      MultiTest='unpaired t-test')
head(test_3UTRmuti,2)
table(test_3UTRmuti$APAreg)

test_IPAmuti=APAdiff(sampleTable1,
                     IPA_OUTraw, 
                     conKET='P0',
                     trtKEY='P10',
                     PAS='IPA',
                     CUTreads=1,
                     p_adjust_methods="none",
                     MultiTest='unpaired t-test')
head(test_IPAmuti,2)
table(test_IPAmuti$APAreg)

APAVolcano(test_3UTRmuti, PAS='3UTR', Pcol = "pvalue", top=15, plot_title='3UTR APA',y_cutoff = 0.05,alpha = 0.3)

APAVolcano(test_IPAmuti, PAS='IPA', Pcol = "pvalue", top=5, plot_title='IPA APA',alpha = 0.5)

APABox(test_3UTRmuti, xlab = "APAreg", ylab = "RED", plot_title = NULL)

test_3UTRmuti$APA="3'UTR"
test_IPAmuti$APA="IPA"
dfplot=rbind(test_3UTRmuti[,c('RED','APA')],test_IPAmuti[,c('RED','APA')])

library(ggplot2)
###violin
ggplot(dfplot, aes(x = APA, y = RED)) + 
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.2)+ theme_bw() + 
  geom_hline(yintercept=0, linetype="dashed", color = "red")

###CDF
ggplot(dfplot, aes( x = RED, color = APA)) + 
  stat_ecdf(geom = "step") +
  ylab("cumulative fraction")+ 
  geom_vline(xintercept=0, linetype="dashed", color = "gray")+ theme_bw() + 
  geom_hline(yintercept=0.5, linetype="dashed", color = "gray")

suppressMessages(library("GenomicFeatures"))
suppressMessages(library("org.Mm.eg.db"))
extpath = system.file("extdata", "mm9.chr19.refGene.R.DB", package="APAlyzer")
txdb=loadDb(extpath, packageName='GenomicFeatures')
IDDB = org.Mm.eg.db
CDSdbraw=REFCDS(txdb,IDDB)

DFGENEraw=GENEXP_CDS(CDSdbraw, flsall, Strandtype="NONE")

## Extract 3 prime most alignment of a paired-end 
## bam file and saved into a new bam file
Bamfile='/path/to/inputdir/input.bam'
Outdir='/path/to/outdir/'  
StrandType="forward-reverse"    ## "forward-reverse",  or "reverse-forward" or "NONE"   
ThreeMostPairBam (BamfilePath=Bamfile, 
                  OutDirPath=Outdir, 
                  StrandType='forward-reverse')

ThreemostBamfile="test.3most.bam"
IPA_OUTraw=PASEXP_IPA(dfIPA, dfLE, ThreemostBamfile, Strandtype="forward", SeqType="ThreeMostPairEnd")