#! /usr/lib/R/bin/Rscript --vanilla

library(APAlyzer)

suppressMessages(library("Rsamtools"))
flsall = c("/disk1/home/user_09/TC1_m6A/04_bam_raw/TC1_P0_input_rep1.bam","/disk1/home/user_09/TC1_m6A/04_bam_raw/TC1_P0_input_rep2.bam","/disk1/home/user_09/TC1_m6A/04_bam_raw/TC1_P5_input_rep1.bam","/disk1/home/user_09/TC1_m6A/04_bam_raw/TC1_P5_input_rep2.bam","/disk1/home/user_09/TC1_m6A/04_bam_raw/TC1_P10_input_rep1.bam","/disk1/home/user_09/TC1_m6A/04_bam_raw/TC1_P10_input_rep2.bam","/disk1/home/user_09/TC1_m6A/04_bam_raw/TC1_rP2_input_rep1.bam","/disk1/home/user_09/TC1_m6A/04_bam_raw/TC1_rP2_input_rep2.bam")
names(flsall) <- c("P0_rep1","P0_rep2","P5_rep1","P5_rep2","P10_rep1","P10_rep2","rP2_rep1","rP2_rep2")

library(repmis)

for (i in 1:8){
  Bamfile <- flsall[i]
  Outdir='/disk1/home/user_09/TC1_m6A/04_bam_raw_3most/'  
  StrandType="reverse-forward"    ## "forward-reverse",  or "reverse-forward" or "NONE"   
  ThreeMostPairBam (BamfilePath=Bamfile, 
                    OutDirPath=Outdir, 
                    StrandType=StrandType)
}