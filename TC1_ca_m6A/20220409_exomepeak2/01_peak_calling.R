library(exomePeak2)

GENE_ANNO_GTF = "~/reference/annotation/mm10/gencode.vM25.annotation.gtf"

f1 = "/disk/user_09/Data/03_TC1_caRNA/11_bam_merge/raw/p0_ip_rep1.sorted.bam"
f2 = "/disk/user_09/Data/03_TC1_caRNA/11_bam_merge/raw/p0_ip_rep2.sorted.bam"

IP_BAM = c(f1,f2)

f1 = "/disk/user_09/Data/03_TC1_caRNA/11_bam_merge/raw/p0_input_rep1.sorted.bam"
f2 = "/disk/user_09/Data/03_TC1_caRNA/11_bam_merge/raw/p0_input_rep2.sorted.bam"

INPUT_BAM = c(f1,f2)

f1 = "/disk/user_09/Data/03_TC1_caRNA/11_bam_merge/raw/p10_ip_rep1.sorted.bam"
f2 = "/disk/user_09/Data/03_TC1_caRNA/11_bam_merge/raw/p10_ip_rep2.sorted.bam"

TREATED_IP_BAM = c(f1,f2)

f1 = "/disk/user_09/Data/03_TC1_caRNA/11_bam_merge/raw/p10_input_rep1.sorted.bam"
f2 = "/disk/user_09/Data/03_TC1_caRNA/11_bam_merge/raw/p10_input_rep2.sorted.bam"

TREATED_INPUT_BAM = c(f1,f2)

# MeRIP_Seq_Alignment <- scanMeripBAM(
#   bam_ip = IP_BAM,
#   bam_input = INPUT_BAM,
#   paired_end = TRUE
# )
# 
# library(BSgenome)
# library(GenomicFeatures)
# SummarizedExomePeaks <- exomePeakCalling(merip_bams = MeRIP_Seq_Alignment,
#                                          gff_dir = GENE_ANNO_GTF,
#                                          fragment_length = 150,
#                                          genome) 
exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           bam_treated_input = TREATED_INPUT_BAM,
           bam_treated_ip = TREATED_IP_BAM,
           gff_dir = GENE_ANNO_GTF,
           fragment_length = 150,
           parallel = 20,
           save_plot_analysis = TRUE,
           peak_calling_mode = "full_tx",
           genome = 'mm10',
           paired_end = T)
