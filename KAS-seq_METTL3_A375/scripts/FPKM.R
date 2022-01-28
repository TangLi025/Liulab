hg_anno <- rtracklayer::import('~/reference/annotation/hg19/gencode.v19.annotation.gene_coding.gtf')
hg_anno <- rtracklayer::as.data.frame(hg_anno)
hg_transcript <- hg_anno[hg_anno["type"] == "transcript",]


colnames_FPKM <- c("chrom","st","end","accession","mRNA_size","gene_strand","Frag_count","FPM","FPKM")
FPKM_1 <- read.table("~/KAS/RNA-seq/05_FPKM_counts/HEK_rep1.FPKM.xls",col.names = colnames_FPKM)

hg_anno_in <- hg_transcript[which(hg_transcript$transcript_id %in% FPKM_1$accession),]

FPKM_gene_coding <- FPKM_1[which(FPKM_1$accession %in% hg_transcript$transcript_id),]
FPKM_gene_coding_dedup <- FPKM_gene_coding[FPKM_gene_coding$chrom == "chrY",]

FPKM_dup <- FPKM_gene_coding[duplicated(FPKM_gene_coding$accession),]
FPKM_dup2 <- FPKM_gene_coding[which(FPKM_gene_coding$accession %in% FPKM_dup$accession),]

hg_transcript_dup <- hg_transcript[which(hg_transcript$transcript_id %in% FPKM_dup$accession),]
hg_transcript_Y <- hg_transcript[which(hg_transcript$seqnames == "chrY"),]
