#!/bin/bash

cd /disk/user_09/user_08_TC1/10_mRNA/05_bam_hisat2/01_bam_sorted/11_bam_merge

/disk1/home/user_09/anaconda3/envs/subread/bin/featureCounts \
  -a /disk/user_08/Data/annotation/mm39/gencode.vM28.annotation.gtf \
  -o 05_gene_expression_tl/CDS_express_counts.tab \
  -F GTF -t 'CDS' -g "transcript_id" --extraAttributes "gene_id" \
  -O --fracOverlap 0.5 -s 2 -p -B -P -d 40 -D 2000 -C \
  -T 50 p0*.bam p5*.bam p10*.bam rp2*.bam
  