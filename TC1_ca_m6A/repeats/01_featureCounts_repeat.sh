#!/bin/bash

cd /disk/user_08/Data/TC1-planB/05_bam_change_index_hisat2/01_bam_sorted/11_bam_merge/

/disk1/home/user_09/anaconda3/envs/subread/bin/featureCounts \
  -a /disk/user_09/reference/annotation/mm19/repeats/mm19_repeats_family.saf \
  -o 06_repeats/repeats_family_counts.tab \
  -F SAF --fracOverlap 0.5 -M -s 0 -p -B -P -T 50 p0*.bam p5*.bam p10*.bam rp2*.bam
  