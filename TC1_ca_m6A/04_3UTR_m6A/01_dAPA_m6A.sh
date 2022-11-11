#!/bin/bash
for i in {p0,p5,p10,rp2}
do
bedtools intersect -a /disk/user_09/Data/03_TC1_caRNA/12_dapars/dAPA_region.bed \
    -b /disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/10_bed_merge/00_common_peaks/${i}_rep1_rep2_common_peaks.bed \
    -wa -c -f 0.5 -F 0.9 -e > /disk/user_09/Data/03_TC1_caRNA/12_dapars/01_PDUI_dAPA_m6A_peak/PDUI_dAPA_peak_${i}.bed
done