#!/bin/bash

sort -k1,1 -k2,2n /disk/user_09/Data/03_TC1_caRNA/12_dapars/02_pAPA_loci/APA_loci.bed > \
    /disk/user_09/Data/03_TC1_caRNA/12_dapars/02_pAPA_loci/APA_loci.sorted.bed

for i in {p0,p5,p10,rp2}
do
sort -k1,1 -k2,2n /disk/user_09/Data/03_TC1_caRNA/10_bed_merge/01_peak_center/${i}_peak_center.bed > \
    /disk/user_09/Data/03_TC1_caRNA/10_bed_merge/01_peak_center/${i}_peak_center.sorted.bed
bedtools closest -a /disk/user_09/Data/03_TC1_caRNA/12_dapars/02_pAPA_loci/APA_loci.sorted.bed \
    -b /disk/user_09/Data/03_TC1_caRNA/10_bed_merge/01_peak_center/${i}_peak_center.sorted.bed \
    -s -D a -k 3 > /disk/user_09/Data/03_TC1_caRNA/10_bed_merge/02_closest/${i}_pAPA_closest3.bed

done  