#!/bin/bash

for i in {p0,p5,p10,rp2}; do bedtools intersect -a ${i}_rep1_rep2_common_peaks.bed -b /disk1/home/user_09/reference/annotation/mm19/dapars/mm39_gencode_dapars_longest_protein.bed -wa -e -f 0.9 -F 0.9 > peak_3UTR/${i}_common_3UTR_f0.1.bed; done

for i in {p0,p5,p10,rp2}; do /disk1/home/user_09/anaconda3/envs/LinLong/bin/findMotifsGenome.pl ${i}_common_3UTR_e_f0.9_F0.9.bed /disk1/home/user_09/reference/genome/mm/GRCm39.genome.fa 02_motif_0.9/${i} -rna -p 5 -len 5,6,7 > ${i}_common_3UTR_e_f0.9_F0.9_motif.log 2>&1 & done