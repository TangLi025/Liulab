#!/bin/bash

SAMPLE=("Control" "KDMETTL3")
TREATMENT=("input" "IP")
REP=("rep1" "rep2")

GROUP=("A375")

for group in ${GROUP[@]}
do
for sample in ${SAMPLE[@]}
do
  for treatment in ${TREATMENT[@]}
  do
    for rep in ${REP[@]}
    do
      #samtools view -c ~/KAS-METTL/METTL3_2/04_bam_rmdup/KAS-seq_METTL3_2_${sample}_${treatment}_${rep}.bam >> ~/KAS-METTL/METTL3_2/04_bam_rmdup/bam_summary.txt
      samtools view -F 0x04 ~/KAS-METTL/${group}/04_bam_rmdup/${group}_${sample}_${treatment}_${rep}.bam | wc -l >> ~/KAS-METTL/${group}/05_bedtools/bedGraph/bam_summary.txt
    done
  done
done
done