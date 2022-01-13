#!/bin/bash

SAMPLE=("CTRL" "KO")
TREATMENT=("input" "IP")
REP=("rep1" "rep2")



for sample in ${SAMPLE[@]}
do
  for treatment in ${TREATMENT[@]}
  do
    for rep in ${REP[@]}
    do
      #samtools view -c ~/KAS-METTL/METTL3_2/04_bam_rmdup/KAS-seq_METTL3_2_${sample}_${treatment}_${rep}.bam >> ~/KAS-METTL/METTL3_2/04_bam_rmdup/bam_summary.txt
      samtools view -F 0x04 ~/KAS-METTL/METTL3_2/04_bam_rmdup/KAS-seq_METTL3_2_${sample}_${treatment}_${rep}.bam | wc -l >> ~/KAS-METTL/METTL3_2/04_bam_rmdup/bam_summary.txt
    done
  done
done