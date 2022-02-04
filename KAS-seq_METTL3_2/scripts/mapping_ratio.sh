#!/bin/bash

SAMPLE=("Ctrl" "KO")
TREATMENT=("input" "IP")
REP=("rep1" "rep2")

GROUP=("KAS-seq_ALKBH3" "KAS-seq_ALKBH5")

for group in ${GROUP[@]}
do
for sample in ${SAMPLE[@]}
do
  for treatment in ${TREATMENT[@]}
  do
    for rep in ${REP[@]}
    do
      echo ~/KAS-METTL/${group}/04_bam_rmdup/${group}_${sample}_${treatment}_${rep}.bam
      echo -n "${group}_${sample}_${treatment}_${rep}_ext.bg " >> ~/KAS-METTL/${group}/05_bedtools/bedGraph/KAS-seq_file.txt
      samtools idxstats ~/KAS-METTL/${group}/04_bam_rmdup/${group}_${sample}_${treatment}_${rep}.bam | awk '{sum+=$3} END {print sum"\n"}' >> ~/KAS-METTL/${group}/05_bedtools/bedGraph/bam_summary.txt
      #samtools view -F 0x04 ~/KAS-METTL/${group}/04_bam_rmdup/${group}_${sample}_${treatment}_${rep}.bam | wc -l >> ~/KAS-METTL/${group}/05_bedtools/bedGraph/bam_summary.txt
    done
  done
done
done