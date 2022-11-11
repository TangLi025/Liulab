#!/bin/bash

SAMPLE=("Lysate" "Result")
TREATMENT=("input" "IP")
REP=("rep1" "rep2")

for sample in ${SAMPLE[@]}
do
  for treatment in ${TREATMENT[@]}
  do
    for rep in ${REP[@]}
    do
      echo ~/LinLong/04_bam_raw/${sample}.${rep}.${treatment}.bam
      #echo -n "${group}_${sample}_${treatment}_${rep}_ext.bg " >> ~/KAS-METTL/${group}/05_bedtools/bedGraph/KAS-seq_file.txt
      samtools idxstats ~/LinLong/04_bam_raw/${sample}.${rep}.${treatment}.bam | awk '{sum+=$3} END {print sum"\n"}' >> ~/LinLong/04_bam_raw/bam_summary.txt
     done
  done
done