#!/bin/bash

SAMPLE=("Lysate" "Result")
TREATMENT=("input" "IP")
REP=("rep1" "rep2")

DUP=("raw" "dedup")

for dup in ${DUP[@]}
do
for sample in ${SAMPLE[@]}
do
  for treatment in ${TREATMENT[@]}
  do
    for rep in ${REP[@]}
    do
      echo ~/LinLong/04_bam_${dup}/${sample}.${rep}.${treatment}.bam
      echo -n "${sample}.${rep}.${treatment}_ext.bg " >> ~/LinLong/05_bedtools/${dup}/bedGraph/MeRIP-seq_file.txt
      samtools idxstats ~/LinLong/04_bam_${dup}/${sample}.${rep}.${treatment}.bam | awk '{sum+=$3} END {print sum"\n"}' >> ~/LinLong/05_bedtools/${dup}/bedGraph/bam_summary.txt
    done
  done
done
done