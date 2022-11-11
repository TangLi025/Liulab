#!/bin/bash

THREAD=50

DIRS=( \
      "/disk1/home/user_09/KAS-seq_human/H3K27ac" \
      "/disk1/home/user_09/KAS-seq_human/H3K4me1"
      )

for DIR in ${DIRS[@]}
do
  echo ${DIR}
  echo "snakemake -s fastqc.py -c ${THREAD} -d ${DIR}"
  /disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s fastqc.py -c ${THREAD} -d ${DIR}
  echo "snakemake -s bowtie2_mapping.py -c ${THREAD} -d ${DIR}"
  /disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s bowtie2_mapping.py -c ${THREAD} -d ${DIR}
  #echo "snakemake -s bedtools.py -c ${THREAD} -d ${DIR}"
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s bedtools.py -c ${THREAD} -d ${DIR}
  #echo "snakemake -s deeptools.py -c ${THREAD} -d ${DIR}"
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s deeptools.py -c ${THREAD} -d ${DIR}
done