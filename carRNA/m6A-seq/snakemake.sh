#!/bin/bash

THREAD=50

DIRS=( \
      "/disk1/home/user_09/LinLong" \
      )

for DIR in ${DIRS[@]}
do
  /disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s fastqc.py -c ${THREAD} -d ${DIR}
  /disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s bowtie2_mapping.py -c ${THREAD} -d ${DIR}
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s bedtools.py -c ${THREAD} -d ${DIR}
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s deeptools.py -c ${THREAD} -d ${DIR} 
done