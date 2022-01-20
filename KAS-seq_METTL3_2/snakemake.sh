#!/bin/bash
# creatation: 2022-1-13
# Author: Candy Lee tangli025@pku.edu.cn)

# Stop on error
set -e

THREAD=50

DIRS=( \
      "/disk1/home/user_09/KAS/" \
      )

for DIR in ${DIRS[@]}
do
  echo ${DIR}
  #echo "snakemake -s fastqc.py -c ${THREAD} -d ${DIR}"
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s fastqc.py -n -d ${DIR}
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s fastqc.py -c ${THREAD} -d ${DIR}
  #echo "snakemake -s bowtie2_mapping.py -c ${THREAD} -d ${DIR}"
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s bowtie2_mapping.py -n -d ${DIR}
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s bowtie2_mapping.py -c ${THREAD} -d ${DIR}
  echo "snakemake -s bam2bg.py -c ${THREAD} -d ${DIR}"
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s bam2bg.py -n -d ${DIR}
  /disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s bam2bg.py -c ${THREAD} -d ${DIR}
  echo "snakemake -s macs2_callpeak.py -c ${THREAD} -d ${DIR}"
  /disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s macs2_callpeak.py -c ${THREAD} -d ${DIR}
  echo "snakemake -s plotFingerprint.py -c ${THREAD} -d ${DIR}"
  /disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s plotFingerprint.py -c ${THREAD} -d ${DIR}
done