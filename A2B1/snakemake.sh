#!/bin/bash
# creatation: 2022-1-13
# Author: Candy Lee tangli025@pku.edu.cn)

# Stop on error
set -e

THREAD=50

DIRS=( \
      "/disk1/home/user_09/LinLong" \
      )

for DIR in ${DIRS[@]}
do
  echo ${DIR}
  echo "snakemake -s cutadapter.py -c ${THREAD} -d ${DIR}"
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s cutadapter.py -n -d ${DIR}
  /disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s cutadapter.py -c ${THREAD} -d ${DIR}
  #echo "snakemake -s bam2bg.py -c ${THREAD} -d ${DIR}"
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s bam2bg.py -n -d ${DIR}
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s bam2bg.py -c ${THREAD} -d ${DIR}
  echo "snakemake -s macs2_callpeak.py -c ${THREAD} -d ${DIR}"
  /disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s macs2_callpeak.py -c ${THREAD} -d ${DIR}
  #echo "snakemake -s plotFingerprint.py -c ${THREAD} -d ${DIR}"
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s plotFingerprint.py -c ${THREAD} -d ${DIR}
done