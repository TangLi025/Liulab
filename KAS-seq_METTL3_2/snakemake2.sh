#!/bin/bash
# creatation: 2022-1-13
# Author: Candy Lee tangli025@pku.edu.cn)

# Stop on error
set -e

THREAD=20

DIRS=( \
      "/disk1/home/user_09/KAS-METTL/METTL3_2" \
      )

for DIR in ${DIRS[@]}
do
  echo ${DIR}
  echo "snakemake -s bg2bw.py -c ${THREAD} -d ${DIR}"
  /disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s bg2bw.py -c ${THREAD} -d ${DIR}
  #echo "snakemake -s plotFingerprint.py -c ${THREAD} -d ${DIR}"
  #/disk1/home/user_09/anaconda3/envs/snakemake/bin/snakemake -s plotFingerprint.py -c ${THREAD} -d ${DIR}
done