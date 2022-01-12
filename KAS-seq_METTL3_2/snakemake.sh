#!/bin/bash
# creatation: 2020-1-14
# Author: Ruitu Lyu (lvruitu@gmail.com)

# Stop on error
set -e
###
### setup.sh - This script is used to setup the KAS-seq analysis pipeline.
###
### make sure you have conda installed in your server or PCs.
###
### you can follow the user guide to accomplish the anaconda installation: https://docs.conda.io/projects/conda/en/latest/user-guide/install/.
###
### Usage: ./setup.sh	
###
### -h or --help Print the help.
###

# Help message for shell scripts

help() { 
    sed -rn 's/^### ?//;T;p' "$0"
}

THREAD=50

DIRS=( \
      "/disk1/home/user_09/LinLong/" \
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