#!/bin/bash

set -e

KAS_seq_file=`ls ~/KAS-METTL/METTL3_2/05_bedtools/bedGraph/*.bg`
ratio=`cat ~/KAS-METTL/METTL3_2/04_bam_rmdup/bam_ratio.txt`

number_of_samples=8
number_of_ratios=8
echo "$number_of_samples"

if [[ ${number_of_samples} != ${number_of_ratios} ]]
then
   echo "error:the number of ratios isn't consistent with the number of samples"
   exit
fi


# parameters for KAS-seq_normalization.sh

#normalize bedGraph files
mkdir -p ~/KAS-METTL/METTL3_2/05_bedtools/bedGraph_nor
for ((i=1; i<=${number_of_samples}; i++))
do
    samples_selected=$(awk '{print $'$i' }' $KAS_seq_file)
    echo $i
    ratio_selected=$(awk '{print $'$i' }' $ratio)
    KAS_seq_basename=$(basename ${samples_selected} .bg) 
    echo "${KAS_seq_basename}"
    echo awk -v ratios=$ratio_selected '{printf("%s\t%d\t%d\t%.2f\n",$1,$2,$3,$4*ratios)}' $samples_selected > ~/KAS-METTL/METTL3_2/05_bedtools/bedGraph_nor/${KAS_seq_basename}.bg 
done