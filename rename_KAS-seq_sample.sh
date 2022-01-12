TREATMENT=("input" "IP")
REP=("rep1" "rep2")
#GROUP=("METTL3_2" "METTL3_3")
GROUP=("METTL3_1" "METTL14")
SAMPLE=("CTRL" "KO")

input=(`seq 19 26`)
i=0
j=1
for group in ${GROUP[@]}
do
  for sample in ${SAMPLE[@]}
  do
    for rep in ${REP[@]}
    do
      echo mv ../KAS-METTL/CHe-Jliu-36S-K${input[i]}_S${input[i]}_R1_001.fastq.gz ../KAS-METTL/KAS-seq_${group}_${sampl_IP_${rep}.fastq.gz
      i=`expr $i + 1`
    done
  done
done
