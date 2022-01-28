GROUP=["METTL3_2","METTL3_3"]
SAMPLE=["CTRL","KO"]

TREATMENT=["IP","input"]
REP=["rep1","rep2"]


BLACKLIST="/disk1/home/user_09/reference/annotation/mm19/mm19.blacklist.bed"

FPKM=["high","medium","low","silent"]
HIGH="/disk1/home/user_09/KAS/fpkm_m6A_Liu/fpkms_high.gtf"
MEDIUM="/disk1/home/user_09/KAS/fpkm_m6A_Liu/fpkms_medium.gtf"
LOW="/disk1/home/user_09/KAS/fpkm_m6A_Liu/fpkms_low.gtf"
SILENT="/disk1/home/user_09/KAS/fpkm_m6A_Liu/fpkms_silent.gtf"

rule all:
  input:
    expand("{group}/07_deeptools/plotProfile_FPKM/{group}_fpkm_liu.png",group=GROUP)

rule computeMatrix_distribution_FPKM:
  input:
    bw=expand("{group}/05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw",group=r'{group}',sample=SAMPLE,treatment=TREATMENT,rep=REP)
  output:
    mat="{group}/07_deeptools/computeMatrix_FPKM/{group}_fpkm_liu.mat.gz",
  log:
    "{group}/logs/computeMatrix_distribution_FPKM/{group}_fpkm_liu.log"
  params:
    gtf=[HIGH,MEDIUM,LOW,SILENT],
    blacklist=BLACKLIST
  threads: 100
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix scale-regions \
      --regionsFileName {params.gtf} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --regionBodyLength 6000 \
      --beforeRegionStartLength 3000  \
      --afterRegionStartLength 3000 \
      --binSize 10 \
      --skipZeros \
      --blackListFileName {params.blacklist} \
      --smartLabels \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"

rule plotProfile_FPKM:
  input:
    mat="{group}/07_deeptools/computeMatrix_FPKM/{group}_fpkm_liu.mat.gz"
  output:
    png="{group}/07_deeptools/plotProfile_FPKM/{group}_fpkm_liu.png"
  log:
    "{group}/logs/plotProfile_FPKM/{group}_fpkm_liu.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --perGroup \
      --matrixFile {input.mat} \
      --outFileName {output.png} \
       > {log} 2>&1"
