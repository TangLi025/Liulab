GROUP=["METTL3_2","METTL3_3"]
SAMPLE=["CTRL","KO"]

TREATMENT=["IP"]
REP=["rep1","rep2"]


BLACKLIST="/disk1/home/user_09/reference/annotation/mm19/mm19.blacklist.bed"

HIGH="/disk1/home/user_09/KAS/fpkms_high.gtf"
MEDIUM="/disk1/home/user_09/KAS/fpkms_medium.gtf"
LOW="/disk1/home/user_09/KAS/fpkms_low.gtf"
SILENT="/disk1/home/user_09/KAS/fpkms_silent.gtf"

rule all:
  input:
    expand("{group}/07_deeptools/plotProfile_FPKM/{group}_liu.png",group=GROUP)

rule computeMatrix_distribution_FPKM:
  input:
    bw="{group}/05_bedtools/bigWig/KAS-seq_{group}_CTRL_IP_rep1_ext.bw"
  output:
    mat="{group}/07_deeptools/computeMatrix_FPKM/{group}_liu.mat.gz",
  log:
    "{group}/logs/computeMatrix_distribution_FPKM/{group}_liu.log"
  params:
    gtf=[HIGH,MEDIUM,LOW,SILENT],
    blacklist=BLACKLIST
  threads: 45
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
    mat="{group}/07_deeptools/computeMatrix_FPKM/{group}_liu.mat.gz"
  output:
    png="{group}/07_deeptools/plotProfile_FPKM/{group}_liu.png"
  log:
    "{group}/logs/plotProfile_FPKM/{group}_liu.log"
  params:
    labels=["high","medium","low","silent"]
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --matrixFile {input.mat} \
      --regionsLabel {params.labels} \
      --outFileName {output.png} \
       > {log} 2>&1"
