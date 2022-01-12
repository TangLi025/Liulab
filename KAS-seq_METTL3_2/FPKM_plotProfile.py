BLACKLIST="/disk1/home/user_09/reference/annotation/hg19/hg19_blacklist.bed"


rule all:
  input:
    "07_deeptools/plotProfile_FPKM/DMSO.png"

rule computeMatrix_distribution_FPKM:
  input:
    bw="07_deeptools/bamCoverage/DMSO_IP_merge.bw"
  output:
    mat="07_deeptools/computeMatrix_FPKM/{sample}.mat.gz",
  log:
    "logs/computeMatrix_distribution_FPKM/{sample}.log"
  params:
    gtf=["~/KAS/RNA-seq/05_FPKM_counts/fpkms_high.gtf","~/KAS/RNA-seq/05_FPKM_counts/fpkms_medium.gtf","~/KAS/RNA-seq/05_FPKM_counts/fpkms_low.gtf","~/KAS/RNA-seq/05_FPKM_counts/fpkms_silent.gtf"],
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
    mat="07_deeptools/computeMatrix_FPKM/{sample}.mat.gz"
  output:
    png="07_deeptools/plotProfile_FPKM/{sample}.png"
  log:
    "logs/plotProfile_FPKM/{sample}.log"
  params:
    labels=["high","medium","low","silent"]
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --matrixFile {input.mat} \
      --regionsLabel {params.labels} \
      --outFileName {output.png} \
       > {log} 2>&1"
