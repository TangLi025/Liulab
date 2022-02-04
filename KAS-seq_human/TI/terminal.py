BED_TERMINAL="/disk1/home/user_09/reference/annotation/hg19/termination_peak_overlap.bed"
BLACKLIST="/disk1/home/user_09/reference/annotation/hg19/hg19_blacklist.bed"

rule all:
  input:
    "07_deeptools/computeMatrix/DMSO_IP_merge_terminal.tab",
    "07_deeptools/plotProfile_terminal/terminal.png"
    
rule computeMatrix_distribution_terminal:
  input:
    bw="07_deeptools/bamCoverage/DMSO_IP_merge.bw"
  output:
    mat="07_deeptools/computeMatrix/DMSO_IP_merge_terminal.mat.gz",
    tab="07_deeptools/computeMatrix/DMSO_IP_merge_terminal.tab"
  log:
    "logs/computeMatrix_distribution/DMSO_IP_merge_terminal.log"
  params:
    gtf=BED_TERMINAL
  threads: 45
  shell:
   "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix scale-regions \
      --regionsFileName {params.gtf} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --regionBodyLength 10000 \
      --beforeRegionStartLength 0  \
      --afterRegionStartLength 0 \
      --binSize 500 \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"
       
rule computeMatrix_distribution:
  input:
    bw="07_deeptools/bamCoverage/DMSO_IP_merge.bw"
  output:
    mat="07_deeptools/computeMatrix_terminal/terminal.mat.gz",
  log:
    "logs/computeMatrix_distribution_terminal/terminal.log"
  params:
    gtf=["~/KAS/TI_KAS-seq/long.gtf","~/KAS/TI_KAS-seq/medium.gtf","~/KAS/TI_KAS-seq/short.gtf"],
    blacklist=BLACKLIST
  threads: 45
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix reference-point \
      --regionsFileName {params.gtf} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --beforeRegionStartLength 2000  \
      --afterRegionStartLength 10000 \
      --binSize 10 \
      --skipZeros \
      --blackListFileName {params.blacklist} \
      --smartLabels \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"

rule plotProfile_terminal:
  input:
    mat="07_deeptools/computeMatrix_terminal/terminal.mat.gz"
  output:
    png="07_deeptools/plotProfile_terminal/terminal.png"
  log:
    "logs/plotProfile_terminal/terminal.log"
  params:
    labels=["long","medium","short"]
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotProfile \
      --matrixFile {input.mat} \
      --regionsLabel {params.labels} \
      --outFileName {output.png} \
       > {log} 2>&1"
