SAMPLE=["P0","P7","P10","rP2o","rP2n"]

REP=["rep1","rep2"]

GENOME="/disk1/home/user_09/reference/genome/GRCm39.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCm39.genome.fa.fai"

 
BLACKLIST="/disk1/home/user_09/reference/annotation/mm19/KAS-seq/mm19.blacklist.bed"

rule all:
  input:
    "07_deeptools/plotCorrelation/AP.png",
    "07_deeptools/plotCorrelation/AP_heatmap.pdf",
    "07_deeptools/plotPCA_m6A/AP.png"

rule multiBamSummary:
  input:
    expand("04_bam_raw/TC1-AP_{sample}_rep1.bam",sample=SAMPLE),
    expand("04_bam_raw/TC1-AP_{sample}_rep2.bam",sample=SAMPLE)
  output:
    npz="07_deeptools/multiBamSummary/AP.npz"
  params:
    blacklist=BLACKLIST
  log:
    "logs/multiBamSummary/AP.log"
  threads: 10
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/multiBamSummary bins \
      -b {input} \
      --outFileName {output.npz} \
      --binSize 10000 \
      --numberOfProcessors {threads} \
      --blackListFileName {params.blacklist} \
       > {log} 2>&1"

rule plotCorrelation:
  input:
    npz="07_deeptools/multiBamSummary/AP.npz"
  output:
    png="07_deeptools/plotCorrelation/AP.png",
    tab="07_deeptools/plotCorrelation/AP.tab"
  log:
    "logs/plotCorrelation/AP.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotCorrelation \
      --corData {input.npz} \
      --corMethod pearson \
      --whatToPlot scatterplot \
      --plotFile {output.png} \
      --outFileCorMatrix {output.tab} \
      --removeOutliers \
       > {log} 2>&1"
       
rule plotCorrelation_heatmap:
  input:
    npz="07_deeptools/multiBamSummary/AP.npz"
  output:
    pdf="07_deeptools/plotCorrelation/AP_heatmap.pdf",
    tab="07_deeptools/plotCorrelation/AP_heatmap.tab"
  log:
    "logs/plotCorrelation/AP_heatmap.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotCorrelation \
      --corData {input.npz} \
      --corMethod pearson \
      --whatToPlot heatmap \
      --plotFile {output.pdf} \
      --outFileCorMatrix {output.tab} \
      --removeOutliers \
       > {log} 2>&1"
       
rule plotPCA:
  input:
    "07_deeptools/multiBamSummary/AP.npz"
  output:
    "07_deeptools/plotPCA_m6A/AP.png"
  log:
    "logs/plotPCA_m6A/AP.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/plotPCA \
      -in {input} \
      --plotHeight 15 --plotWidth 10 \
      --plotFile {output} \
      > {log} 2>&1"
