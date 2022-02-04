TREATMENT=["input","IP"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]
SAMPLE=["DMSO"]
#SAMPLE=["DMSO","DRB","TRIP"]
MERGE=["","_merge"]

GENOME_2bit="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.2bit"
GENOME="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa.fai"
GTF="/disk1/home/user_09/reference/annotation/hg19/gencode.v19.annotation.gtf"
GTF_protein_coding="/disk1/home/user_09/reference/annotation/hg19/gencode.v19.annotation.protein_coding.chr.gtf"
BED_GB="/disk1/home/user_09/reference/annotation/hg19/gencode.v19.annotation.protein_coding.chr.gene.GB.bed"
BED_TSS="/disk1/home/user_09/reference/annotation/hg19/gencode.v19.annotation.protein_coding.chr.gene.TSS.bed"
BLACKLIST="/disk1/home/user_09/reference/annotation/hg19/hg19_blacklist.bed"

rule all:
  input:
    #"07_deeptools/computeMatrix/DMSO_IP_merge.tab",
    "07_deeptools/computeMatrix/DMSO_IP_merge_TSS.tab",
    "07_deeptools/computeMatrix/DMSO_IP_merge_body_1000.tab",
    "07_deeptools/computeMatrix/DMSO_IP_merge_body_10000.tab"

rule computeMatrix_distribution_TSS:
  input:
    bw="07_deeptools/bamCoverage/DMSO_IP_merge.bw"
  output:
    mat="07_deeptools/computeMatrix/DMSO_IP_merge_TSS.mat.gz",
    tab="07_deeptools/computeMatrix/DMSO_IP_merge_TSS.tab",
    bed="07_deeptools/computeMatrix/DMSO_IP_merge_TSS.bed"
  log:
    "logs/computeMatrix_distribution/DMSO_IP_merge_TSS.log"
  params:
    gtf=BED_TSS,
    blacklist=BLACKLIST,
  threads: 45
  shell:
   "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix scale-regions \
      --regionsFileName {params.gtf} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --regionBodyLength 600 \
      --beforeRegionStartLength 0  \
      --afterRegionStartLength 0 \
      --binSize 50 \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"
       
rule computeMatrix_distribution_gene_body_1000:
  input:
    bw="07_deeptools/bamCoverage/DMSO_IP_merge.bw"
  output:
    mat="07_deeptools/computeMatrix/DMSO_IP_merge_body_1000.mat.gz",
    tab="07_deeptools/computeMatrix/DMSO_IP_merge_body_1000.tab",
    bed="07_deeptools/computeMatrix/DMSO_IP_merge_body_1000.bed"
  log:
    "logs/computeMatrix_distribution/DMSO_IP_merge_body_1000.log"
  params:
    gtf=BED_GB,
    blacklist=BLACKLIST,
  threads: 45
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix scale-regions \
      --regionsFileName {params.gtf} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --regionBodyLength 1000 \
      --beforeRegionStartLength 0  \
      --afterRegionStartLength 0 \
      --binSize 50 \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"
       
rule computeMatrix_distribution_gene_body_10000:
  input:
    bw="07_deeptools/bamCoverage/DMSO_IP_merge.bw"
  output:
    mat="07_deeptools/computeMatrix/DMSO_IP_merge_body_10000.mat.gz",
    tab="07_deeptools/computeMatrix/DMSO_IP_merge_body_10000.tab",
    bed="07_deeptools/computeMatrix/DMSO_IP_merge_body_10000.bed"
  log:
    "logs/computeMatrix_distribution/DMSO_IP_merge_body_10000.log"
  params:
    gtf=BED_GB,
    blacklist=BLACKLIST,
  threads: 45
  shell:
    "/disk1/home/user_09/anaconda3/envs/deeptools/bin/computeMatrix scale-regions \
      --regionsFileName {params.gtf} \
      --scoreFileName {input.bw} \
      --outFileName {output.mat} \
      --outFileNameMatrix {output.tab} \
      --outFileSortedRegions {output.bed} \
      --regionBodyLength 10000 \
      --beforeRegionStartLength 0  \
      --afterRegionStartLength 0 \
      --binSize 50 \
      --numberOfProcessors {threads} \
      --missingDataAsZero \
       > {log} 2>&1"
