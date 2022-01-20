GROUP=["METTL3_1","METTL3_2","METTL3_3","METTL14","YTHDC1"]
TREATMENT=["input","IP"]
REP=["rep1","rep2"]
SAMPLE=["CTRL","KO"]

GENOME_2bit="/disk1/home/user_09/reference/genome/GRCm39.genome.2bit"
GENOME="/disk1/home/user_09/reference/genome/GRCm39.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCm39.genome.fa.fai"
GTF="/disk1/home/user_09/reference/annotation/mm19/gencode.vM28.annotation.gtf"
GTF_protein_coding="/disk1/home/user_09/reference/annotation/mm19/gencode.vM28.annotation.protein_coding.chr.gtf"
BED_GB="/disk1/home/user_09/reference/annotation/mm19/mm19_Refseq.GB.bed"
BED_TSS="/disk1/home/user_09/reference/annotation/mm19/mm19_Refseq.TSS.bed"
BLACKLIST="/disk1/home/user_09/reference/annotation/mm19/mm19.blacklist.bed"

rule all:
  input:
    expand("{group}/07_deeptools/computeMatrix/{sample}_{treatment}_{rep}_TSS.tab",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP),
    expand("{group}/07_deeptools/computeMatrix/{sample}_{treatment}_{rep}_body.tab",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP)

rule computeMatrix_distribution_TSS:
  input:
    bw="{group}/05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw"
  output:
    mat="{group}/07_deeptools/computeMatrix/{sample}_{treatment}_{rep}_TSS.mat.gz",
    tab="{group}/07_deeptools/computeMatrix/{sample}_{treatment}_{rep}_TSS.tab",
    bed="{group}/07_deeptools/computeMatrix/{sample}_{treatment}_{rep}_TSS.bed"
  log:
    "{group}/logs/computeMatrix_distribution/{sample}_{treatment}_{rep}_TSS.log"
  params:
    gtf=BED_TSS,
    blacklist=BLACKLIST,
  threads: 10
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
       
rule computeMatrix_distribution_gene_body:
  input:
    bw="{group}/05_bedtools/bigWig/KAS-seq_{group}_{sample}_{treatment}_{rep}_ext.bw"
  output:
    mat="{group}/07_deeptools/computeMatrix/{sample}_{treatment}_{rep}_body.mat.gz",
    tab="{group}/07_deeptools/computeMatrix/{sample}_{treatment}_{rep}_body.tab",
    bed="{group}/07_deeptools/computeMatrix/{sample}_{treatment}_{rep}_body.bed"
  log:
    "{group}/logs/computeMatrix_distribution/{sample}_{treatment}_{rep}_body.log"
  params:
    gtf=BED_GB,
    blacklist=BLACKLIST,
  threads: 10
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
