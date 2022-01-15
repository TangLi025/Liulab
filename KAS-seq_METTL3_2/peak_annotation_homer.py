GROUP=["METTL3_2","METTL3_1","METTL3_3","METTL14","YTHDC1"]
SAMPLE=["CTRL","KO"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

MODE=["regular","broad"]

GENOME="/disk1/home/user_09/reference/genome/GRCm39.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCm39.genome.fa.fai"
GTF="/disk1/home/user_09/reference/annotation/mm19/gencode.vM28.annotation.gtf"

rule all:
  input:
    expand("{group}/06_macs2/peak_annotation/broad/KAS-seq_{group}_{sample}_{rep}_annotation_homer.xls",group=GROUP,sample=SAMPLE,rep=REP)

rule peak_annotation_homer_raw:
  input:
    "{group}/06_macs2/broad/KAS-seq_{group}_{sample}_{rep}_peaks.broadPeak"
  output:
    "{group}/06_macs2/peak_annotation/broad/KAS-seq_{group}_{sample}_{rep}_annotation_homer.xls"
  log:
    "logs/peak_annotation_homer/broad/KAS-seq_{group}_{sample}_{rep}_annotation_homer.log"
  threads:2
  params:
    genome=GENOME,
    gtf=GTF
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/annotatePeaks.pl {input} \
      {params.genome} -gtf {params.gtf} \
      -cpu {threads} \
      1> {output} 2> {log}"
