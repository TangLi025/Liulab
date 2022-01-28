

GROUP=["METTL3_2","METTL3_1","METTL3_3","METTL14","YTHDC1"]
SAMPLE=["CTRL","KO"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

MODE=["regular","broad"]

GENOME="/disk1/home/user_09/reference/genome/GRCm38.p6.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCm38.p6.genome.fa"
GTF="/disk1/home/user_09/reference/annotation/mm10/gencode.vM19.annotation.gtf"

rule all:
  input:
    "/disk1/home/user_09/carRNA_science/m6A-seq_narrowPeak/Mettl3_Control_commen_peaks_annotation_homer.xls"

rule peak_annotation_homer_raw:
  input:
    "/disk1/home/user_09/carRNA_science/m6A-seq_narrowPeak/Mettl3_Control_commen_peaks.narrowPeak"
  output:
    "/disk1/home/user_09/carRNA_science/m6A-seq_narrowPeak/Mettl3_Control_commen_peaks_annotation_homer.xls"
  log:
    "/disk1/home/user_09/carRNA_science/logs/peak_annotation_homer/Mettl3_Control_commen_peaks_annotation_homer.log"
  threads:2
  params:
    genome=GENOME,
    gtf=GTF
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/annotatePeaks.pl {input} \
      {params.genome} -gtf {params.gtf} \
      -cpu {threads} \
      1> {output} 2> {log}"
