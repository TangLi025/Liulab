TREATMENT=["input","IP"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]
SHUFFLE=["","_shuffled"]
#SAMPLE=["DMSO","DRB","TRIP"]
SAMPLE=["DMSO"]

GENOME="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa.fai"
GTF="/disk1/home/user_09/reference/annotation/hg19/gencode.v19.annotation.gtf"
CHROM_SIZE="/disk1/home/user_09/reference/annotation/hg19/hg19.chrom.sizes"

rule all:
  input:
    expand("06_macs2/peak_annotation/{sample}_{rep}_annotation_homer_raw.xls",sample=SAMPLE,rep=REP),
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_DMSO_raw.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_DMSO_merge.pdf"
    
    
rule bedtools_shuffle:
  input:
    "06_macs2/{sample}_rep1_peaks.broadPeak",
    GENOME_INDEX
  output:
    "06_macs2/{sample}_rep1_peaks_shuffled.broadPeak"
  log:
    "logs/bedtools_shuffle/{sample}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools shuffle -i {input[0]} \
      -g {input[1]} 1> {output} 2> {log}"
  
rule peak_annotation_homer_raw:
  input:
    "06_macs2/{sample}_rep1_peaks.broadPeak"
  output:
    "06_macs2/peak_annotation/{sample}_{rep}_annotation_homer_raw.xls"
  log:
    "logs/peak_annotation_homer/{sample}_{rep}_annotation_homer_raw.log"
  threads:10
  params:
    genome=GENOME,
    gtf=GTF
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/annotatePeaks.pl {input} \
      {params.genome} -gtf {params.gtf} \
      -cpu {threads} \
      1> {output} 2> {log}"
      
rule peak_annotation_homer_raw_shuffled:
  input:
    "06_macs2/{sample}_rep1_peaks_shuffled.broadPeak"
  output:
    "06_macs2/peak_annotation/{sample}_rep1_annotation_homer_raw_shuffled.xls"
  log:
    "logs/peak_annotation_homer/{sample}_rep1_annotation_homer_raw_shuffled.log"
  threads:10
  params:
    genome=GENOME,
    gtf=GTF
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/annotatePeaks.pl {input} \
      {params.genome} -gtf {params.gtf} \
      -cpu {threads} \
      1> {output} 2> {log}"

rule Distribution_of_m6A_peaks_raw:
  input:
    "06_macs2/peak_annotation/{sample}_rep1_annotation_homer_raw_shuffled.xls",
    "06_macs2/peak_annotation/{sample}_rep1_annotation_homer_raw.xls"
  output:
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_{sample}_raw.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_bar_{sample}_raw.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_percent_{sample}_raw.txt"
  log:
    "logs/Distribution_of_m6A_peaks/Distribution_of_m6A_peaks_{sample}_raw.log" 
  params:
    sample=["Random","KAS-seq"]
  threads:4
  script:
    "scripts/Distribution_of_m6A_peaks.R"

rule Distribution_of_m6A_peaks_modified:
  input:
    "06_macs2/peak_annotation/{sample}_rep1_annotation_homer_modified_shuffled.xls",
    "06_macs2/peak_annotation/{sample}_rep1_annotation_homer_modified.xls"
  output:
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_{sample}_modified.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_bar_{sample}_modified.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_percent_{sample}_modified.txt"
  log:
    "logs/Distribution_of_m6A_peaks/Distribution_of_m6A_peaks_{sample}_modified.log" 
  params:
    sample=["Random","KAS-seq"]
  threads:4
  script:
    "scripts/Distribution_of_m6A_peaks.R"
    
rule Distribution_of_m6A_peaks_merge:
  input:
    "06_macs2/peak_annotation/{sample}_rep1_annotation_homer_merge_shuffled.xls",
    "06_macs2/peak_annotation/{sample}_rep1_annotation_homer_merge.xls"
  output:
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_{sample}_merge.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_bar_{sample}_merge.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_percent_{sample}_merge.txt"
  log:
    "logs/Distribution_of_m6A_peaks/Distribution_of_m6A_peaks_{sample}_merge.log" 
  params:
    sample=["Random","KAS-seq"]
  threads:4
  script:
    "scripts/Distribution_of_m6A_peaks.R"


  
