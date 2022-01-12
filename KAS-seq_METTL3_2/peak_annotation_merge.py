TREATMENT=["input","IP"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]
SHUFFLE=["","_shuffled"]
SAMPLE=["DMSO","DRB","TRIP"]
#SAMPLE=["DMSO"]

GENOME="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa.fai"
GTF="/disk1/home/user_09/reference/annotation/hg19/gencode.v19.annotation.gtf"
CHROM_SIZE="/disk1/home/user_09/reference/annotation/hg19/hg19.chrom.sizes"

rule all:
  input:
    expand("06_macs2/bed_merge/{sample}_common_peaks.broadPeak",sample=SAMPLE)
    #"06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_DMSO_merge.pdf",

rule bed_merge:
  input:
    "06_macs2/{sample}_rep1_peaks.broadPeak",
    "06_macs2/{sample}_rep2_peaks.broadPeak"
  output:
    "06_macs2/bed_merge/{sample}_rep1_unique_peaks.broadPeak",
    "06_macs2/bed_merge/{sample}_rep2_unique_peaks.broadPeak",
    "06_macs2/bed_merge/{sample}_rep1_common_peaks.broadPeak",
    "06_macs2/bed_merge/{sample}_rep2_common_peaks.broadPeak",
    "06_macs2/bed_merge/{sample}_common_peaks.broadPeak"
  log:
    "logs/bed_merge/{sample}.log"
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools intersect -a {input[0]} -b {input[1]} -v \
    > {output[0]}
    echo "{output[0]}:" > {log} 2>&1
    cat {output[0]} | wc -l >> {log} 2>&1
  
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools intersect -a {input[1]} -b {input[0]} -v \
    > {output[1]}
    echo "{output[1]}:" >> {log} 2>&1
    cat {output[1]} | wc -l >> {log} 2>&1
  
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools intersect -a {input[0]} -b {output[0]} -v \
    > {output[2]}
    echo "{output[2]}:" >> {log} 2>&1
    cat {output[2]} | wc -l >> {log} 2>&1
    
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools intersect -a {input[1]} -b {output[1]} -v \
    > {output[3]}
    echo "{output[3]}:" >> {log} 2>&1
    cat {output[3]} | wc -l >> {log} 2>&1
  
    cat {output[2]} {output[3]} \
    | sort -k1,1 -k2,2n \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools merge \
    | awk -F "\t" -v OFS="\t" '{{print $1,$2,$3,NR,$5,$6,$7,$8,$9}}' \
    > {output[4]}
    echo "{output[4]}:" >> {log} 2>&1
    cat {output[4]} | wc -l >> {log} 2>&1
    """

rule bedtools_shuffle_merge:
  input:
    "06_macs2/bed_merge/{sample}_common_peaks.broadPeak",
    GENOME_INDEX
  output:
    "06_macs2/bed_merge/{sample}_common_peaks_shuffled.broadPeak"
  log:
    "logs/bedtools_shuffle/{sample}_merge.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools shuffle -i {input[0]} \
      -g {input[1]} 1> {output} 2> {log}"
  
rule peak_annotation_homer_merge:
  input:
    "06_macs2/bed_merge/{sample}_common_peaks.broadPeak"
  output:
    "06_macs2/peak_annotation/{sample}_annotation_homer_merge.xls"
  log:
    "logs/peak_annotation_homer/{sample}_annotation_homer_merge.log"
  threads:10
  params:
    genome=GENOME,
    gtf=GTF
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/annotatePeaks.pl {input} \
      {params.genome} -gtf {params.gtf} \
      -cpu {threads} \
      1> {output} 2> {log}"
      
rule peak_annotation_homer_merge_shuffled:
  input:
    "06_macs2/bed_merge/{sample}_common_peaks_shuffled.broadPeak"
  output:
    "06_macs2/peak_annotation/{sample}_annotation_homer_merge_shuffled.xls"
  log:
    "logs/peak_annotation_homer/{sample}_annotation_homer_merge_shuffled.log"
  threads:10
  params:
    genome=GENOME,
    gtf=GTF
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/annotatePeaks.pl {input} \
      {params.genome} -gtf {params.gtf} \
      -cpu {threads} \
      1> {output} 2> {log}"

rule Distribution_of_m6A_peaks_rep1:
  input:
    "06_macs2/peak_annotation/{sample}_annotation_homer_merge_shuffled.xls",
    "06_macs2/peak_annotation/{sample}_annotation_homer_merge.xls"
  output:
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_{sample}_merge.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_bar_{sample}_merge.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_percent_{sample}_merge.txt"
  log:
    "logs/Distribution_of_m6A_peaks/Distribution_of_m6A_peaks_{sample}_rep1.log" 
  params:
    sample=["Random","KAS-seq"]
  threads:4
  script:
    "scripts/Distribution_of_m6A_peaks.R"
