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
    #expand("05_bedtools/{sample}_{treatment}_{rep}_ext.bw",sample=SAMPLE,treatment=TREATMENT,rep=REP),
    #expand("06_macs2/{sample}_{rep}_peaks.xls",sample=SAMPLE,rep=REP),
    #expand("06_macs2/{sample}_{rep}_peaks.broadPeak",sample=SAMPLE,rep=REP),
    #expand("06_macs2/{sample}_{rep}_peaks_shuffled.broadPeak",sample=SAMPLE,rep=REP),
    #expand("06_macs2/peak_annotation/{sample}_{rep}{shuffle}_annotation_homer.xls",sample=SAMPLE,rep=REP,shuffle=SHUFFLE),
    #expand("06_macs2/peak_annotation/{sample}_{rep}_annotation_homer_raw.xls",sample=SAMPLE,rep=REP,shuffle=SHUFFLE),
    #expand("06_macs2/bed_merge/{sample}_common_peaks.bed",sample=SAMPLE),
    #expand("06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_{sample}.pdf",sample=SAMPLE),
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_DMSO_raw.pdf",
    #expand("06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_{sample}_{rep}.pdf",sample=SAMPLE,rep=REP),
    GENOME_INDEX

rule genome_index:
  input:
    GENOME
  output:
    GENOME_INDEX
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/samtools faidx {input}"

rule bam2bed:
  input:
    "04_bam_rmdup/{sample}_{treatment}_{rep}.bam"
  output:
    "05_bedtools/{sample}_{treatment}_{rep}_ext.bed"
  log:
    "logs/bam2bed/{sample}_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools bamtobed -i {input} \
      | sort -k 1,1 1> {output} 2> {log}"

#rule bed_extend:
#  input:
#    "05_bedtools/{sample}_{treatment}_{rep}.bed",
#    CHROM_SIZE
#  output:
#    "05_bedtools/{sample}_{treatment}_{rep}_ext.bed"
#  log:
#    "logs/bed_extend/{sample}_{treatment}_{rep}.log"
##  script:
#    "script/bed_extend.R"



rule bed2bedGraph:
  input:
    "05_bedtools/{sample}_{treatment}_{rep}_ext.bed",
    GENOME_INDEX
  output:
    "05_bedtools/{sample}_{treatment}_{rep}_ext.bg"
  log:
    "logs/bed2bedGraph/{sample}_{treatment}_{rep}.log"
  params:
    genome_index=GENOME_INDEX
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools genomecov -i {input[0]} \
      -bg -g {input[1]} \
      | /disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools sort -i 1> {output} 2> {log}"

rule bedGraphToBigWig:
  input:
    "05_bedtools/{sample}_{treatment}_{rep}_ext.bg",
    GENOME_INDEX
  output:
    "05_bedtools/{sample}_{treatment}_{rep}_ext.bw"
  log:
    "logs/bedGraphToBigWig/{sample}_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedGraphToBigWig {input[0]}\
      {input[1]} {output} > {log} 2>&1"

rule macs2_callpeak:
  input:
    "05_bedtools/{sample}_IP_{rep}_ext.bed",
    "05_bedtools/{sample}_input_{rep}_ext.bed"
  output:
    "06_macs2/{sample}_{rep}_peaks.broadPeak",
    "06_macs2/{sample}_{rep}_peaks.xls"
  log:
    "logs/macs2_callpeak/{sample}_{rep}.log"
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/macs2 callpeak \
      -t {input[0]} -c {input[1]} \
      -n {wildcards.sample}_{wildcards.rep} \
      --verbose 3 --extsize 150 \
      -B \
      -g hs -q 0.01 --broad --broad-cutoff 0.01 \
      --outdir 06_macs2 > {log} 2>&1"

rule bed_merge:
  input:
    "06_macs2/{sample}_rep1_peaks.broadPeak",
    "06_macs2/{sample}_rep2_peaks.broadPeak"
  output:
    "06_macs2/bed_merge/{sample}_rep1_unique_peaks.bed",
    "06_macs2/bed_merge/{sample}_rep2_unique_peaks.bed",
    "06_macs2/bed_merge/{sample}_rep1_common_peaks.bed",
    "06_macs2/bed_merge/{sample}_rep2_common_peaks.bed",
    "06_macs2/bed_merge/{sample}_common_peaks.bed"
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
    | awk -F "\t" -v OFS="\t" '{{print $1,$2,$3,NR,".","."}}' \
    > {output[4]}
    echo "{output[4]}:" >> {log} 2>&1
    cat {output[4]} | wc -l >> {log} 2>&1
    """
    
rule bedtools_shuffle:
  input:
    "06_macs2/{sample}_rep1_peaks.broadPeak",
    GENOME_INDEX
  output:
    "06_macs2/{sample}_rep1_peaks_shuffled.broakPeak"
  log:
    "logs/bedtools_shuffle/{sample}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools shuffle -i {input[0]} \
      -g {input[1]} 1> {output} 2> {log}"
      
rule peak_annotation_homer_raw:
  input:
    "06_macs2/{sample}_{rep}_peaks.broadPeak"
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

rule peak_annotation_homer_merge:
  input:
    "06_macs2/bed_merge/{sample}_rep1_common_peaks.bed"
  output:
    "06_macs2/peak_annotation/{sample}_rep1_annotation_homer_merge.xls"
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

rule peak_shuffled_annotation_homer:
  input:
    "06_macs2/bed_merge/{sample}_rep1_common_peaks_shuffled.bed"
  output:
    "06_macs2/peak_annotation/{sample}_rep1_annotation_homer_merge_shuffled.xls"
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
      
rule Distribution_of_m6A_peaks:
  input:
    "06_macs2/peak_annotation/{sample}_rep1_annotation_homer_merge_shuffled.xls",
    "06_macs2/peak_annotation/{sample}_rep1_annotation_homer_merge.xls"
  output:
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_{sample}.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_bar_{sample}.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_percent_{sample}.txt"
  log:
    "logs/Distribution_of_m6A_peaks/Distribution_of_m6A_peaks_{sample}.log" 
  params:
    sample=["Random","KAS-seq"]
  threads:4
  script:
    "scripts/Distribution_of_m6A_peaks.R"
