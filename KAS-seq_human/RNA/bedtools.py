TREATMENT=["input","IP"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]
SHUFFLE=["","_shuffled"]

GENOME="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCh37.p13.genome.fa.fai"
GTF="/disk1/home/user_09/reference/annotation/gencode.v19.annotation.gtf"

rule all:
  input:
    expand("05_bedtools/HEK_{treatment}_{rep}.bw",treatment=TREATMENT,rep=REP),
    expand("06_macs2/HEK_{rep}_peaks.xls",rep=REP),
    expand("06_macs2/HEK_{rep}_peaks.broadPeak",rep=REP),
    expand("06_macs2/HEK_{rep}_peaks_shuffled.broadPeak",rep=REP),
    expand("06_macs2/peak_annotation/HEK_{rep}{shuffle}_annotation_homer.xls",rep=REP,shuffle=SHUFFLE),
    expand("06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_{rep}.pdf",rep=REP),
    GENOME_INDEX
    
rule bam2bed:
  input:
    "04_bam_rmdup/HEK_{treatment}_{rep}.bam"
  output:
    "05_bedtools/HEK_{treatment}_{rep}.bed"
  log:
    "logs/bam2bed/HEK_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools bamtobed -i {input} \
      | sort -k 1,1 1> {output} 2> {log}"

rule genome_index:
  input:
    GENOME
  output:
    GENOME_INDEX
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/samtools faidx {input}"

rule bed2bedGraph:
  input:
    "05_bedtools/HEK_{treatment}_{rep}.bed",
    GENOME_INDEX
  output:
    "05_bedtools/HEK_{treatment}_{rep}.bg"
  log:
    "logs/bed2bedGraph/HEK_{treatment}_{rep}.log"
  params:
    genome_index=GENOME_INDEX
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools genomecov -i {input[0]} \
      -bg -g {input[1]} \
      | /disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools sort -i 1> {output} 2> {log}"

rule bedGraphToBigWig:
  input:
    "05_bedtools/HEK_{treatment}_{rep}.bg",
    GENOME_INDEX
  output:
    "05_bedtools/HEK_{treatment}_{rep}.bw"
  log:
    "logs/bedGraphToBigWig/HEK_{treatment}_{rep}.log"
  params:
    genome_index=GENOME_INDEX
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedGraphToBigWig {input[0]}\
      {input[1]} {output} > {log} 2>&1"

rule macs2_callpeak:
  input:
    "05_bedtools/HEK_IP_{rep}.bed",
    "05_bedtools/HEK_input_{rep}.bed"
  output:
    "06_macs2/HEK_{rep}_peaks.xls",
    "06_macs2/HEK_{rep}_peaks.broadPeak"
  log:
    "logs/macs2_callpeak/HEK_{rep}.log"
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/macs2 callpeak \
      -t {input[0]} -c {input[1]} \
      -n HEK_{wildcards.rep} \
      --verbose 3 --extsize 150 \
      -B \
      -g hs -q 0.01 --broad --broad-cutoff 0.01 \
      --outdir 06_macs2 > {log} 2>&1"

rule bedtools_shuffle:
  input:
    "06_macs2/HEK_{rep}_peaks.broadPeak",
    GENOME_INDEX
  output:
    "06_macs2/HEK_{rep}_peaks_shuffled.broadPeak"
  log:
    "logs/bedtools_shuffle/HEK_{rep}.log"
  params:
    genome_index=GENOME_INDEX
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools shuffle -i {input[0]} \
      -g {input[1]} 1> {output} 2> {log}"
      
rule peak_annotation_homer:
  input:
    "06_macs2/HEK_{rep}_peaks.broadPeak"
  output:
    "06_macs2/peak_annotation/HEK_{rep}_annotation_homer.xls"
  log:
    "logs/peak_annotation_homer/HEK_{rep}_annotation_homer.log"
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
    "06_macs2/HEK_{rep}_peaks_shuffled.broadPeak"
  output:
    "06_macs2/peak_annotation/HEK_{rep}_shuffled_annotation_homer.xls"
  log:
    "logs/peak_annotation_homer/HEK_{rep}_shuffled_annotation_homer.log"
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
    "06_macs2/peak_annotation/HEK_{rep}_annotation_homer.xls",
    "06_macs2/peak_annotation/HEK_{rep}_shuffled_annotation_homer.xls"
  output:
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie_{rep}.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_bar_{rep}.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_percent_{rep}.txt"
  log:
    "logs/Distribution_of_m6A_peaks/Distribution_of_m6A_peaks_{rep}.log" 
  threads:4
  script:
    "scripts/Distribution_of_m6A_peaks.R"
