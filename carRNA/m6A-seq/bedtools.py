TREATMENT=["input","IP"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]
SHUFFLE=["","_shuffled"]
SAMPLE=["Lysate","Result"]
STRAND=["neg","pos"]

GENOME="/disk1/home/user_09/reference/genome/NCBIM37.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/NCBIM37.genome.fa.fai"
GTF="/disk1/home/user_09/reference/annotation/mm10/gencode.vM1.annotation.gtf"

rule all:
  input:
    expand("05_bam_{dup}_separated/{sample}_{treatment}_{rep}_{strand}.bam.bai",dup=DUP,sample=SAMPLE,treatment=TREATMENT,rep=REP,strand=STRAND),
    #expand("05_bedtools/{sample}_{treatment}_{rep}.bw",treatment=TREATMENT,rep=REP),
    #expand("06_macs2/{sample}_{rep}_peaks.xls",sample=SAMPLE,rep=REP),
    #expand("06_macs2/{sample}_{rep}_peaks.broadPeak",rep=REP),
    #expand("06_macs2/{sample}_{rep}_peaks_shuffled.broadPeak",rep=REP),
    #expand("06_macs2/peak_annotation/{sample}_{rep}{shuffle}_annotation_homer.xls",rep=REP,shuffle=SHUFFLE),
    #expand("06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie.pdf",rep=REP),
    GENOME_INDEX
    
rule bam_separate:
  input:
    "04_bam_{dup}/{sample}_{treatment}_{rep}.bam",
    "04_bam_{dup}/{sample}_{treatment}_{rep}.bam.bai"
  output:
    temp("05_bam_{dup}_separated/{sample}_{treatment}_{rep}_83.bam"),
    temp("05_bam_{dup}_separated/{sample}_{treatment}_{rep}_163.bam"),
    temp("05_bam_{dup}_separated/{sample}_{treatment}_{rep}_99.bam"),
    temp("05_bam_{dup}_separated/{sample}_{treatment}_{rep}_147.bam"),
    "05_bam_{dup}_separated/{sample}_{treatment}_{rep}_pos.bam",
    "05_bam_{dup}_separated/{sample}_{treatment}_{rep}_neg.bam"
  threads: 1
  shell:
    """
  	/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -b -f 83 {input[0]} 1> {output[0]}
  	/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -b -f 163 {input[0]} 1> {output[1]}
  	/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -b -f 99 {input[0]} 1> {output[2]}
  	/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -b -f 147 {input[0]} 1> {output[3]}
  	/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools merge -@ {threads} {output[4]} {output[0]} {output[1]} 
  	/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools merge -@ {threads} {output[5]} {output[2]} {output[3]}
  	"""
    
rule bam_separated_index:
  input:
    "05_bam_{dup}_separated/{sample}_{treatment}_{rep}_{strand}.bam"
  output:
    "05_bam_{dup}_separated/{sample}_{treatment}_{rep}_{strand}.bam.bai"
  log:
    "logs/05_bam_{dup}_separated/bam_{strand}_index_{sample}_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"

rule bam2bed:
  input:
    "04_bam_rmdup/{smaple}_{treatment}_{rep}.bam"
  output:
    "05_bedtools/{smaple}_{treatment}_{rep}.bed"
  log:
    "logs/bam2bed/{smaple}_{treatment}_{rep}.log"
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
    "05_bedtools/{smaple}_{treatment}_{rep}.bed",
    GENOME_INDEX
  output:
    "05_bedtools/{smaple}_{treatment}_{rep}.bg"
  log:
    "logs/bed2bedGraph/{smaple}_{treatment}_{rep}.log"
  params:
    genome_index=GENOME_INDEX
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools genomecov -i {input[0]} \
      -bg -g {input[1]} \
      | /disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools sort -i 1> {output} 2> {log}"

rule bedGraphToBigWig:
  input:
    "05_bedtools/{smaple}_{treatment}_{rep}.bg",
    GENOME_INDEX
  output:
    "05_bedtools/{smaple}_{treatment}_{rep}.bw"
  log:
    "logs/bedGraphToBigWig/{smaple}_{treatment}_{rep}.log"
  params:
    genome_index=GENOME_INDEX
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedGraphToBigWig {input[0]}\
      {input[1]} {output} > {log} 2>&1"

rule macs2_callpeak:
  input:
    "05_bedtools/{sample}_IP_{rep}.bed",
    "05_bedtools/{sample}_input_{rep}.bed"
  output:
    "06_macs2/{sample}_{rep}_summits.bed",
    "06_macs2/{sample}_{rep}_peaks.xls",
    "06_macs2/{sample}_{rep}_peaks.narrowPeak"
  log:
    "logs/macs2_callpeak/{sample}_{rep}.log"
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/macs2 callpeak \
      -t {input[0]} -c {input[1]} \
      -n {wildcards.sample}_{wildcards.rep} \
      --verbose 3 --extsize 150 \
      -g 1.3e8 -q 0.01 \
      -B \
      --outdir 06_macs2 > {log} 2>&1"

rule bed_merge:
  input:
    "06_macs2/{sample}_rep1_peaks.narrowPeak",
    "06_macs2/{sample}_rep2_peaks.narrowPeak"
  output:
    "06_macs2/{sample}_rep1_unique_peaks.bed",
    "06_macs2/{sample}_rep2_unique_peaks.bed",
    "06_macs2/{sample}_rep1_common_peaks.bed",
    "06_macs2/{sample}_rep2_common_peaks.bed",
    "06_macs2/{sample}_common_peaks.bed"
  log:
    "logs/bed_merge/{sample}.log"
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools intersect -a {input[0]} -b {input[1]} -v -s \
    > {output[0]}
    echo "{output[0]}:" > {log} 2>&1
    cat {output[0]} | wc -l >> {log} 2>&1
  
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools intersect -a {input[1]} -b {input[0]} -v -s \
    > {output[1]}
    echo "{output[1]}:" >> {log} 2>&1
    cat {output[1]} | wc -l >> {log} 2>&1
  
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools intersect -a {input[0]} -b {output[0]} -v -s \
    > {output[2]}
    echo "{output[2]}:" >> {log} 2>&1
    cat {output[2]} | wc -l >> {log} 2>&1
  
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools intersect -a {input[1]} -b {output[1]} -v -s \
    > {output[3]}
    echo "{output[3]}:" >> {log} 2>&1
    cat {output[3]} | wc -l >> {log} 2>&1
  
    cat {output[2]} {output[3]} \
    | sort -k1,1 -k2,2n \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools merge -s -c 6 -o distinct \
    | awk -F "\t" -v OFS="\t" '{{print $1,$2,$3,NR,".",$4}}' \
    > {output[4]}
    echo "{output[4]}:" >> {log} 2>&1
    cat {output[4]} | wc -l >> {log} 2>&1
    """
    
rule peak_annotation_homer:
  input:
    "06_macs2/{sample}_rep1_common_peaks.bed"
  output:
    "06_macs2/peak_annotation/{sample}_annotation_homer.xls"
  log:
    "logs/peak_annotation_homer/{sample}_annotation_homer.log"
  threads:5
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
    "06_macs2/peak_annotation/Lysate_annotation_homer.xls",
    "06_macs2/peak_annotation/Result_annotation_homer.xls"
  output:
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_pie.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_bar.pdf",
    "06_macs2/peak_annotation/Distribution_of_m6A_peaks_percent.txt"
  log:
    "logs/Distribution_of_m6A_peaks/Distribution_of_m6A_peaks.log" 
  threads:4
  script:
    "scripts/Distribution_of_m6A_peaks.R"
