SAMPLE=["Lysate","Result"]
TREATMENT=["input","IP"]
REP=["rep1","rep2"]
READ=["1","2"]

DUP=["raw","dedup"]
STRAND=["pos","neg"]

GENOME="/disk1/home/user_09/reference/genome/GRCm39.genome.fa"
GTF="/disk1/home/user_09/reference/annotation/mm19/gencode.vM28.annotation.gtf"

rule all:
  input:
    expand("10_bed_merge/{dup}/Distribution_of_m6A_peaks_percent.txt",dup=DUP),
    expand("10_bed_merge/{dup}/annotation_output/{sample}_annotation_homer_rearranged.xls",sample=SAMPLE,dup=DUP),
    expand("10_bed_merge/{dup}/Metagene_profiles_of_m6A_peak_density.png",dup=DUP),
    

rule Metagene_profile_of_m6A_peak_density:
  input:
    "10_bed_merge/{dup}/Lysate_common_peaks.bed",
    "10_bed_merge/{dup}/Result_common_peaks.bed"
  output:
    "10_bed_merge/{dup}/Metagene_profiles_of_m6A_peak_density.png"
  log:
    "logs/10_bed_merge/{dup}/Metagene_profiles_of_m6A_peak_density.log"
  threads:1
  script:
    "scripts/Metagene-profiles-of-m6A-peak-density.R"

rule reproducibility:
  input:
    "10_bed_merge/{dup}/{sample}_common_peaks.bed",
    "11_bam_merge/{dup}/{sample}_input_rep1.bam",
    "11_bam_merge/{dup}/{sample}_IP_rep1.bam",
    "11_bam_merge/{dup}/{sample}_input_rep2.bam",
    "11_bam_merge/{dup}/{sample}_IP_rep2.bam"
  output:
    "11_bam_merge/{dup}/reproducibility/{sample}_reproducibility.pdf",
    "10_bed_merge/{dup}/{sample}_zero_deprived.bed"
  log:
    "logs/11_bam_merge/{dup}/reproducibility/{sample}_reproducibility.log"
  threads:10
  script:
    "scripts/reproducibility.R"
    
rule peak_annotation_homer:
  input:
    "10_bed_merge/{dup}/{sample}_zero_deprived.bed"
  output:
    "10_bed_merge/{dup}/{sample}_zero_deprived_rearranged.bed",
    "10_bed_merge/{dup}/annotation_output/{sample}_annotation_homer_rearranged.xls"
  log:
    "logs/10_bed_merge/{dup}/annotation_output/{sample}_annotation_homer_rearranged.log"
  params:
    genome=GENOME,
    gtf=GTF
  threads:10
  shell:
    """
    cat {input[0]} \
      | awk -F "\t" -v OFS="\t" '{{$5=".";print $4,$1,$2,$3,$6,$5}}' \
      > {output[0]}
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/annotatePeaks.pl {output[0]} \
      {params.genome} \
      -gtf {params.gtf} \
      -cpu {threads} \
      1> {output[1]} 2> {log}
    """
    
rule Distribution_of_m6A_peaks:
  input:
    "10_bed_merge/{dup}/annotation_output/Lysate_annotation_homer_rearranged.xls",
    "10_bed_merge/{dup}/annotation_output/Result_annotation_homer_rearranged.xls"
  output:
    "10_bed_merge/{dup}/Distribution_of_m6A_peaks_pie.pdf",
    "10_bed_merge/{dup}/Distribution_of_m6A_peaks_bar.pdf",
    "10_bed_merge/{dup}/Distribution_of_m6A_peaks_percent.txt"
  log:
    "logs/Distribution_of_m6A_peaks/{dup}/Distribution_of_m6A_peaks.log" 
  threads:4
  script:
    "scripts/Distribution_of_m6A_peaks.R"
