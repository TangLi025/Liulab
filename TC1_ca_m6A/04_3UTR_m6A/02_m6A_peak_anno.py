SAMPLE=["p0p5","p0p10","rp2p5","rp2p10"]

GENOME="/disk/user_09/reference/genome/mm/GRCm39.genome.fa"
GTF="/disk/user_09/reference/annotation/mm19/gencode.vM28.annotation.gtf"



rule all:
  input:
    expand("10_bed_merge/peak_annotation_homer/{sample}_annotation_homer_rearranged.xls",sample=SAMPLE),
    expand("10_bed_merge/mm10/peak_annotation_homer/{sample}_annotation_homer_rearranged.xls",sample=SAMPLE)

rule peak_annotation_homer:
  input:
    "/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/10_bed_merge/02_merge_peaks/{sample}.merge.bed"
  output:
    "10_bed_merge/{sample}_zero_deprived_rearranged.bed",
    "10_bed_merge/peak_annotation_homer/{sample}_annotation_homer_rearranged.xls"
  log:
    "logs/10_bed_merge/peak_annotation_homer/{sample}_annotation_homer_rearranged.log"
  params:
    genome=GENOME,
    gtf=GTF
  threads:10
  shell:
    """
    cat {input[0]} \
      | awk -F "\t" -v OFS="\t" '{{$5=".";print $4,$1,$2,$3,$6,$5}}' \
      > {output[0]}
    /disk/user_09/anaconda3/envs/LinLong/bin/annotatePeaks.pl {output[0]} \
      {params.genome} \
      -gtf {params.gtf} \
      -cpu {threads} \
      1> {output[1]} 2> {log}
    """

rule peak_annotation_homer_mm10:
  input:
    "/disk/user_09/Data/user_08_TC1/05_bam_change_index_hisat2/01_bam_sorted/10_bed_merge/02_merge_peaks/mm39ToMm10/{sample}.merge.bed"
  output:
    "10_bed_merge/mm10/{sample}_zero_deprived_rearranged.bed",
    "10_bed_merge/mm10/peak_annotation_homer/{sample}_annotation_homer_rearranged.xls"
  log:
    "logs/10_bed_merge/peak_annotation_homer_mm10/{sample}_annotation_homer_rearranged.log"
  params:
    genome=GENOME,
    gtf=GTF
  threads:10
  shell:
    """
    cat {input[0]} \
      | awk -F "\t" -v OFS="\t" '{{$5=".";print $4,$1,$2,$3,$6,$5}}' \
      > {output[0]}
    /disk/user_09/anaconda3/envs/LinLong/bin/annotatePeaks.pl {output[0]} \
      mm10 \
      -cpu {threads} \
      1> {output[1]} 2> {log}
    """