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
    expand("10_bed_merge/{dup}/Lysate_Result_merge_peaks.bed",dup=DUP),
    expand("10_bed_merge/{dup}/annotation_output/Lysate_Result_merge_peaks_annotation_homer_rearranged.xls",dup=DUP),
    
    
rule bed_merge:
  input:
    "10_bed_merge/{dup}/Lysate_common_peaks.bed",
    "10_bed_merge/{dup}/Result_common_peaks.bed"
  output:
    "10_bed_merge/{dup}/Lysate_Result_merge_peaks.bed"
  log:
    "logs/10_bed_merge/{dup}/merge.log"
  shell:
    """
    cat {input[0]} {input[1]} \
        | sort -k1,1 -k2,2n \
        | /disk1/home/user_09/anaconda3/envs/LinLong/bin/bedtools merge -s -c 6 -o distinct \
        | awk -F "\t" -v OFS="\t" '{{print $1,$2,$3,NR,".",$4}}' \
        > {output}
        echo "{output}:" >> {log} 2>&1
        cat {output} | wc -l >> {log} 2>&1
    """

rule peak_annotation_homer:
  input:
    "10_bed_merge/{dup}/Lysate_Result_merge_peaks.bed"
  output:
    "10_bed_merge/{dup}/Lysate_Result_merge_peaks_rearanged.bed",
    "10_bed_merge/{dup}/annotation_output/Lysate_Result_merge_peaks_annotation_homer_rearranged.xls"
  log:
    "logs/10_bed_merge/{dup}/annotation_output/Lysate_Result_merge_peaks_annotation_homer_rearranged.log"
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
    
