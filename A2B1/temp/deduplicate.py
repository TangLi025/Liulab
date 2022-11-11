SAMPLE=["Lysate","Result"]
TREATMENT=["input","IP"]
REP=["rep1","rep2"]
READ=["1","2"]

DUP=["raw","dedup"]
STRAND=["pos","neg"]

INDEX="/disk1/home/user_09/reference/index/hisat2/mm19_hisat2_ss_exon/mm19"

rule all:
  input:
    #expand("02_trim_galore/{sample}_{treatment}_{rep}_1_val_1.fq.gz",sample=SAMPLE,treatment=TREATMENT,rep=REP),
    #expand("02_trim_galore/{sample}_{treatment}_{rep}_2_val_2.fq.gz",sample=SAMPLE,treatment=TREATMENT,rep=REP),
    expand("04_bam_{dup}/{sample}.{rep}.{treatment}.bam.bai",dup=DUP,sample=SAMPLE,treatment=TREATMENT,rep=REP),
    #"01_fastqc1/multiqc_report.html",
    #"03_fastqc2/multiqc_report.html"

rule deduplicate:
  input:
    "04_bam_raw/{sample}.{rep}.{treatment}.bam"
  output:
    "04_bam_dedup/dedup_record/{sample}.{rep}.{treatment}.log",
    "04_bam_dedup/{sample}.{rep}.{treatment}.bam"
  log:
    "logs/bam_dedup/deduplicate_{sample}.{rep}.{treatment}.log"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/picard MarkDuplicates REMOVE_DUPLICATES=true \
      I={input} M={output[0]} O={output[1]} > {log} 2>&1"

rule bam_index:
  input:
    "04_bam_{dup}/{sample}.{rep}.{treatment}.bam"
  output:
    "04_bam_{dup}/{sample}.{rep}.{treatment}.bam.bai"
  log:
    "logs/bam_index/{sample}.{rep}.{treatment}_{dup}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"
    
