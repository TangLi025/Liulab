GROUP=["METTL3_2"]
SAMPLE=["CTRL","KO"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

DUP=["raw","rmdup"]

INDEX="~/reference/index/bowtie2/mm19/mm19"

rule all:
  input:
    expand("04_bam_{dup}/KAS-seq_{group}_{sample}_{treatment}_{rep}.bam.bai",dup=DUP,group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP)
    
rule bowtie2_mapping:
  input:
    "02_trim_galore/KAS-seq_{group}_{sample}_{treatment}_{rep}_trimmed.fq.gz"
  output:
    temp("04_bam_raw/KAS-seq_{group}_{sample}_{treatment}_{rep}.sam"),
    "04_bam_raw/KAS-seq_{group}_{sample}_{treatment}_{rep}.bam",
    "04_bam_raw/KAS-seq_{group}_{sample}_{treatment}_{rep}.summary.txt"
  params:
    index=INDEX
  log:
    "logs/bowtie2_mapping/KAS-seq_{group}_{sample}_{treatment}_{rep}.log"
  threads: 10
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/bowtie2/bin/bowtie2 -x {params.index} \
      -p {threads} -U {input} -S {output[0]} > {output[2]} 2>&1
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools sort -@ {threads} {output[0]} -o {output[1]} > {log} 2>&1
    """
  
rule samtools_rmdup:
  input:
    "04_bam_raw/KAS-seq_{group}_{sample}_{treatment}_{rep}.bam"
  output:
    "04_bam_rmdup/KAS-seq_{group}_{sample}_{treatment}_{rep}.bam"
  log:
    "logs/samtools_rmdup/KAS-seq_{group}_{sample}_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools rmdup -s {input} {output} > {log} 2>&1"
    
rule bam_index:
  input:
    "04_bam_{dup}/KAS-seq_{group}_{sample}_{treatment}_{rep}.bam"
  output:
    "04_bam_{dup}/KAS-seq_{group}_{sample}_{treatment}_{rep}.bam.bai"
  log:
    "logs/bam_index/{dup}/KAS-seq_{group}_{sample}_{treatment}_{rep}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"
    

