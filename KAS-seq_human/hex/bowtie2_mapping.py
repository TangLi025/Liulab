TREATMENT=["input","IP"]
REP=["rep1","rep2"]
DUP=["raw","rmdup"]

SAMPLE=["T00","T05","T15","T30","T60"]


INDEX="~/reference/index/bowtie2/hg19/hg19"

rule all:
  input:
    expand("04_bam_{dup}/{sample}_{treatment}_{rep}.bam.bai",dup=DUP,sample=SAMPLE,treatment=TREATMENT,rep=REP)
    
rule bowtie2_mapping:
  input:
    "02_trim_galore/{sample}_{treatment}_{rep}_trimmed.fq"
  output:
    temp("04_bam_raw/{sample}_{treatment}_{rep}.sam"),
    "04_bam_raw/{sample}_{treatment}_{rep}.bam",
    "04_bam_raw/{sample}_{treatment}_{rep}.summary.txt"
  params:
    index=INDEX
  log:
    "logs/bowtie2_mapping/{sample}_{treatment}_{rep}.log"
  threads: 10
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/bowtie2/bin/bowtie2 -x {params.index} \
    -p {threads} -U {input} -S {output[0]} > {output[2]} 2>&1
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -q 20 -bS {output[0]} \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools sort -@ {threads} -o {output[1]} > {log} 2>&1
    """
  
rule samtools_rmdup:
  input:
    "04_bam_raw/{sample}_{treatment}_{rep}.bam"
  output:
    "04_bam_rmdup/{sample}_{treatment}_{rep}.bam"
  log:
    "logs/samtools_rmdup/{sample}_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools rmdup -s {input} {output} > {log} 2>&1"
    
rule bam_index:
  input:
    "04_bam_{dup}/{sample}_{treatment}_{rep}.bam"
  output:
    "04_bam_{dup}/{sample}_{treatment}_{rep}.bam.bai"
  log:
    "logs/bam_index/{dup}/{sample}_{treatment}_{rep}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"
    

