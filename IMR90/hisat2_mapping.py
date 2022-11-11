SAMPLE=["CTRL","SENE"]

REP=["rep1","rep2"]

DUP=["raw"]


INDEX="/disk1/home/user_09/reference/index/hisat2/hg38/hg38"

rule all:
  input:
    expand("04_bam_{dup}/{sample}_{rep}.bam.bai",dup=DUP,sample=SAMPLE,rep=REP)
    
rule hisat2_mapping:
  input:
    "02_trim_galore/{sample}_{rep}_trimmed.fq"
  output:
    temp("04_bam_raw/{sample}_{rep}.sam"),
    "04_bam_raw/{sample}_{rep}.bam",
    summary="04_bam_raw/{sample}_{rep}.summary.txt"
  params:
    index=INDEX
  log:
    "logs/bowtie2_mapping/{sample}_{rep}.log"
  threads: 10
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/m6A/bin/hisat2 -x {params.index} \
      --summary-file {output.summary} \
      -p {threads} -U {input} -S {output[0]} > {log} 2>&1
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -bS -q 20  {output[0]} | 
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools sort -@ {threads} -o {output[1]} >> {log} 2>&1
    """
    
rule bam_index:
  input:
    "04_bam_{dup}/{sample}_{rep}.bam"
  output:
    "04_bam_{dup}/{sample}_{rep}.bam.bai"
  log:
    "logs/bam_index/{dup}/{sample}_{rep}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"
    

