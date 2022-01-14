SAMPLE=["Lysate","Result"]
DUP=["raw","rmdup"]
TREATMENT=["input","IP"]
REP=["rep1","rep2"]
READ=["1","2"]
STRAND=["neg","pos"]


INDEX="/disk1/home/user_09/reference/index/hisat2/mm10/mm10"

rule all:
  input:
     expand("04_bam_{dup}/{sample}_{treatment}_{rep}.bam.bai",dup=DUP,sample=SAMPLE,treatment=TREATMENT,rep=REP),
     
rule hisat2_mapping:
  input:
    "02_fastp/{sample}_{treatment}_{rep}_1.fastq",
    "02_fastp/{sample}_{treatment}_{rep}_2.fastq"
  output:
    bam="04_bam_raw/{sample}_{treatment}_{rep}.bam",
    summary="04_bam_raw/{sample}_{treatment}_{rep}.summary.txt"
  params:
    index=INDEX
  log:
    "logs/hisat2_mapping/{sample}_{treatment}_{rep}.log"
  threads: 10
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/hisat2 -x {params.index} --rna-strandness RF \
    --summary-file {output.summary} \
    -p {threads} -1 {input[0]} -2 {input[1]} \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -q 20 -bS \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools sort -@ {threads} -o {output.bam} 
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
    "logs/bam_index/{sample}_{treatment}_{rep}_{dup}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"
    
