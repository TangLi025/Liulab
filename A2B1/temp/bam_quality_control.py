SAMPLE=["Lysate","Result"]
TREATMENT=["input","IP"]
REP=["rep1","rep2"]
READ=["1","2"]

DUP=["raw","rmdup"]
STRAND=["pos","neg"]

INDEX="/disk1/home/user_09/reference/index/hisat2/mm19_hisat2_ss_exon/mm19"

rule all:
  input:
    expand("04_bam_unique/{sample}_{treatment}_{rep}.bam.bai",sample=SAMPLE,treatment=TREATMENT,rep=REP),


rule unique:
  input:
    "04_bam_raw/{sample}_{treatment}_{rep}.bam"
  output:
    "04_bam_unique/{sample}_{treatment}_{rep}.bam",
  log:
    "logs/unique/{sample}_{treatment}_{rep}.log"
  threads: 10
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -q 20 -bS {input} \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools sort -@ {threads} -o {output} > {log} 2>&1
    """

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
    
