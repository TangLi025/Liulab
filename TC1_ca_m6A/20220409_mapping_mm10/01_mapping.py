SAMPLE=["p0","p5","p10","rp2"]
TREATMENT=["input","ip"]
REP=["rep1","rep2"]
READ=["1","2"]

DUP=["raw","dedup"]
STRAND=["pos","neg"]

INDEX="/disk/user_09/reference/index/hisat2/mm10/mm10"

rule all:
  input:
    expand("01_bam_{dup}/{sample}_{treatment}_{rep}.bam.bai",dup=DUP,sample=SAMPLE,treatment=TREATMENT,rep=REP)
    
rule hisat2_mapping:
  input:
    "/disk/user_08/Data/TC1-planB/04_bam_rm_rRNA_bowtie2/{sample}_{treatment}_{rep}.1.derRNA.fq.gz",
    "/disk/user_08/Data/TC1-planB/04_bam_rm_rRNA_bowtie2/{sample}_{treatment}_{rep}.2.derRNA.fq.gz"
  output:
    bam="01_bam_raw/{sample}_{treatment}_{rep}.bam",
    summary="01_bam_raw/{sample}_{treatment}_{rep}.summary.txt"
  params:
    index=INDEX
  log:
    "logs/hisat2_mapping/{sample}_{treatment}_{rep}.log"
  threads: 25
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/hisat2 -x {params.index} --rna-strandness RF \
    --summary-file {output.summary} \
    -p {threads} -1 {input[0]} -2 {input[1]} \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -bS \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools sort -@ {threads} -o {output.bam} \
    1> {log} 2>&1
    """
    
rule deduplicate:
  input:
    "01_bam_raw/{sample}_{treatment}_{rep}.bam"
  output:
    "01_bam_dedup/dedup_record/{sample}_{treatment}_{rep}.log",
    "01_bam_dedup/{sample}_{treatment}_{rep}.bam"
  log:
    "logs/bam_dedup/deduplicate_{sample}_{treatment}_{rep}.log"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/picard MarkDuplicates REMOVE_DUPLICATES=true \
      I={input} M={output[0]} O={output[1]} 1> {log} 2>&1"

rule bam_index:
  input:
    "01_bam_{dup}/{sample}_{treatment}_{rep}.bam"
  output:
    "01_bam_{dup}/{sample}_{treatment}_{rep}.bam.bai"
  log:
    "logs/bam_index/{sample}_{treatment}_{rep}_{dup}.log"
  threads: 2
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} 1> {log} 2>&1"
    
