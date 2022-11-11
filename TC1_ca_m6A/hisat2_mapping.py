SAMPLE=["TC"]
REP=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"]
READ=["1","2"]

DUP=["raw","dedup"]
STRAND=["pos","neg"]

INDEX="/disk1/home/user_09/reference/index/hisat2/grcm38/genome"

rule all:
  input:
    expand("04_bam_{dup}/{sample}{rep}.bam.bai",dup=DUP,sample=SAMPLE,rep=REP)

rule hisat2_mapping:
  input:
    "02_trim_galore/{sample}{rep}_1_val_1.fq.gz",
    "02_trim_galore/{sample}{rep}_2_val_2.fq.gz"
  output:
    bam="04_bam_raw/{sample}{rep}.bam",
    summary="04_bam_raw/{sample}{rep}.summary.txt"
  params:
    index=INDEX
  log:
    "logs/hisat2_mapping/{sample}{rep}.log"
  threads: 10
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/hisat2 -x {params.index} --rna-strandness RF \
    --summary-file {output.summary} \
    -p {threads} -1 {input[0]} -2 {input[1]} \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -bS \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools sort -@ {threads} -o {output.bam} 
    """
    
rule deduplicate:
  input:
    "04_bam_raw/{sample}{rep}.bam"
  output:
    "04_bam_dedup/dedup_record/{sample}{rep}.log",
    "04_bam_dedup/{sample}{rep}.bam"
  log:
    "logs/bam_dedup/deduplicate_{sample}{rep}.log"
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/picard MarkDuplicates REMOVE_DUPLICATES=true \
      I={input} M={output[0]} O={output[1]} > {log} 2>&1"

rule bam_index:
  input:
    "04_bam_{dup}/{sample}{rep}.bam"
  output:
    "04_bam_{dup}/{sample}{rep}.bam.bai"
  log:
    "logs/bam_index/{sample}{rep}_{dup}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"
