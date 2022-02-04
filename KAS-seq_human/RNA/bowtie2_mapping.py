TREATMENT=["input","IP"]
REP=["rep1","rep2"]

INDEX="~/reference/index/hisat2/hg19/hg19"

rule all:
  input:
    expand("04_bam_raw/HEK_{rep}.bam.bai",rep=REP)
    
rule hisat2_mapping:
  input:
    "02_trim_galore/HEK_{rep}_trimmed.fq"
  output:
    temp("04_bam_raw/HEK_{rep}.sam"),
    bam="04_bam_raw/HEK_{rep}.bam",
    summary="04_bam_raw/HEK_{rep}.summary.txt"
  params:
    index=INDEX
  log:
    "logs/hisat2_mapping/HEK_{rep}.log"
  threads: 10
  shell:
    """
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/hisat2 -x {params.index} \
    --summary-file {output.summary} \
    -p {threads} -U {input} 1> {output[0]} 2> {log}
    /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -q 20 -bS {output[0]} \
    | /disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools sort -@ {threads} -o {output.bam} 
    """
    
rule bam_index:
  input:
    "04_bam_raw/HEK_{rep}.bam"
  output:
    "04_bam_raw/HEK_{rep}.bam.bai"
  log:
    "logs/bam_index/HEK_{rep}.log"
  threads: 4
  shell:
    "/disk1/home/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"
    

