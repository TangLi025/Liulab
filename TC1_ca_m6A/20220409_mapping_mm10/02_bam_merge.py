SAMPLE=["p0","p5","p10","rp2"]
TREATMENT=["input","ip"]
REP=["rep1","rep2"]
READ=["1","2"]

DUP=["raw","dedup"]
STRAND=["pos","neg"]

rule all:
  input:
    expand("03_bam_merge/{dup}/{sample}_{treatment}_{rep}.bam",dup=DUP,sample=SAMPLE,treatment=TREATMENT,rep=REP)
    
rule bam_separate:
  input:
    "01_bam_{dup}/{sample}_{treatment}_{rep}.bam",
    "01_bam_{dup}/{sample}_{treatment}_{rep}.bam.bai"
  output:
    temp("02_bam_{dup}_separated/{sample}_{treatment}_{rep}_83.bam"),
    temp("02_bam_{dup}_separated/{sample}_{treatment}_{rep}_163.bam"),
    temp("02_bam_{dup}_separated/{sample}_{treatment}_{rep}_99.bam"),
    temp("02_bam_{dup}_separated/{sample}_{treatment}_{rep}_147.bam"),
    "02_bam_{dup}_separated/{sample}_{treatment}_{rep}_pos.bam",
    "02_bam_{dup}_separated/{sample}_{treatment}_{rep}_neg.bam"
  threads: 10
  shell:
    """
  	/disk/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -q 20 -b -f 83 {input[0]} 1> {output[0]}
  	/disk/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -q 20 -b -f 163 {input[0]} 1> {output[1]}
  	/disk/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -q 20 -b -f 99 {input[0]} 1> {output[2]}
  	/disk/user_09/anaconda3/envs/LinLong/bin/samtools view -@ {threads} -q 20 -b -f 147 {input[0]} 1> {output[3]}
  	/disk/user_09/anaconda3/envs/LinLong/bin/samtools merge -@ {threads} {output[4]} {output[0]} {output[1]} 
  	/disk/user_09/anaconda3/envs/LinLong/bin/samtools merge -@ {threads} {output[5]} {output[2]} {output[3]}
  	"""
    
rule bam_separated_index:
  input:
    "02_bam_{dup}_separated/{sample}_{treatment}_{rep}_{strand}.bam"
  output:
    "02_bam_{dup}_separated/{sample}_{treatment}_{rep}_{strand}.bam.bai"
  log:
    "logs/02_bam_{dup}_separated/bam_{strand}_index_{sample}_{treatment}_{rep}.log"
  threads: 2
  shell:
    "/disk/user_09/anaconda3/envs/LinLong/bin/samtools index -@ {threads} {input} > {log} 2>&1"

rule bam_merge:
  input:
    "02_bam_{dup}_separated/{sample}_{treatment}_{rep}_pos.bam",
    "02_bam_{dup}_separated/{sample}_{treatment}_{rep}_neg.bam"
  output:
    "03_bam_merge/{dup}/{sample}_{treatment}_{rep}.bam"
  log:
    "logs/bam_merge/{dup}/{sample}_{treatment}_{rep}.log"
  threads:10
  shell:
    "/disk/user_09/anaconda3/envs/LinLong/bin/samtools merge -f -@ {threads} {output} {input[0]} {input[1]} > {log} 2>&1"

