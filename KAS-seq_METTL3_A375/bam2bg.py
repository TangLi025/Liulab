GROUP=["A375"]
SAMPLE=["Control","KDMETTL3"]

TREATMENT=["input","IP"]
REP=["rep1","rep2"]

GENOME="/disk1/home/user_09/reference/genome/GRCh38.p13.genome.fa"
GENOME_INDEX="/disk1/home/user_09/reference/genome/GRCh38.p13.genome.fa.fai"

rule all:
  input:
    expand("{group}/05_bedtools/bedGraph/{group}_{sample}_{treatment}_{rep}_ext.bg",group=GROUP,sample=SAMPLE,treatment=TREATMENT,rep=REP)
  
rule genome_index:
  input:
    GENOME
  output:
    GENOME_INDEX
  shell:
    "/disk1/home/user_09/anaconda3/envs/m6A/bin/samtools faidx {input}"

rule bam2bed:
  input:
    "{group}/04_bam_rmdup/{group}_{sample}_{treatment}_{rep}.bam"
  output:
    "{group}/05_bedtools/bed/{group}_{sample}_{treatment}_{rep}.bed"
  log:
    "{group}/logs/bam2bed/{group}_{sample}_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools bamtobed -i {input} \
      | sort -k 1,1 1> {output} 2> {log}"
      
rule bed_extend:
  input:
    "{group}/05_bedtools/bed/{group}_{sample}_{treatment}_{rep}.bed"
  output:
    "{group}/05_bedtools/bed_extend/{group}_{sample}_{treatment}_{rep}_ext.bed"
  log:
    "{group}/logs/bed_extend/{group}_{sample}_{treatment}_{rep}.log"
  shell:
    r"""
    awk '$3-150>0 {{if ($6~"+") printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$2,$2+150,$4,$5,$6); else if ($6~"-") printf("%s\t%d\t%d\t%s\t%d\t%s\n",$1,$3-150,$3,$4,$5,$6)}}' \
      {input} 1> {output} 2> {log}
    """

rule bed2bedGraph:
  input:
    "{group}/05_bedtools/bed_extend/{group}_{sample}_{treatment}_{rep}_ext.bed",
    GENOME_INDEX
  output:
    "{group}/05_bedtools/bedGraph/{group}_{sample}_{treatment}_{rep}_ext.bg"
  log:
    "{group}/logs/bed2bedGraph/{group}_{sample}_{treatment}_{rep}.log"
  threads: 1
  shell:
    "/disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools genomecov -i {input[0]} \
      -bg -g {input[1]} \
      | /disk1/home/user_09/anaconda3/envs/bedtools/bin/bedtools sort -i 1> {output} 2> {log}"
